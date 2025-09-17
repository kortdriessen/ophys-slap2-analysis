import sys
from pathlib import Path
import numpy as np
import torch
import scipy.io as spio
import scipy.ndimage as ndimage
from scipy.interpolate import RectBivariateSpline, PchipInterpolator
import os
import time
import multiprocessing as mp
from functools import partial
import tkinter as tk
from tkinter import filedialog, ttk, messagebox
import datetime
import importlib
import skimage.io as skimio
import h5py
from scipy import stats, signal
import slap2_utils
import pandas as pd

from tqdm import tqdm 

sys.path.append(str(Path(__file__).parent.parent))
import reconstruct

import cv2

def nearest_interp(x, xp, yp):
    if len(xp) == 1:
        return yp
    
    x_bds = xp[:-1] / 2.0 + xp[1:] / 2.0
    idx = np.searchsorted(x_bds, x, side='left')
    idx = np.clip(idx, 0, len(xp) - 1)

    return yp[idx]

def fast_dilation(mask, kernel, iterations=1):
    if kernel is None:
        kernel = np.ones((7, 7), np.uint8)  # Structuring element for XY dilation
    out = np.empty_like(mask)
    
    # Iterate over Z slices (first dimension)
    for z in range(mask.shape[0]):
        out[z] = cv2.dilate(mask[z].astype(np.uint8), kernel, iterations=iterations)
    
    return out.astype(bool)

def _movmean_nan(x: np.ndarray, window: int) -> np.ndarray:
    """
    Centered moving mean that ignores NaNs, similar to MATLAB smoothdata(...,'movmean','omitnan').
    Uses partial windows at the edges (like MATLAB), returns NaN where the window has no valid data.
    """
    window = max(int(window), 1)
    if window == 1:
        return x.astype(float, copy=True)

    kernel = np.ones(window, dtype=float)
    valid = ~np.isnan(x)
    x_filled = np.where(valid, x, 0.0)

    num = np.convolve(x_filled, kernel, mode="same")
    den = np.convolve(valid.astype(float), kernel, mode="same")

    out = num / np.where(den == 0.0, np.nan, den)
    return out

def compute_f0(Fin, denoise_window: int, hull_window: int):
    """
    Python port of:

      function [F0] = computeF0(Fin, denoiseWindow, hullWindow)
      % 1) median filter to reduce noise
      % 2) rolling convex-hull-like min envelope on decimated grids
      % 3) discard doubtful samples near NaNs; smooth/fill; PCHIP back to full grid

    Args
    ----
    Fin : array_like, shape (T, ...) with time along axis 0 (NaNs allowed)
    denoise_window : int
        Median filter window (time samples).
    hull_window : int
        Window that controls the “convex hull”-like operation.

    Returns
    -------
    F0 : ndarray, same shape as Fin
    """
    F = np.asarray(Fin)
    orig_shape = F.shape
    if F.ndim == 1:
        F = F[:, None]
    else:
        F = F.reshape(F.shape[0], -1)
    
    T, C = F.shape
    if T < 4:
        return np.ones_like(Fin, dtype=float) * np.nanmean(Fin, axis=0, keepdims=True)

    hull_window = int(min(hull_window, T//4))
    delta_des = max(4.0, denoise_window / 6.0)
    
    sample_times = np.rint(np.linspace(0, T-1, num=int(np.ceil(T/delta_des)+1))).astype(int)
    n_samps_in_hull = int(np.ceil(hull_window / delta_des))

    F0 = pd.DataFrame(F).rolling(window=denoise_window, center=True, min_periods=1).median().to_numpy()
    
    origsz = F0.shape
    F0 = F0.reshape(origsz[0], -1)
    for cix in range(F0.shape[1]):
        if np.all(np.isnan(F0[:, cix])):
            continue
        
        F00 = np.full((sample_times.shape[0], n_samps_in_hull), np.nan)
        for dix in range(n_samps_in_hull, 0, -1):
            xi = sample_times[dix-1::n_samps_in_hull]
            F00[:,dix-1] = np.interp(sample_times,xi, F0[xi,cix], left=np.nan, right=np.nan)
        FF = np.nanmin(F00,axis=1)

        doubt = np.sum(~np.isnan(F00),axis=1) < int(np.ceil(n_samps_in_hull/2))
        if np.sum(~doubt)>2:
            FF[doubt] = np.nan

        win = 2*int(np.ceil(n_samps_in_hull/2.0)) + 1
        fill = _movmean_nan(FF, win)
        nan_mask = np.isnan(FF)
        FF[nan_mask] = fill[nan_mask]
        FF = _movmean_nan(FF, win)

        pchip = PchipInterpolator(sample_times, FF, extrapolate=True)
        F0[:,cix] = pchip(np.arange(T))

    return F0.reshape(orig_shape)

def get_trial_data(trial_info, DMDix, params, sampFreq, refStack, fastZ2RefZ, allSuperPixelIDs, dr, trialTable):
    trialIx, keepTrial = trial_info
    
    if not keepTrial:
        return [], [], []

    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
    numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]

    nPixels = dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs

    source_fn = trialTable['filename'][DMDix,trialIx][0]
    firstLine = trialTable['firstLine'][DMDix,trialIx]
    lastLine = trialTable['lastLine'][DMDix,trialIx]

    importlib.reload(slap2_utils)
    hDataFile = slap2_utils.DataFile(os.path.join(dr, source_fn))

    linesPerCycle = hDataFile.header['linesPerCycle']

    dt = 1/sampFreq/hDataFile.metaData.linePeriod_s

    DSframes = np.ceil(np.arange(firstLine, lastLine+1, dt))
    nDSframes= len(DSframes)

    # Pre-compute time windows for all frames
    dtRead = max(3 * dt, linesPerCycle)
    timeWindows = [np.arange(max(1,np.floor(DSframes[i]-dtRead)), 
                            min(np.ceil(DSframes[i]+dtRead),hDataFile.numCycles*linesPerCycle)+1) 
                    for i in range(nDSframes)]

    # Pre-compute line and cycle indices for all time windows
    lineIndices_all = [(tw - 1) % linesPerCycle + 1 for tw in timeWindows]
    cycleIndices_all = [np.floor((tw - 1) / linesPerCycle) + 1 for tw in timeWindows]

    # Collect all unique line-cycle combinations across all frames
    all_line_cycles = set()
    for DSframeIx in range(nDSframes):
        for li, ci in zip(lineIndices_all[DSframeIx], cycleIndices_all[DSframeIx]):
            all_line_cycles.add((int(li), int(ci)))

    # Get all line data at once
    all_lines = np.array([lc[0] for lc in all_line_cycles])
    all_cycles = np.array([lc[1] for lc in all_line_cycles])

    # Create a mapping from (line, cycle) to index in the cache
    line_cycle_to_idx = {(int(li), int(ci)): i for i, (li, ci) in enumerate(zip(all_lines, all_cycles))}

    start_line_data_time = time.time()
    print(f"Getting line data for {len(all_lines)} lines", end="")
    # Get all line data at once
    all_line_data = hDataFile.getLineData(all_lines, all_cycles, params['activityChannel'])
    print(f" - completed in {time.time() - start_line_data_time:.3f} sec")

    data = np.zeros((numSuperPixels,nDSframes))
    dataCt = np.zeros((numSuperPixels,nDSframes))
    # Initialize timing variables
    start_time = time.time()

    # Process each frame
    for DSframeIx in range(nDSframes):
        if DSframeIx % 100 == 0:
            avg_time = (time.time() - start_time) / DSframeIx if DSframeIx > 0 else 0
            print(f"{DSframeIx} of {nDSframes}, Average time per frame: {avg_time:.3f} sec")
        
        weights = np.exp(-np.abs(DSframes[DSframeIx] - timeWindows[DSframeIx]) / dt)

        # Get the line and cycle indices for this frame
        frame_line_indices = lineIndices_all[DSframeIx]
        frame_cycle_indices = cycleIndices_all[DSframeIx]
        
        # Process only valid lines
        valid_lines = [i for i, li in enumerate(frame_line_indices) 
                      if hDataFile.lineDataNumElements[int(li)-1] != 0]
        
        for i in valid_lines:
            line_idx = int(frame_line_indices[i])
            cycle_idx = int(frame_cycle_indices[i])
            
            # Get the cached line data
            cache_idx = line_cycle_to_idx[(line_idx, cycle_idx)]
            line_data = all_line_data[cache_idx]
            
            # Get positions and z-index
            positions = hDataFile.lineSuperPixelIDs[line_idx-1]
            zIdx = hDataFile.lineFastZIdxs[line_idx-1]
            
            # Compute lookup values and matches
            lookup_values = positions * 100 + zIdx
            matching_mask = np.isin(allSuperPixelIDs[f'DMD{DMDix+1}'], lookup_values)
            matching_indices = np.where(matching_mask)[0]
            
            if len(matching_indices) > 0:
                # Create a mapping from lookup values to their indices
                value_to_pos = dict(zip(lookup_values.astype(np.uint32), range(len(lookup_values))))
                # Get the positions in lineData for each matching index
                matched_positions = [value_to_pos[int(allSuperPixelIDs[f'DMD{DMDix+1}'][idx])] for idx in matching_indices]
                
                weight = weights[i]
                data[matching_indices, DSframeIx] += line_data[matched_positions, 0] * weight
                dataCt[matching_indices, DSframeIx] += weight

    aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]
    aData['DSframes'] = aData['DSframes'] * hDataFile.metaData.linePeriod_s
    
    return data/100, dataCt, aData, DSframes * hDataFile.metaData.linePeriod_s

def get_high_res_traces(trial_info, DMDix, params, sampFreq, refStack, subsampleMatrixInds, fastZ2RefZ, sparseHInds, sparseHVals, 
                allSuperPixelIDs, dr, trialTable, A_final, uniqueMotionDS, motIndsToKeepDS, psf):
    trialIx, keepTrial, backgroundDS = trial_info
    
    if not isinstance(A_final, torch.Tensor):
        A_final = torch.tensor(A_final, dtype=torch.float32)

    nSources = A_final.shape[1]

    if not keepTrial:
        return np.full((0,nSources),np.nan,dtype=np.float32), \
            np.full((0,nSources),np.nan,dtype=np.float32), \
                np.full((0,),np.nan), \
                    np.full((0,),np.nan), \
                        np.full((0,),np.nan), \
                            (np.full((0,),np.nan), np.full((0,),np.nan), np.full((0,),np.nan)), \
                                (np.full((0,),0,dtype=np.int16), np.full((0,),0,dtype=np.int16), np.full((0,),0,dtype=np.int16))

    data_file = os.path.join(dr, f'trial_data_DMD{DMDix+1}_trial{trialIx}.npz')
    if os.path.exists(data_file):
        print(f'Loading existing trial data from {data_file}')
        data_arrays = np.load(data_file)
        dataNonNorm = data_arrays['dataNonNorm']
        dataCt = data_arrays['dataCt']
        frames = data_arrays['DSframes']

        aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]
        aData['DSframes'] = data_arrays['aData_DSframes']
    else:
        dataNonNorm, dataCt, aData, frames = get_trial_data(trial_info[:2], DMDix, params, sampFreq, refStack, fastZ2RefZ, allSuperPixelIDs, dr, trialTable)
        # Save data arrays
        data_arrays = {
            'dataNonNorm': dataNonNorm,
            'dataCt': dataCt,
            'DSframes': frames,
            'aData_DSframes': aData['DSframes']
        }
        np.savez(data_file, **data_arrays)
        print(f'Saved trial data to {data_file}')
    
    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
    numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]

    data = dataNonNorm / dataCt

    motionR = np.interp(frames, aData['DSframes'][0], aData['motionDSr'].T[0])
    motionC = np.interp(frames, aData['DSframes'][0], aData['motionDSc'].T[0])
    motionZ = np.interp(frames, aData['DSframes'][0], aData['motionDSz'].T[0])
    background = np.array([np.interp(frames, aData['DSframes'][0], backgroundDS[i]) for i in range(backgroundDS.shape[0])])

    onlineYshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineYshift'].T[0])
    onlineXshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineXshift'].T[0])
    onlineZshifts = nearest_interp(frames, aData['DSframes'][0], aData['onlineZshift'].T[0])

    uniqueMotion, motInds = np.unique(np.round(np.stack((motionR, motionC, motionZ), axis=1)), axis=0, return_inverse=True)
    
    motIndsToKeep = np.zeros_like(motIndsToKeepDS, dtype=np.int64)
    motIndsToKeep[:] = -1

    for i, motion_idx_DS in enumerate(motIndsToKeepDS):
        matches = np.all(uniqueMotion == uniqueMotionDS[motion_idx_DS,:],axis=1).nonzero()[0]
        if len(matches) > 0:
            motIndsToKeep[i] = matches[0]
    
    # background_spatial_components = background_spatial_components[:, motIndsToKeep != -1]
    motIndsToKeep = motIndsToKeep[motIndsToKeep != -1]

    uniqueMotionYX = np.unique(uniqueMotion[motIndsToKeep,:2],axis=0)

    refPixs = torch.from_numpy(subsampleMatrixInds[:,0])
    refD = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor').int()
    refC = torch.div((refPixs - refD * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor').int()
    refR = (refPixs % dmdPixelsPerColumn).int()

    selPixMask = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), dtype=bool)
    for i in range(uniqueMotionYX.shape[0]):
        selPixMask[refD,refR + int(uniqueMotionYX[i,0]),refC + int(uniqueMotionYX[i,1])] = True
    selPixMask = ndimage.binary_dilation(selPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
    selPixIdxs = np.flatnonzero(selPixMask)

    if not isinstance(A_final, torch.Tensor):
        A_final = torch.tensor(A_final, dtype=torch.float32)

    nSources = A_final.shape[1]
    phi = torch.full((data.shape[1], nSources), np.nan, dtype=torch.float32)
    F0 = torch.full((data.shape[1], nSources), np.nan, dtype=torch.float32)

    residual = data - background
    
    for i, motion_idx in enumerate(motIndsToKeep):

        # Extract data for the most common motion mode
        motion_frames = (motInds == motion_idx).nonzero()[0]
        
        # data_tensor = torch.from_numpy(data[:, motion_frames].astype(np.float32))
        data_tensor = torch.from_numpy(residual[:, motion_frames].astype(np.float32))

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)
        sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

        # project image space (A) into superpixel space (X)
        X = torch.sparse.mm(H, A_final[selPixIdxs,:])

        # add background spatial component
        # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

        XtX = X.T @ X
        Xtd = X.T @ data_tensor  # This gives all time points at once
        
        # Add small regularization to ensure stability
        regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
        
        # Solve the system for all time points at once
        # We need to solve (X^T * X) * phi = X^T * data for each column of data
        phi[motion_frames,:] = torch.linalg.solve(
            regularized_XtX,
            Xtd
        ).T

        Xtbackground = X.T @ torch.from_numpy(background[:,motion_frames].astype(np.float32))
        F0[motion_frames,:] = torch.linalg.solve(
            XtX + 1e-10 * torch.eye(XtX.shape[0]),
            Xtbackground
        ).T

        # Xtbackground = X[:,:-1].T @ (X[:,-1].unsqueeze(-1) * phi[motion_frames,-1].unsqueeze(-1).T)
        # F0[motion_frames,:] = torch.linalg.solve(
        #     XtX[:-1,:-1] + 1e-10 * torch.eye(XtX.shape[0]-1),
        #     Xtbackground
        # ).T

    globalF = np.sum(data, axis=0)

    return phi.numpy(), F0.numpy(), frames, selPixIdxs, globalF, (motionR, motionC, motionZ), (onlineYshifts, onlineXshifts, onlineZshifts)

def create_parameter_gui():
    root = tk.Tk()
    root.title("SLAP2 Analysis Parameters")
    
    # Create main frame
    main_frame = ttk.Frame(root, padding="10")
    main_frame.grid(row=0, column=0, sticky=(tk.W, tk.E, tk.N, tk.S))
    
    # Parameter variables
    param_vars = {}
    
    # Create parameter inputs
    row = 0
    
    # # Discard Initial (s)
    # ttk.Label(main_frame, text="Discard Initial (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    # param_vars['discardInitial_s'] = tk.DoubleVar(value=0.1)
    # ttk.Entry(main_frame, textvariable=param_vars['discardInitial_s'], width=15).grid(row=row, column=1, pady=2)
    # row += 1
    
    # Analyze Hz
    ttk.Label(main_frame, text="Analyze Hz:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['analyzeHz'] = tk.IntVar(value=100)
    ttk.Entry(main_frame, textvariable=param_vars['analyzeHz'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Activity Channel
    ttk.Label(main_frame, text="Activity Channel:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['activityChannel'] = tk.IntVar(value=1)
    ttk.Entry(main_frame, textvariable=param_vars['activityChannel'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Decay Tau (s)
    ttk.Label(main_frame, text="Decay Tau (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['decayTau_s'] = tk.DoubleVar(value=0.05)
    ttk.Entry(main_frame, textvariable=param_vars['decayTau_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Baseline Window (s)
    ttk.Label(main_frame, text="Baseline Window (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['baselineWindow_s'] = tk.DoubleVar(value=4)
    ttk.Entry(main_frame, textvariable=param_vars['baselineWindow_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Denoise Window (s)
    ttk.Label(main_frame, text="Denoise Window (s):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['denoiseWindow_s'] = tk.DoubleVar(value=1)
    ttk.Entry(main_frame, textvariable=param_vars['denoiseWindow_s'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # dXY
    ttk.Label(main_frame, text="dXY:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['dXY'] = tk.IntVar(value=5)
    ttk.Entry(main_frame, textvariable=param_vars['dXY'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Sparse Factor (log scale)
    ttk.Label(main_frame, text="Sparse Factor (log):").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['sparse_fac_log'] = tk.DoubleVar(value=-3.0)
    ttk.Entry(main_frame, textvariable=param_vars['sparse_fac_log'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Operator
    ttk.Label(main_frame, text="Operator:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['operator'] = tk.StringVar(value='Maria Goeppert Mayer')
    ttk.Entry(main_frame, textvariable=param_vars['operator'], width=15).grid(row=row, column=1, pady=2)
    row += 1

    # Max Workers
    ttk.Label(main_frame, text="Max Workers:").grid(row=row, column=0, sticky=tk.W, pady=2)
    param_vars['max_workers'] = tk.IntVar(value=6)
    ttk.Entry(main_frame, textvariable=param_vars['max_workers'], width=15).grid(row=row, column=1, pady=2)
    row += 1
    
    # Add some spacing
    row += 1
    
    # Result variable
    result = {'params': None, 'cancelled': True}
    
    def on_ok():
        try:
            result['params'] = {
                # 'discardInitial_s': param_vars['discardInitial_s'].get(),
                'analyzeHz': param_vars['analyzeHz'].get(),
                'activityChannel': param_vars['activityChannel'].get(),
                'decayTau_s': param_vars['decayTau_s'].get(),
                'baselineWindow_s': param_vars['baselineWindow_s'].get(),
                'dXY': param_vars['dXY'].get(),
                'sparse_fac': torch.exp(torch.tensor(param_vars['sparse_fac_log'].get())),
                'denoiseWindow_s': param_vars['denoiseWindow_s'].get(),
                'operator': param_vars['operator'].get(),
                'max_workers': param_vars['max_workers'].get()
            }
            result['cancelled'] = False
            root.destroy()
        except Exception as e:
            messagebox.showerror("Error", f"Invalid parameter values: {e}")
    
    def on_cancel():
        result['cancelled'] = True
        root.destroy()
    
    # Buttons
    button_frame = ttk.Frame(main_frame)
    button_frame.grid(row=row, column=0, columnspan=2, pady=10)
    
    ttk.Button(button_frame, text="OK", command=on_ok).pack(side=tk.LEFT, padx=5)
    ttk.Button(button_frame, text="Cancel", command=on_cancel).pack(side=tk.LEFT, padx=5)
    
    # Update window to calculate required size
    root.update_idletasks()
    
    # Get the required size
    width = main_frame.winfo_reqwidth()  # Add padding
    height = main_frame.winfo_reqheight()  # Add padding
    
    # Set window size to fit contents
    root.geometry(f"{width}x{height}")
    
    # Center the window
    root.transient()
    root.grab_set()
    root.mainloop()
    
    return result

def main():
    # load SLAP2 data folder
    dr = filedialog.askdirectory(initialdir = 'Z:\\scratch\\ophys\\Michael', \
                                        title = "Select data directory")
    # dr = 'Z:\\scratch\\ophys\\Michael\\slap2_integration+raster\\slap2_760268_2024-11-05_12-35-49\\fov1\\experiment1'
    print(dr)

    # Show GUI and get parameters
    gui_result = create_parameter_gui()
    
    if gui_result['cancelled']:
        print("Parameter selection cancelled by user")
        return
    
    # Extract metadata parameters
    params = gui_result['params']
    print("selected params:")
    print(params)

    # find trialTable.mat file in data directory
    trialTableFile = os.path.join(dr, 'trialTable.mat')
    if not os.path.exists(trialTableFile):
        raise FileNotFoundError(f"Trial table file not found at: {trialTableFile}")
    print(trialTableFile)

    ''' trialTable structure
    trialTable = 

    struct with fields:

                    refStack: {[1×1 struct]  [1×1 struct]}
                    filename: {2×16 cell}
                    firstLine: [2×16 double]
                    lastLine: [2×16 double]
            trialEndTimeFromPC: [7.3956e+05 7.3956e+05 … ] (1×16 double)
        trialStartTimeInferred: [7.3956e+05 7.3956e+05 … ] (1×16 double)
                trueTrialIx: [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]
                        epoch: [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1]
                    lookupFile: 'Z:\scratch\ophys\Michael\slap2_integration+raster\slap2_760268_2024-11-05_12-35-49\fov1\experiment1\\integrationRegLookupTable.mat'
                    fnRegDSInt: {2×16 cell}
                    fnAdataInt: {2×16 cell}
                alignParamsInt: [1×1 struct]
                    fnRegDS: {2×16 cell}
                    fnAdata: {2×16 cell}
                alignParams: [1×1 struct]
    '''

    # Get the struct
    trialTable = spio.loadmat(trialTableFile)['trialTable'][0,0]

    nDMDs = trialTable['filename'].shape[0]
    nTrials = trialTable['filename'].shape[1]

    # Verify all required files exist
    keepTrials = np.ones((nDMDs, nTrials), dtype=bool)

    for trialIx in range(nTrials-1, -1, -1):
        for DMDix in range(nDMDs):
            # Check registration TIFF files
            # tiff_fn = os.path.splitext(os.path.basename(trialTable['fnRegDSInt'][DMDix,trialIx][0]))[0]
            # if not os.path.exists(os.path.join(dr, tiff_fn + '.tif')):
            #     print(f'Missing tiff file: {tiff_fn}')
            #     keepTrials[DMDix,trialIx] = False
                
            # Check alignment data files
            align_fn = os.path.splitext(os.path.basename(trialTable['fnAdataInt'][DMDix,trialIx][0]))[0]
            if not os.path.exists(os.path.join(dr, align_fn + '.mat')):
                print(f'Missing alignData file: {align_fn}')
                keepTrials[DMDix,trialIx] = False
                
            # Check source data files
            source_fn = trialTable['filename'][DMDix,trialIx][0]
            if not os.path.exists(os.path.join(dr, source_fn)):
                print(f'Missing source data file: {source_fn}')
                keepTrials[DMDix,trialIx] = False

    if not np.all(keepTrials):
        print(f'Files were missing for {np.sum(~keepTrials)} recordings; likely failed alignments. Proceeding without them.')

    if not np.any(keepTrials):
        raise RuntimeError('All trials were rejected due to missing alignment files!')

    firstValidTrial = np.where(keepTrials)[1][0]

    # Create ExperimentSummary directory if it doesn't exist
    savedr = os.path.join(dr, 'ExperimentSummary')
    if not os.path.exists(savedr):
        os.makedirs(savedr)
    output_h5_filename = os.path.join(savedr, f'experiment_summary_{datetime.datetime.now().strftime("%Y%m%d-%H%M%S")}.h5')

    # Load aData file
    aData = spio.loadmat(trialTable['fnAdataInt'][0,firstValidTrial][0])['aData'][0,0]
    
    params['numChannels'] = aData['numChannels'][0,0]
    params['alignHz'] = aData['alignHz'][0,0]

    # Get the lookup file path
    lookupFile = trialTable['lookupFile'][0]

    ''' Lookup table structure
    lookupTable = 

    struct with fields:

        likelihood_means: {[5-D single]  [5-D single]}
        allSuperPixelIDs: {[3786×1 uint32]  [2699×1 uint32]}
        sparseMaskInds: {[34102×2 uint32]  [24255×2 uint32]}
                    xPre: 25
                xPost: 25
                    yPre: 25
                yPost: 25
                    zPre: 10
                zPost: 10
    '''

    with h5py.File(lookupFile, 'r') as f:
        lt = f['lookupTable']
        refs = f['#refs#']
        
        # see if this can be made compatible for variable nDMDs
        # allSuperPixelIDs
        allSuperPixelIDs_refs = lt['allSuperPixelIDs'][:].flat
        allSuperPixelIDs = {'DMD1': refs[allSuperPixelIDs_refs[0]][:].T.astype(np.int32),
                            'DMD2': refs[allSuperPixelIDs_refs[1]][:].T.astype(np.int32)}
        
        # sparseMaskInds
        sparseMaskInds_refs = lt['sparseMaskInds'][:].flat
        sparseMaskInds = {'DMD1': refs[sparseMaskInds_refs[0]][:].T.astype(np.int32),
                        'DMD2': refs[sparseMaskInds_refs[1]][:].T.astype(np.int32)}
        
        fastZ2RefZ_refs = lt['fastZ2RefZ'][:].flat
        fastZ2RefZ = {'DMD1': refs[fastZ2RefZ_refs[0]][:].T.astype(np.int32),
                    'DMD2': refs[fastZ2RefZ_refs[1]][:].T.astype(np.int32)}

    del sparseMaskInds_refs, allSuperPixelIDs_refs, fastZ2RefZ_refs

    # Print shapes to verify
    print("Shapes:")
    for DMDix in range(nDMDs):
        print(f"allSuperPixelIDs DMD{DMDix+1}: {allSuperPixelIDs['DMD'+str(DMDix+1)].shape}")
        print(f"sparseMaskInds DMD{DMDix+1}: {sparseMaskInds['DMD'+str(DMDix+1)].shape}")

    # get integration reference stack
    refStack = {}
    for DMDix in range(nDMDs):
        pattern = f"**/*DMD{DMDix+1}_CONFIG2-REFERENCE*"
        matching_files = list(Path(dr).glob(pattern))
        if matching_files:
            first_file = str(matching_files[0])
            print(f"DMD{DMDix+1} ref stack file: {first_file}")
            numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
            refStackTmp = skimio.imread(first_file) / 100
            refStackTmp = refStackTmp.reshape(-1, numChannels, refStackTmp.shape[1], refStackTmp.shape[2]).transpose(1,0,2,3)
            refStack[f'DMD{DMDix+1}'] = refStackTmp
        else:
            pattern = f"**/*DMD{DMDix+1}-REFERENCE*"
            matching_files = list(Path(dr).glob(pattern))
            if matching_files:
                first_file = str(matching_files[0])
                print(f"DMD{DMDix+1} ref stack file: {first_file}")
                numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
                refStackTmp = skimio.imread(first_file) / 100
                refStackTmp = refStackTmp.reshape(-1, numChannels, refStackTmp.shape[1], refStackTmp.shape[2]).transpose(1,0,2,3)
                refStack[f'DMD{DMDix+1}'] = refStackTmp
            else:
                print(f"No matching files found for DMD{DMDix+1}")

    # Get dilation 17 PSF or load in from file
    pattern = f"**/*DMD*-REFERENCE*"
    matching_files = list(Path(dr).glob(pattern))
    if matching_files:
        first_file = str(matching_files[0])
        first_file_name = os.path.splitext(os.path.basename(first_file))[0].split('_DMD')[0]
        psf_path = os.path.join(dr, first_file_name + '-PSF.tif')
    else:
        psf_path = os.path.join(dr, 'refStack-PSF.tif')
        print(f"No ref stack found")

    if not os.path.exists(psf_path):
        psf = {}
        for DMDix in range(nDMDs):
            print(f"Calculating PSF for DMD{DMDix+1}")
            psf[f'DMD{DMDix+1}'] = skimio.imread(os.path.join(str(Path(__file__).parent), 'psfs', 'dil-17.tif')) # todo: make dilation a parameter
        combined_dims = [0,0]
        for DMDix in range(nDMDs):
            combined_dims[0] = max(psf[f'DMD{DMDix+1}'].shape[0], combined_dims[0])
            combined_dims[1] = max(psf[f'DMD{DMDix+1}'].shape[1], combined_dims[1])

        psf_combined = np.zeros((nDMDs, combined_dims[0], combined_dims[1]))
        for DMDix in range(nDMDs):
            psf_combined[DMDix] = np.pad(psf[f'DMD{DMDix+1}'], (((combined_dims[0] - psf[f'DMD{DMDix+1}'].shape[0])//2,(combined_dims[0] - psf[f'DMD{DMDix+1}'].shape[0])//2), ((combined_dims[1] - psf[f'DMD{DMDix+1}'].shape[1])//2,(combined_dims[1] - psf[f'DMD{DMDix+1}'].shape[1])//2)), constant_values = np.min(psf[f'DMD{DMDix+1}']))

        skimio.imsave(psf_path, psf_combined.astype(np.float32))
    else:
        psf_combined = skimio.imread(psf_path)
    
    psf_combined[0][psf_combined[0] < np.max(psf_combined[0])*np.exp(-3)] = 0
    psf_combined[1][psf_combined[1] < np.max(psf_combined[1])*np.exp(-3)] = 0

    psf = {
        'DMD1': psf_combined[0],
        'DMD2': psf_combined[1]
    }

    del psf_combined

    # Crop PSFs to remove boundary zeros
    for DMDix in range(nDMDs):
        # Find non-zero rows and columns
        non_zero_rows = np.any(psf[f'DMD{DMDix+1}'] != 0, axis=1)
        non_zero_cols = np.any(psf[f'DMD{DMDix+1}'] != 0, axis=0)
        
        # Get indices of first and last non-zero rows/cols
        row_start, row_end = np.where(non_zero_rows)[0][[0, -1]]
        col_start, col_end = np.where(non_zero_cols)[0][[0, -1]]
        
        # Crop the PSF
        psf[f'DMD{DMDix+1}'] = psf[f'DMD{DMDix+1}'][row_start:row_end+1, col_start:col_end+1]

    for DMDix in range(nDMDs): 
        # range(nDMDs-1, -1, -1):
        print(f'Processing DMD{DMDix+1}')

        dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
        dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
        numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
        numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]
        numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]

        nPixels = dmdPixelsPerColumn * dmdPixelsPerRow

        subsampleMatrixInds = np.zeros((numSuperPixels,2))
        for spIdx in range(numSuperPixels):
            currSpInds = np.where(sparseMaskInds[f'DMD{DMDix+1}'][:,1] == spIdx+1)[0]
            currSpOpenPixs = sparseMaskInds[f'DMD{DMDix+1}'][currSpInds,0]-1
            spRefPix = currSpOpenPixs[np.floor(len(currSpOpenPixs)/2).astype(int)]
            subsampleMatrixInds[spIdx,0] = spRefPix
            subsampleMatrixInds[spIdx,1] = spIdx+1

        refPixs = torch.from_numpy(subsampleMatrixInds[:,0])
        refD = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor').int()
        refC = torch.div((refPixs - refD * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor').int()
        refR = (refPixs % dmdPixelsPerColumn).int()

        filterSize = psf[f'DMD{DMDix+1}'].shape[0]*psf[f'DMD{DMDix+1}'].shape[1]
        
        sparseHInds = np.zeros((2,numSuperPixels*filterSize))
        sparseHVals = np.zeros((numSuperPixels*filterSize,))
        for spIdx in range(subsampleMatrixInds.shape[0]):
            tmpMap = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow))
            tmpMap[:] = np.nan

            tmpMap[refD[spIdx],
                    refR[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[0]//2+1,
                    refC[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

            sparseHInds[0,spIdx*filterSize:(spIdx+1)*filterSize] = subsampleMatrixInds[spIdx,1] - 1
            sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize] = np.where(~np.isnan(tmpMap.flatten()))[0]

            sparseHVals[spIdx*filterSize:(spIdx+1)*filterSize] = tmpMap.flatten()[sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize].astype(np.uint32)]
        non_zero_mask = sparseHVals != 0
        sparseHVals = sparseHVals[non_zero_mask]
        sparseHInds = sparseHInds[:, non_zero_mask]

        convMatrix = np.zeros((numSuperPixels,numSuperPixels))
        for spIdx in range(subsampleMatrixInds.shape[0]):
            tmpMap = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow))
            tmpMap[refD[spIdx],
                    refR[spIdx]-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx]+psf[f'DMD{DMDix+1}'].shape[0]//2+1,
                    refC[spIdx]-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx]+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

            convMatrix[spIdx,:] = tmpMap[refD,refR,refC]

        trial_info = [(i, keepTrials[DMDix,i]) for i in range(nTrials)]

        get_trial_data_partial = partial(
            get_trial_data,
            DMDix=DMDix,
            params=params,
            sampFreq=params['alignHz'],
            refStack=refStack,
            fastZ2RefZ=fastZ2RefZ,
            allSuperPixelIDs=allSuperPixelIDs,
            dr=dr,
            trialTable=trialTable
        )

        data_file = os.path.join(dr, f'lowres_data_DMD{DMDix+1}.npz')
        
        if os.path.exists(data_file):
            print(f'Loading existing low resolution data from {data_file}')
            data_arrays = np.load(data_file)
            lowResData = data_arrays['lowResData']
            lowResDataCt = data_arrays['lowResDataCt'] 
            lowResMotionR = data_arrays['lowResMotionR']
            lowResMotionC = data_arrays['lowResMotionC']
            lowResMotionZ = data_arrays['lowResMotionZ']
            lowResTrialID = data_arrays['lowResTrialID']
        else:
            with mp.Pool(processes=min(params['max_workers'],min(mp.cpu_count(),len(trial_info)))) as pool:
                results = list(pool.imap(get_trial_data_partial, trial_info))

            lowResData = np.concatenate([r[0] for r in results], axis=1)
            lowResDataCt = np.concatenate([r[1] for r in results], axis=1)
            # lowResBaseline = np.concatenate([r[2] for r in results], axis=1)
            lowResMotionR = np.concatenate([r[2]['motionDSr'] for r in results], axis=0)
            lowResMotionC = np.concatenate([r[2]['motionDSc'] for r in results], axis=0)
            lowResMotionZ = np.concatenate([r[2]['motionDSz'] for r in results], axis=0)
            # lowResBrightness = np.concatenate([r[2]['brightnessDS'] for r in results], axis=0)
            lowResTrialID = np.concatenate([np.ones_like(r[3])*i for i,r in enumerate(results)], axis=0)

            # Save data arrays
            data_arrays = {
                'lowResData': lowResData,
                'lowResDataCt': lowResDataCt,
                'lowResMotionR': lowResMotionR,
                'lowResMotionC': lowResMotionC,
                'lowResMotionZ': lowResMotionZ,
                'lowResTrialID': lowResTrialID
            }
            np.savez(data_file, **data_arrays)
            print(f'Saved low resolution data to {data_file}')

        uniqueMotion, motInds = np.unique(np.round(np.concatenate((lowResMotionR,lowResMotionC,lowResMotionZ),axis=1)),axis=0,return_inverse=True)
        motIndsToKeep = (np.bincount(motInds) > 100).nonzero()[0]

        uniqueMotionYX = np.unique(uniqueMotion[motIndsToKeep,:2],axis=0)
        uniqueMotionZ = np.unique(uniqueMotion[motIndsToKeep,2],axis=0)

        selPixMask = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), dtype=bool)
        for i in range(uniqueMotionYX.shape[0]):
            selPixMask[refD,refR + int(uniqueMotionYX[i,0]),refC + int(uniqueMotionYX[i,1])] = True
        selPixMask = ndimage.binary_dilation(selPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
        selPixIdxs = np.flatnonzero(selPixMask)

        # Get pixel coordinates for selected pixels
        pixel_coords = np.zeros((len(selPixIdxs), 3))  # [z, y, x] coordinates
        for i, pixel_idx in enumerate(selPixIdxs):
            # Convert linear index to 3D coordinates
            z = pixel_idx // (dmdPixelsPerColumn * dmdPixelsPerRow)
            remainder = pixel_idx % (dmdPixelsPerColumn * dmdPixelsPerRow)
            y = remainder // dmdPixelsPerRow
            x = remainder % dmdPixelsPerRow
            pixel_coords[i] = [z, y, x]
        
        # Convert to torch tensor
        pixel_coords_tensor = torch.tensor(pixel_coords, dtype=torch.float32)

        # seed sources
        psf_shape = psf[f'DMD{DMDix+1}'].shape
        psf_center = (psf_shape[0] // 2, psf_shape[1] // 2)

        # Convert selected pixel indices to 2D coordinates (vectorized)
        sel_pixels_2d = np.column_stack([
            (selPixIdxs % (dmdPixelsPerColumn * dmdPixelsPerRow)) // dmdPixelsPerRow,  # rows
            (selPixIdxs % (dmdPixelsPerColumn * dmdPixelsPerRow)) % dmdPixelsPerRow    # cols
        ])

        # Pre-allocate result matrix
        D = torch.zeros((len(selPixIdxs), len(selPixIdxs)), dtype=torch.float32)

        # Vectorized computation using broadcasting
        src_rows = sel_pixels_2d[:, 0][np.newaxis, :]  # Shape: (1, n_sources)
        src_cols = sel_pixels_2d[:, 1][np.newaxis, :]  # Shape: (1, n_sources)
        tgt_rows = sel_pixels_2d[:, 0][:, np.newaxis]  # Shape: (n_targets, 1)
        tgt_cols = sel_pixels_2d[:, 1][:, np.newaxis]  # Shape: (n_targets, 1)

        # Calculate relative positions
        rel_rows = tgt_rows - src_rows + psf_center[0]  # Shape: (n_targets, n_sources)
        rel_cols = tgt_cols - src_cols + psf_center[1]  # Shape: (n_targets, n_sources)

        # Create mask for valid PSF indices
        valid_mask = ((rel_rows >= 0) & (rel_rows < psf_shape[0]) & 
                    (rel_cols >= 0) & (rel_cols < psf_shape[1]))

        # Convert PSF to torch tensor and extract values where valid
        psf_tensor = torch.from_numpy(psf[f'DMD{DMDix+1}']).float()
        psf_tensor = psf_tensor / torch.sum(psf_tensor)
        D[valid_mask] = psf_tensor[rel_rows[valid_mask], rel_cols[valid_mask]]
        
        # Create expanded PSF with 5x resolution, maintaining center alignment
        ex_fac = 5
        psf_tensor_expanded = torch.zeros((psf_shape[0]*ex_fac, psf_shape[1]*ex_fac), dtype=torch.float32)

        # Calculate center points
        orig_center_y = (psf_shape[0] - 1) / 2
        orig_center_x = (psf_shape[1] - 1) / 2
        expanded_center_y = (psf_shape[0]*ex_fac - 1) / 2
        expanded_center_x = (psf_shape[1]*ex_fac - 1) / 2

        # Create coordinate grids centered around zero
        orig_y = torch.linspace(-expanded_center_y, expanded_center_y, psf_shape[0])
        orig_x = torch.linspace(-expanded_center_x, expanded_center_x, psf_shape[1])
        expanded_y = torch.linspace(-expanded_center_y, expanded_center_y, psf_shape[0]*ex_fac)
        expanded_x = torch.linspace(-expanded_center_x, expanded_center_x, psf_shape[1]*ex_fac)

        # Create meshgrids
        orig_Y, orig_X = torch.meshgrid(orig_y, orig_x, indexing='ij')
        expanded_Y, expanded_X = torch.meshgrid(expanded_y, expanded_x, indexing='ij')

        # Interpolate PSF values while maintaining center alignment
        interp_spline = RectBivariateSpline(orig_y.numpy(), orig_x.numpy(), psf_tensor.numpy())
        psf_tensor_expanded = torch.from_numpy(interp_spline(expanded_y.numpy(), expanded_x.numpy())).float()

        # Normalize expanded PSF
        psf_tensor_expanded = psf_tensor_expanded / torch.sum(psf_tensor_expanded)

        D_expanded = torch.zeros_like(D)

        psf_shape_expanded = psf_tensor_expanded.shape
        psf_center_expanded = (psf_shape_expanded[0] // 2, psf_shape_expanded[1] // 2)

        rel_rows_expanded = tgt_rows - src_rows + psf_center_expanded[0]
        rel_cols_expanded = tgt_cols - src_cols + psf_center_expanded[1]

        valid_mask_expanded = ((rel_rows_expanded >= 0) & (rel_rows_expanded < psf_shape_expanded[0]) & 
                            (rel_cols_expanded >= 0) & (rel_cols_expanded < psf_shape_expanded[1]))

        D_expanded[valid_mask_expanded] = psf_tensor_expanded[rel_rows_expanded[valid_mask_expanded], rel_cols_expanded[valid_mask_expanded]]

        lowResDataNorm = lowResData / lowResDataCt

        def build_interp_data(data, refR, refC, uniqueMotion, motInds,
                        H=800, W=1280, dtype=np.float32):
            # ---- normalize dtypes once
            refR = np.asarray(refR, dtype=np.int32)
            refC = np.asarray(refC, dtype=np.int32)
            uniqueMotion = np.asarray(uniqueMotion, dtype=np.int32)

            # ---- compute windows once (mirrors your logic)
            windowR = [max(0, refR.min()) + uniqueMotion[:,0].min(),
                    min(H, refR.max()) + uniqueMotion[:,0].max()]
            windowC = [max(0, refC.min()) + uniqueMotion[:,1].min(),
                    min(W, refC.max()) + uniqueMotion[:,1].max()]
            r0, r1 = int(windowR[0]), int(windowR[1])
            c0, c1 = int(windowC[0]), int(windowC[1])

            out = np.full((sel_pixels_2d.shape[0],data.shape[1]), np.nan, dtype=dtype)

            # ---- frames per motion bucket (avoid (motInds==idx) every time)
            frames_by_motion = {idx: np.flatnonzero(motInds == idx) for idx in np.unique(motInds)}

            for m_idx in np.unique(motInds):
                frames = frames_by_motion.get(m_idx)
                if frames is None or frames.size == 0:
                    continue

                # shift once per motion group
                dR, dC = int(uniqueMotion[m_idx, 0]), int(uniqueMotion[m_idx, 1])
                sR = refR + dR
                sC = refC + dC

                # columns that actually land in our window and have any samples
                in_win = (sR >= r0) & (sR <= r1) & (sC >= c0) & (sC <= c1)
                if not np.any(in_win):
                    continue
                cols = np.unique(sC[in_win])

                # ---- precompute per-column row positions and 1D data indices (sorted by row)
                rows_by_col = {}
                idxs_by_col = {}
                selPixIdxs_by_col = {}
                single_point = []  # track columns with only one sample (no interpolation)
                remaining_cols = np.unique(sel_pixels_2d[:,1])
                for c in cols:
                    mask = (sC == c)
                    nnz = mask.sum()
                    if nnz < 2:
                        if nnz == 1:
                            single_point.append(int(c))
                        continue
                    rows = sR[mask].astype(np.int32)
                    data_ix = np.flatnonzero(mask).astype(np.int32)
                    order = rows.argsort(kind="mergesort")
                    rows_by_col[int(c)] = rows[order]
                    idxs_by_col[int(c)]  = data_ix[order]
                    selPixIdxs_by_col[int(c)] = np.flatnonzero(sel_pixels_2d[:,1] == c)
                    remaining_cols = np.delete(remaining_cols, np.flatnonzero(remaining_cols == c))
                
                closest_col = [None] * len(remaining_cols)
                for i, c in enumerate(remaining_cols):
                    closest_col[i] = cols[np.argmin(np.abs(cols - c))]
                    selPixIdxs_by_col[int(c)] = np.flatnonzero(sel_pixels_2d[:,1] == c)
                    selPixIdxs_by_col[int(c)] = selPixIdxs_by_col[int(c)][np.isin(sel_pixels_2d[selPixIdxs_by_col[int(c)],0],sel_pixels_2d[selPixIdxs_by_col[int(closest_col[i])],0])]

                # ---- fill frames in this motion bucket
                for f in tqdm(frames):
                    yvec = data[:, f].astype(dtype, copy=False)

                    # interpolate only columns that have >= 2 samples
                    for c, rows in rows_by_col.items():
                        vals = yvec[idxs_by_col[c]]
                        colSelPixIdxs = selPixIdxs_by_col[c]
                        # interpolate just within our row window; outside stays NaN
                        out[colSelPixIdxs, f] = np.interp(sel_pixels_2d[colSelPixIdxs,0], rows, vals, left=vals[0], right=vals[-1])

                    # keep single points (no interpolation), like your original code
                    for c in single_point:
                        # (exact row/val for that one sample)
                        pos = np.flatnonzero((sC == c))[0]
                        rr = sR[pos]
                        if r0 <= rr <= r1:
                            out[np.flatnonzero((sel_pixels_2d[:,1] == c) & (sel_pixels_2d[:,0] == rr)), f] = yvec[pos]

                    for i, c in enumerate(remaining_cols):
                        closest_col_idx_mask = np.isin(sel_pixels_2d[selPixIdxs_by_col[closest_col[i]],0],sel_pixels_2d[selPixIdxs_by_col[c],0])
                        out[selPixIdxs_by_col[c], f] = out[selPixIdxs_by_col[closest_col[i]][closest_col_idx_mask], f]

            # out has shape: (n_frames, r1-r0+1, c1-c0+1), same as your bg_frame slice
            return out, (r0, r1), (c0, c1)
        
        interp_data, (r0,r1), (c0,c1) = build_interp_data(lowResDataNorm, refR, refC, uniqueMotion, motInds)
        # validFrames = ~np.all(np.isnan(interp_data),axis=0)
        # interp_data_valid = interp_data[:,validFrames]
        # validPix = np.all(~np.isnan(interp_data_valid),axis=1)
        # interp_data_valid = interp_data_valid[validPix]

        baseline_window = int(params['alignHz']*params['baselineWindow_s'])
        # nan_mask = np.isnan(interp_data)
        # Fill NaN pixels with median value of each row (vectorized)
        # row_medians = np.nanmedian(interp_data, axis=1, keepdims=True)
        # interp_data = np.where(nan_mask, row_medians, interp_data)
        interp_data_background = ndimage.uniform_filter1d(interp_data, size=baseline_window, axis=1, mode='nearest')
        # interp_data_background = ndimage.percentile_filter(interp_data, size=(baseline_window,1), mode='nearest', percentile=40)
        # interp_data_background = np.where(nan_mask,np.nan,interp_data_background)

        del interp_data
        # svd = TruncatedSVD(n_components=1)
        # U = svd.fit_transform(interp_data_valid)
        # V = svd.components_

        # U_im = np.zeros((800*1280,))
        # U_im[:] = np.nan

        # U_im[selPixIdxs] = U[:,0]
        # U_im = U_im.reshape((800,1280))

        # background_spatial_components = np.zeros((numSuperPixels,len(motIndsToKeep)))
        # for i in range(len(motIndsToKeep)):
        #     motion_idx = motIndsToKeep[i]
        #     dR, dC = int(uniqueMotion[motion_idx, 0]), int(uniqueMotion[motion_idx, 1])
        #     sR = refR + dR
        #     sC = refC + dC
        #     background_spatial_components[:,i] = U_im[sR,sC]
            # motion_frames = (motInds == motion_idx).nonzero()[0]
            # background_spatial_components[:,i] = torch.median(data_for_nmf[:, motion_frames],dim=1)[0]
            # background_spatial_components[:,i] = background_spatial_components[:,i] / torch.norm(background_spatial_components[:,i])

        # precompute H matrices
        H_mots = [None] * len(motIndsToKeep)
        for i, motion_idx in enumerate(motIndsToKeep):
            sparseHIndsShifted = sparseHInds.copy()
            sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)

            sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
            sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
            H_mots[i] = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)
        
        # Pre-compute sparse matrix multiplications with D for all motion indices
        HD_mots = [None] * len(H_mots)
        for i in range(len(H_mots)):
            HD = torch.sparse.mm(H_mots[i], D)
            HD = HD / torch.sum(HD,dim=0,keepdim=True)

            HD_expanded = torch.sparse.mm(H_mots[i], D_expanded)
            HD_expanded = HD_expanded / torch.sum(HD_expanded,dim=0,keepdim=True)

            HD_mots[i] = HD - HD_expanded

            validPixMask = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow), dtype=bool)
            validPixMask[refD,refR + int(uniqueMotion[motIndsToKeep[i],0]),refC + int(uniqueMotion[motIndsToKeep[i],1])] = True
            validPixMask = ndimage.binary_dilation(validPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]), dtype=bool))
            validPixMask = ndimage.binary_erosion(validPixMask, structure=np.ones((1,psf[f'DMD{DMDix+1}'].shape[0],psf[f'DMD{DMDix+1}'].shape[1]*2-1), dtype=bool))
            validPixIdxs = np.flatnonzero(validPixMask)

            HD_mots[i][:,~np.isin(selPixIdxs,validPixIdxs)] = 0
        
        # background = np.zeros_like(lowResDataNorm)
        # background[:] = np.nan
        # bg_trace = np.full((lowResDataNorm.shape[1],),np.nan)
        # for i, motion_idx in enumerate(motIndsToKeep):
        #     motion_frames = (motInds == motion_idx).nonzero()[0]
            
        #     # Get background component as contiguous vector
        #     bg_comp = background_spatial_components[:,i]
            
        #     # Compute background trace more efficiently
        #     bg_norm_sq = np.sum(bg_comp * bg_comp)
        #     bg_trace[motion_frames] = (bg_comp @ lowResDataNorm[:,motion_frames]) / bg_norm_sq

        #     background[:,motion_frames] = np.outer(bg_comp, bg_trace[motion_frames])
            
        background = np.zeros_like(lowResDataNorm)
        background[:] = np.nan
        
        for motion_idx in tqdm(np.unique(motInds)):
            motion_frames = (motInds == motion_idx).nonzero()[0]
            dR, dC = int(uniqueMotion[motion_idx, 0]), int(uniqueMotion[motion_idx, 1])
            sR = refR + dR
            sC = refC + dC

            shifted_indices = sR * 1280 + sC
            
            # Vectorized operation: create mapping from selPixIdxs to shifted positions
            # This avoids creating the full 800*1280 array for each frame
            for f in motion_frames:
                # Create sparse mapping directly
                bg_values = interp_data_background[:, f]
                # Map selected pixel values to their shifted positions
                background[:, f] = bg_values[np.searchsorted(selPixIdxs, shifted_indices)]

        residual = lowResDataNorm - background
        del interp_data_background

        decayTau_frames = params['decayTau_s']*params['alignHz']
        decay_kernel = np.exp(np.linspace(-np.ceil(decayTau_frames*3),0,int(np.ceil(decayTau_frames*3)+1))/decayTau_frames)
        residualFilt = signal.convolve2d(residual,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same')
        # data_for_nmf = torch.from_numpy(data_for_nmf.astype(np.float32))

            
        rho = np.zeros((lowResData.shape[1],len(selPixIdxs)))
        rho[:] = np.nan
        for i, motion_idx in enumerate(motIndsToKeep):
            motion_frames = (motInds == motion_idx).nonzero()[0]
            # Use pre-computed HD matrix
            rho[motion_frames,:] = residualFilt[:,motion_frames].T @ HD_mots[i].numpy()

        # Process frames more efficiently by avoiding unnecessary allocations
        n_frames = rho.shape[0]
        batch_size = min(1000, n_frames)  # Smaller batch size to reduce memory
        num_batches = int(np.ceil(n_frames / batch_size))
        
        # Pre-compute 2D coordinates once
        spatial_coords = np.unravel_index(selPixIdxs, (dmdPixelsPerColumn, dmdPixelsPerRow))
        row_indices = spatial_coords[0][None, :]  # Shape: (1, len(selPixIdxs))
        col_indices = spatial_coords[1][None, :]  # Shape: (1, len(selPixIdxs))
        
        # Pre-compute dilation structure
        dilation_struct = np.ones((1,7,7))
        
        # Initialize output array
        rho_activity = np.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))
        temporal_pad = 1
        
        for batch_idx in tqdm(range(num_batches), desc="Creating activity map"):
            # Calculate batch bounds
            batch_start = batch_idx * batch_size 
            batch_end = min(batch_start + batch_size, n_frames)
            padded_start = max(0, batch_start - temporal_pad)
            padded_end = min(n_frames, batch_end + temporal_pad)
            curr_size = padded_end - padded_start
            
            # Allocate batch arrays - use float32 for memory efficiency
            batch_rho = np.zeros((curr_size, dmdPixelsPerColumn, dmdPixelsPerRow), dtype=np.float32)
            
            # Fill batch data more efficiently using advanced indexing
            time_indices = np.arange(curr_size)[:, None]
            batch_rho[time_indices, row_indices, col_indices] = rho[padded_start:padded_end, :].astype(np.float32)
            batch_rho[batch_rho == 0] = np.nan

            # Handle NaN values
            nan_mask = np.isnan(batch_rho)
            batch_rho[nan_mask] = 0
            dilated_nan_mask = fast_dilation(nan_mask, dilation_struct[0])

            # Calculate local maxima directly without intermediate mask array
            local_maxima = np.zeros_like(batch_rho, dtype=bool)
            local_maxima[1:-1,1:-1,1:-1] = (batch_rho[1:-1,1:-1,1:-1] > batch_rho[:-2,1:-1,1:-1]) & \
                                          (batch_rho[1:-1,1:-1,1:-1] > batch_rho[2:,1:-1,1:-1]) & \
                                          (batch_rho[1:-1,1:-1,1:-1] > batch_rho[1:-1,:-2,1:-1]) & \
                                          (batch_rho[1:-1,1:-1,1:-1] > batch_rho[1:-1,2:,1:-1]) & \
                                          (batch_rho[1:-1,1:-1,1:-1] > batch_rho[1:-1,1:-1,:-2]) & \
                                          (batch_rho[1:-1,1:-1,1:-1] > batch_rho[1:-1,1:-1,2:])
            
            # Combine conditions and compute activity in one step
            valid_maxima = local_maxima & ~dilated_nan_mask
            start_offset = batch_start - padded_start
            end_offset = start_offset + (batch_end - batch_start)
            rho_activity += np.sum(np.where(valid_maxima[start_offset:end_offset], 
                                          batch_rho[start_offset:end_offset]**3, 0), 
                                 axis=0)

        act_im = rho_activity.copy()
        nan_mask = (act_im == 0) | (np.isnan(act_im))
        act_im[nan_mask] = np.nan

        med_act_im = ndimage.generic_filter(act_im, np.nanmedian, size=(15,15))
        act_im = act_im - med_act_im
        act_im[nan_mask] = np.nan
        
        act_im[nan_mask] = np.nanmedian(act_im.flatten())
        act_im_filt = ndimage.gaussian_filter(act_im, sigma=[0.75, 0.75])
        act_im_filt[nan_mask] = np.nan
        act_im[nan_mask] = np.nan

        act_im = act_im_filt

        # act_im_hp = act_im - ndimage.gaussian_filter(act_im, sigma=[10*i for i in psf_shape])

        # Create a local maximum filter with a footprint of 3x3 pixels
        local_max = ndimage.maximum_filter(act_im, size=3)

        maxima_mask = (act_im == local_max) & ~np.isnan(act_im)

        # Get coordinates and values of local maxima
        maxima_coords = np.where(maxima_mask)
        maxima_values = act_im[maxima_mask]

        # Sort maxima by value in descending order
        sort_idx = np.argsort(-maxima_values)
        maxima_coords = (maxima_coords[0][sort_idx], maxima_coords[1][sort_idx])
        maxima_values = maxima_values[sort_idx]

        mad = np.median(np.abs(act_im[~np.isnan(act_im)] - np.median(act_im[~np.isnan(act_im)])))
        thresh = 3 * 1.4826 * mad

        # Remove maxima that are within 3 pixels of nan values, unless they are in top 5%
        valid_maxima = maxima_values > thresh
        maxima_coords = (maxima_coords[0][valid_maxima], maxima_coords[1][valid_maxima])
        maxima_values = maxima_values[valid_maxima]
        source_seeds = np.array(maxima_coords).T

        nSources = source_seeds.shape[0]

        def sel_pix_gaussian_profile(gaussian_params):
            y_means = gaussian_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]
            x_means = gaussian_params[:, 1].unsqueeze(0)  # Shape: [1, nSources]
            y_sigmas = gaussian_params[:, 2].unsqueeze(0)  # Shape: [1, nSources]
            x_sigmas = gaussian_params[:, 3].unsqueeze(0)  # Shape: [1, nSources]
            corr_coef = torch.tanh(gaussian_params[:, 4].unsqueeze(0))  # Shape: [1, nSources]

            # Center the coordinates
            y_centered = (pixel_coords_tensor[:, 1].unsqueeze(1) - y_means)  # Shape: [nPixels, nSources]
            x_centered = (pixel_coords_tensor[:, 2].unsqueeze(1) - x_means)  # Shape: [nPixels, nSources]

            # Compute terms for bivariate Gaussian with correlation
            z_score_y = y_centered / y_sigmas
            z_score_x = x_centered / x_sigmas
            
            # Full bivariate Gaussian formula with correlation
            exponent = (-1 / (2 * (1 - corr_coef**2))) * (
                z_score_y**2 - 
                2 * corr_coef * z_score_x * z_score_y + 
                z_score_x**2
            )

            profile = torch.exp(exponent)
            profile = profile * (torch.sqrt(z_score_y**2 + z_score_x**2) <= 3).float()

            return profile / torch.sum(profile,dim=0,keepdim=True)

        def sel_pix_patch_profile(patch_params):
            y_means = patch_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]
            x_means = patch_params[:, 1].unsqueeze(0)  # Shape: [1, nSources]
            y_radii = patch_params[:, 2].unsqueeze(0)  # Shape: [1, nSources]
            x_radii = patch_params[:, 3].unsqueeze(0)  # Shape: [1, nSources]

            # Center the coordinates
            y_centered = (pixel_coords_tensor[:, 1].unsqueeze(1) - y_means)  # Shape: [nPixels, nSources]
            x_centered = (pixel_coords_tensor[:, 2].unsqueeze(1) - x_means)  # Shape: [nPixels, nSources]

            profile = (torch.abs(y_centered) < y_radii) & (torch.abs(x_centered) < x_radii)

            return profile
        
        source_params = torch.cat([
            torch.tensor(source_seeds, dtype=torch.float32), # x,y means
            torch.ones(nSources, 2, dtype=torch.float32), # x,y sigmas
            torch.zeros(nSources, 1, dtype=torch.float32) # invtanh of correlation / tilt
        ], dim=1)

        A_patches = torch.zeros((nPixels, nSources), dtype=torch.bool)
        A_patches[selPixIdxs,:] = sel_pix_patch_profile(torch.cat([torch.tensor(source_seeds, dtype=torch.float32),
                                                            params['dXY'] * torch.ones(nSources, 2, dtype=torch.float32)],
                                                            dim=1))

        X_support_mots = [None] * len(motIndsToKeep)
        for i, motion_idx in enumerate(motIndsToKeep):
            X_support_mots[i] = torch.sparse.mm(H_mots[i], A_patches[selPixIdxs,:].float()) > 0

        # initialize sources as gaussian blobs
        A = torch.zeros((nPixels, nSources))
        A[selPixIdxs,:] = sel_pix_gaussian_profile(source_params * torch.tensor([1, 1, 1, 1, 1]))
        A[~A_patches] = 0
        sources_total_mass = torch.sum(A, dim=0, keepdim=True)
        A = torch.where(sources_total_mass > 0, A / sources_total_mass, A)

        # NMF parameters
        outer_loop_iters = 10
        
        mult_nmf_max_iters = 10
        nmf_tol = 1e-6

        # Gaussian fitting optimization parameters
        learning_rate = 0.01
        num_epochs = len(motIndsToKeep) * 5
        gd_tol = 1e-4
        
        # phi_lowRes = torch.zeros(lowResData.shape[1], nSources+1, dtype=torch.float32)
        # phi_lowRes[:] = np.nan

        phi_lowRes = torch.full((lowResData.shape[1], nSources), np.nan, dtype=torch.float32)

        X_mots = [None] * len(motIndsToKeep)
        overall_losses = [0] * (outer_loop_iters+1)

        data_for_nmf = torch.from_numpy(residualFilt.astype(np.float32))
        # data_for_nmf = torch.from_numpy(signal.convolve2d(lowResDataNorm,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same').astype(np.float32))
        # background_spatial_components = torch.from_numpy(background_spatial_components.astype(np.float32)) / torch.norm(background_spatial_components)
        
        for outer_loop_iter in range(outer_loop_iters): # outer loop

            # get phi_lowRes for current spatial profiles
            for i, motion_idx in enumerate(motIndsToKeep):
                motion_frames = (motInds == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # print warning for blank spatial components
                norms = torch.norm(X,dim=0,keepdim=True)
                if torch.any(norms == 0):
                    print(f"Warning: {np.nonzero(norms.numpy().squeeze() == 0)[0]} norms are zero for motion {i}")

                # least squares to fit phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T

                overall_losses[outer_loop_iter] += torch.sum((data_for_nmf[:, motion_frames] - X @ phi_lowRes[motion_frames,:].T) ** 2).item()

            shuffled_indices = torch.randperm(len(motIndsToKeep))
            for idx in tqdm(shuffled_indices,desc=f'Multiplicative NMF per motion displacement'): # loop over all motion displacements in random order
                i = idx.item()
                motion_idx = motIndsToKeep[i]
                motion_frames = (motInds == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                X[~X_support_mots[i]] = 0

                # add background spatial component
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # normalize spatial components
                X = X / torch.norm(X,dim=0,keepdim=True)

                # projected least squares to initialize phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]  # This gives all time points at once
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
                phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:], min=0)
            
                # Initialize variables for NMF
                error_values = []
                prev_reconstruction_error = float('inf')

                for iter_idx in range(mult_nmf_max_iters): # iterative fit in sp space with NMF
                    # Multiplicative update for NMF
                    # Update phi (temporal components) using multiplicative update rule
                    # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
                    numerator = X.T @ data_for_nmf[:, motion_frames]
                    denominator = (X.T @ X) @ phi_lowRes[motion_frames,:].T + 1e-10
                    phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:] * (numerator / denominator).T, min=0)
                    
                    # Update X (spatial components) using multiplicative update rule
                    # X = X * (data * phi^T) / (X * phi * phi^T + epsilon)
                    numerator = data_for_nmf[:, motion_frames] @ phi_lowRes[motion_frames,:]
                    denominator = X @ (phi_lowRes[motion_frames,:].T @ phi_lowRes[motion_frames,:]) + 1e-10
                    X = torch.clamp(X * (numerator / denominator), min=0)

                    # Apply spatial support constraint
                    X = X * X_support_mots[i].float()
                    # X[:,:-1] = X[:,:-1] * X_support_mots[i].float()
                    # X[:,-1] = background_spatial_components[:,i]
                    
                    # Normalize X and phi to avoid scaling ambiguity
                    norms = torch.norm(X, dim=0, keepdim=True)
                    X = torch.where(norms > 0, X / norms, X)
                    phi_lowRes[motion_frames,:] = torch.where(norms > 0, phi_lowRes[motion_frames,:] * norms, phi_lowRes[motion_frames,:])
                    
                    # Calculate reconstruction error
                    # reconstruction = X @ phi_lowRes[motion_frames,:].T
                    # current_error = torch.mean((data_for_nmf[:, motion_frames] - reconstruction)**2).item()
                    # error_values.append(current_error)
                    
                    # # Check for convergence
                    # if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    #     break
                        
                    # prev_reconstruction_error = current_error

                    if iter_idx % 3 == 0:   # sparsify step
                        X = X / torch.max(X, dim=0, keepdim=True)[0]
                        X = torch.clamp(X, min=params['sparse_fac']) - params['sparse_fac']

                        # re-normalize X after sparsification
                        norms = torch.norm(X, dim=0, keepdim=True)
                        X = torch.where(norms > 0, X / norms, X)
                
                X_mots[i] = X
            
            # Initialize Adam optimizer
            optim_loc_params = source_params[:,:2].clone().requires_grad_(True)
            optim_scale_params = source_params[:,2:4].clone().requires_grad_(True)
            optim_tilt_params = source_params[:,4].unsqueeze(1).clone().requires_grad_(True)
            optimizer = torch.optim.Adam([{'params': optim_loc_params, 'lr': 10 * learning_rate},
                                            {'params': optim_scale_params, 'lr': 0.1 * learning_rate},
                                            {'params': optim_tilt_params, 'lr': 0.1 * learning_rate}])
            
            # Gradient descent optimization for Gaussian parameters
            losses = []
            for epoch in tqdm(range(num_epochs),desc='Fitting Gaussian spatial profiles'): # gradient descent to fit pixel space profile to sp profile
                # Zero gradients
                optimizer.zero_grad()
                
                # calculate spatial profile
                A_step_sel_pix = sel_pix_gaussian_profile(torch.cat([optim_loc_params, optim_scale_params, optim_tilt_params], dim=1))
                A_step_sel_pix = A_step_sel_pix * A_patches[selPixIdxs,:].float()
                sources_total_mass = torch.sum(A_step_sel_pix, dim=0, keepdim=True)
                A_step_sel_pix = torch.where(sources_total_mass > 0, A_step_sel_pix / sources_total_mass, A_step_sel_pix)
                
                if epoch % len(motIndsToKeep) == 0:
                    shuffled_indices = torch.randperm(len(motIndsToKeep))

                i = shuffled_indices[epoch % len(motIndsToKeep)].item()
                motion_idx = motIndsToKeep[i]

                # Compute superpixel spatial profile
                X_step = torch.sparse.mm(H_mots[i], A_step_sel_pix)
                norms = torch.norm(X_step, dim=0, keepdim=True)
                X_step = torch.where(norms > 0, X_step / norms, X_step)
            
                # Compute loss across all sources
                # loss = torch.sum((X_step - X_mots[i][:,:-1]) ** 2)
                loss = torch.sum((X_step - X_mots[i]) ** 2)

                loss.backward()
                optimizer.step()
                
                # Apply constraints after optimizer step
                with torch.no_grad():
                    # losses.append(loss.item())
                    optim_loc_params.clamp_(min=torch.from_numpy(source_seeds - params['dXY']).float(), max=torch.from_numpy(source_seeds + params['dXY']).float())
                    optim_scale_params[:,0].clamp_(min=0.3, max=5)
                    optim_scale_params[:,1].clamp_(min=0.3, max=5)

                if epoch % len(motIndsToKeep) == len(motIndsToKeep)-1:
                    total_loss = 0
                    for i, motion_idx in enumerate(motIndsToKeep):
                        X_step = torch.sparse.mm(H_mots[i], A_step_sel_pix)
                        norms = torch.norm(X_step, dim=0, keepdim=True)
                        X_step = torch.where(norms > 0, X_step / norms, X_step)

                        # total_loss += torch.sum((X_step - X_mots[i][:,:-1]) ** 2)
                        total_loss += torch.sum((X_step - X_mots[i]) ** 2)

                    # Check for convergence
                    losses.append(total_loss.item())
                    if epoch // len(motIndsToKeep) > 0 and abs(losses[-1] - losses[-2]) < gd_tol:
                        print(f"Converged at epoch {epoch+1} with loss: {loss.item():.8f}")
                        break

            # Update parameters
            source_params = torch.cat([optim_loc_params.detach(), optim_scale_params.detach(), optim_tilt_params.detach()], dim=1)

            A = torch.zeros((nPixels, nSources), dtype=torch.float32)            
            A[selPixIdxs,:] = sel_pix_gaussian_profile(source_params * torch.tensor([1, 1, 1, 1, 1]))
            A[~A_patches] = 0
            sources_total_mass = torch.sum(A, dim=0, keepdim=True)
            A = torch.where(sources_total_mass > 0, A / sources_total_mass, A)

            # get phi_lowRes for current spatial profiles
            for i, motion_idx in enumerate(motIndsToKeep):
                motion_frames = (motInds == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # print warning for blank spatial components
                norms = torch.norm(X,dim=0,keepdim=True)
                if torch.any(norms == 0):
                    print(f"Warning: {np.nonzero(norms.numpy().squeeze() == 0)[0]} norms are zero for motion {i}")

                # least squares to fit phi_lowRes
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
            
            # sort sources by variance
            sortorder = np.argsort(-np.nansum((phi_lowRes[:,:nSources].numpy()-np.nanmean(phi_lowRes[:,:nSources].numpy(),axis=0))**2,axis=0))
            source_params = source_params[sortorder,:]
            source_seeds = source_seeds[sortorder,:]
            A = A[:,sortorder]
            A_patches = A_patches[:,sortorder]
            X_support_mots = [support[:,sortorder] for support in X_support_mots]
            # sortorder = np.append(sortorder,max(sortorder)+1)
            phi_lowRes = phi_lowRes[:,sortorder]

            if (outer_loop_iter+1) % 4 == 3: # prune sources
                residual = torch.zeros_like(data_for_nmf)
                residual[:] = np.nan
                for i in range(len(motIndsToKeep)):
                    motion_frames = (motInds == motIndsToKeep[i]).nonzero()[0]
                    X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                    # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)
                    residual[:,motion_frames] = data_for_nmf[:,motion_frames] - X @ phi_lowRes[motion_frames,:].T

                varExpSource = torch.zeros(nSources)
                varResidual = torch.zeros(nSources)
                for j in range(nSources-1,-1,-1):
                    # find pixels that correspond to half mass of spatial profile
                    # sorted_vals = torch.sort(A[selPixIdxs,j], descending=True)[0]
                    # cumsum_vals = torch.cumsum(sorted_vals, dim=0)
                    # total_mass = cumsum_vals[-1]
                    # half_mass_idx = torch.searchsorted(cumsum_vals, total_mass * 0.5)
                    # thresh = sorted_vals[half_mass_idx]
                    # validPixs = (A[selPixIdxs,j] >= thresh).nonzero()[:,0]
                    # tmp_A = A[selPixIdxs,j].unsqueeze(1).clone()
                    # tmp_A[~validPixs] = 0

                    for i in range(len(motIndsToKeep)):
                        motion_frames = (motInds == motIndsToKeep[i]).nonzero()[0]

                        # X = torch.sparse.mm(H_mots[i], tmp_A)
                        X = torch.sparse.mm(H_mots[i], A[selPixIdxs,j].unsqueeze(1))

                        contributingPixs = (X > 0).nonzero()[:,0]
                        varExpSource[j] += torch.sum(torch.sum(X[contributingPixs] * phi_lowRes[motion_frames,j].unsqueeze(-1).T,dim=0)**2)
                        varResidual[j] += torch.sum(torch.sum(residual[contributingPixs][:,motion_frames],dim=0)**2)
                
                SNR = varExpSource / varResidual
                
                SNR_cut = 1/3
                keepSources = (SNR > SNR_cut).nonzero()[:,0]
                print(f"Keeping {keepSources.shape[0]} of {nSources} sources")
                nSources = keepSources.shape[0]
                source_params = source_params[keepSources,:]
                source_seeds = source_seeds[keepSources,:]
                A = A[:,keepSources]
                A_patches = A_patches[:,keepSources]
                X_support_mots = [support[:,keepSources] for support in X_support_mots]
                # keepSources = np.append(keepSources,phi_lowRes.shape[1]-1)
                phi_lowRes = phi_lowRes[:,keepSources]
        
        # get phi_lowRes for current spatial profiles
        for i, motion_idx in enumerate(motIndsToKeep):
            motion_frames = (motInds == motion_idx).nonzero()[0]

            # project image space (A) into superpixel space (X)
            X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
            # X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

            overall_losses[outer_loop_iters] += torch.sum((data_for_nmf[:, motion_frames] - X @ phi_lowRes[motion_frames,:].T) ** 2).item()

        trial_info = [(i, keepTrials[DMDix,i], background[:,lowResTrialID == i]) for i in range(nTrials)]

        get_high_res_traces_partial = partial(get_high_res_traces,
                                              DMDix=DMDix,
                                              params=params,
                                              sampFreq=params['analyzeHz'],
                                              refStack=refStack,
                                              subsampleMatrixInds=subsampleMatrixInds,
                                              fastZ2RefZ=fastZ2RefZ,
                                              sparseHInds=sparseHInds,
                                              sparseHVals=sparseHVals,
                                              allSuperPixelIDs=allSuperPixelIDs,
                                              dr=dr, 
                                              trialTable=trialTable, 
                                              A_final=A,
                                              uniqueMotionDS=uniqueMotion,
                                              motIndsToKeepDS=motIndsToKeep,
                                              # background=background,
                                              # background_spatial_components=background_spatial_components,
                                              psf=psf)

        with mp.Pool(processes=min(params['max_workers'],min(mp.cpu_count(),len(trial_info)))) as pool:
            results = list(pool.imap(get_high_res_traces_partial, trial_info))

        # Save data to HDF5 file
        with h5py.File(output_h5_filename, 'a') as f:
            # Delete group if it exists
            group_name = f'DMD{DMDix+1}'
            if group_name in f:
                del f[group_name]
            
            # Create group and add datasets
            dmd_group = f.create_group(group_name)
            
            # Create subgroups for trial data
            source_group = dmd_group.create_group('sources')

            spatial_group = source_group.create_group('spatial')

            # spatial_group.create_dataset('source_params', data=source_params.numpy())
            spatial_group.create_dataset('profiles', data=A.numpy().reshape(dmdPixelsPerRow,dmdPixelsPerColumn,-1).transpose(2,0,1))
            spatial_group.create_dataset('coords', data=source_params[:,:2].numpy())

            temporal_group = source_group.create_group('temporal')

            dF = np.concatenate([r[0] for r in results], axis=0)
            F0 = np.concatenate([r[1] for r in results], axis=0)
            F = dF + F0

            F0 = compute_f0(F, int(np.ceil(params['denoiseWindow_s']*params['analyzeHz'])), int(np.ceil(params['baselineWindow_s']*params['analyzeHz'])))
            dF = F - F0
            temporal_group.create_dataset('dF', data=dF)
            temporal_group.create_dataset('F0', data=F0)
            # temporal_group.create_dataset('F', data=F)

            # trial_start_idxs = np.concatenate([[0], np.cumsum([len(r[2]) for r in results])[:-1]])
            trial_num_frames = np.concatenate([[len(r[2])] for r in results])

            frame_group = dmd_group.create_group('frame_info')
            # frame_group.create_dataset('trial_start_idxs', data=trial_start_idxs)
            frame_group.create_dataset('trial_num_frames', data=trial_num_frames)
            frame_group.create_dataset('discard_frames', data=np.any(np.isnan(F), axis=1))

            frame_group.create_dataset('offlineXshifts', data=np.concatenate([r[5][1] for r in results], axis=0))
            frame_group.create_dataset('offlineYshifts', data=np.concatenate([r[5][0] for r in results], axis=0))
            frame_group.create_dataset('offlineZshifts', data=np.concatenate([r[5][2] for r in results], axis=0))

            frame_group.create_dataset('onlineXshifts', data=np.concatenate([r[6][1] for r in results], axis=0))
            frame_group.create_dataset('onlineYshifts', data=np.concatenate([r[6][0] for r in results], axis=0))
            frame_group.create_dataset('onlineZshifts', data=np.concatenate([r[6][2] for r in results], axis=0))

            # frame_group.create_dataset('selPixIdxs', data=[r[2] for r in results])

            globalF = np.concatenate([r[4] for r in results], axis=0)
            global_group = dmd_group.create_group('global')
            global_group.create_dataset('F', data=globalF)

            visualizations = dmd_group.create_group('visualizations')
            visualizations.create_dataset('act_im', data=act_im)

            ref_stack_channels = trialTable['refStack'][0,DMDix]['channels'][0,0].squeeze().reshape((-1,))
            num_ref_stack_channels = len(ref_stack_channels)
            
            ref_stack = trialTable['refStack'][0,DMDix]['IM'][0,0].T
            ref_stack_reshaped = np.zeros((num_ref_stack_channels,ref_stack.shape[0] // num_ref_stack_channels,ref_stack.shape[1], ref_stack.shape[2]))
            for i in range(ref_stack.shape[0]):
                c = i % num_ref_stack_channels
                ref_stack_reshaped[c,i//num_ref_stack_channels,:,:] = ref_stack[i,:,:]
            del ref_stack
            
            ref_stack_dataset = visualizations.create_dataset('ref_stack', data=ref_stack_reshaped)
            ref_stack_dataset.attrs['channels'] = ref_stack_channels

        print(f"Added DMD{DMDix+1} data to {output_h5_filename}")

if __name__ == '__main__':
    main()