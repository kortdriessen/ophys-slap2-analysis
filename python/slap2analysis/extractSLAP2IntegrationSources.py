import sys
from pathlib import Path
import numpy as np
import torch
import scipy.io as spio
import scipy.ndimage as ndimage
from scipy.interpolate import RectBivariateSpline
import os
import time
import multiprocessing as mp
from functools import partial
import tkinter as tk
from tkinter import filedialog
import datetime
import importlib
import skimage.io as skimio
import h5py
from scipy import stats, signal
import slap2_utils

from tqdm import tqdm 

sys.path.append('C:/Users/michael.xie/Documents/ophys-slap2-analysis/python')
# sys.path.append(str(Path(__file__).parent.parent))
import reconstruct

import cv2

def fast_dilation(mask, kernel, iterations=1):
    if kernel is None:
        kernel = np.ones((7, 7), np.uint8)  # Structuring element for XY dilation
    out = np.empty_like(mask)
    
    # Iterate over Z slices (first dimension)
    for z in range(mask.shape[0]):
        out[z] = cv2.dilate(mask[z].astype(np.uint8), kernel, iterations=iterations)
    
    return out.astype(bool)

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
                allSuperPixelIDs, dr, trialTable, A_final, psf):
    trialIx, keepTrial = trial_info
    
    if not keepTrial:
        return [], [], []

    data_file = os.path.join(params['savedr'], f'trial_data_DMD{DMDix+1}_trial{trialIx}.npz')
    if os.path.exists(data_file):
        print(f'Loading existing trial data from {data_file}')
        data_arrays = np.load(data_file)
        dataNonNorm = data_arrays['dataNonNorm']
        dataCt = data_arrays['dataCt']
        DSframes = data_arrays['DSframes']

        aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]
        aData['DSframes'] = data_arrays['aData_DSframes']
    else:
        dataNonNorm, dataCt, aData, DSframes = get_trial_data(trial_info, DMDix, params, sampFreq, refStack, fastZ2RefZ, allSuperPixelIDs, dr, trialTable)
        # Save data arrays
        data_arrays = {
            'dataNonNorm': dataNonNorm,
            'dataCt': dataCt,
            'DSframes': DSframes,
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

    motionR = np.interp(DSframes, aData['DSframes'].squeeze(), aData['motionDSr'].squeeze())
    motionC = np.interp(DSframes, aData['DSframes'].squeeze(), aData['motionDSc'].squeeze())
    motionZ = np.interp(DSframes, aData['DSframes'].squeeze(), aData['motionDSz'].squeeze())

    uniqueMotion, motInds = np.unique(np.round(np.stack((motionR, motionC, motionZ), axis=1)), axis=0, return_inverse=True)
    motIndsToKeep = (np.bincount(motInds) > 100).nonzero()[0]

    uniqueMotionYX = np.unique(uniqueMotion[:,:2],axis=0)

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
    phi = torch.zeros(data.shape[1], nSources+1, dtype=torch.float32)
    phi[:] = np.nan
    
    for motion_idx in motIndsToKeep:

        # Extract data for the most common motion mode
        motion_frames = (motInds == motion_idx).nonzero()[0]
        
        data_tensor = torch.from_numpy(data[:, motion_frames].astype(np.float32))

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)
        sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

        # project image space (A) into superpixel space (X)
        X = torch.sparse.mm(H, A_final[selPixIdxs,:])
        X = X / torch.max(X, dim=0, keepdim=True)[0]
        X = torch.clamp(X, min=torch.exp(torch.tensor(-3.0))) - torch.exp(torch.tensor(-3.0))

        # newR = refR+uniqueMotion[motion_idx,0].astype(int)
        # newC = refC+uniqueMotion[motion_idx,1].astype(int)
        # newD = torch.from_numpy(fastZ2RefZ[f'DMD{DMDix+1}'][refD.numpy()].flatten()-1)+uniqueMotion[motion_idx,2].astype(int)

        # add background spatial component
        # background_spatial_component = torch.as_tensor(refStack[f'DMD{DMDix+1}'][0][newD,newR,newC].reshape((-1,1)))
        background_spatial_component = torch.median(data_tensor,dim=1)[0]
        background_spatial_component = background_spatial_component / torch.norm(background_spatial_component)
        X = torch.concat((X,background_spatial_component.unsqueeze(-1)),dim=1)

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

    globalF = np.sum(data, axis=0)

    return phi.numpy(), DSframes, selPixIdxs, globalF

def main():
    # load SLAP2 data folder
    dr = filedialog.askdirectory(initialdir = 'Z:\\scratch\\ophys\\Michael', \
                                        title = "Select data directory")
    # dr = 'Z:\\scratch\\ophys\\Michael\\slap2_integration+raster\\slap2_760268_2024-11-05_12-35-49\\fov1\\experiment1'
    print(dr)

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
            tiff_fn = os.path.splitext(os.path.basename(trialTable['fnRegDSInt'][DMDix,trialIx][0]))[0]
            if not os.path.exists(os.path.join(dr, tiff_fn + '.tif')):
                print(f'Missing tiff file: {tiff_fn}')
                keepTrials[DMDix,trialIx] = False
                
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

    # Create filename with timestamp
    fnsave = os.path.join(savedr, f'IntegrationSummary-{datetime.datetime.now().strftime("%y%m%d-%H%M%S")}.mat')

    # Load aData file
    aData = spio.loadmat(trialTable['fnAdataInt'][0,firstValidTrial][0])['aData'][0,0]

    # Extract metadata parameters
    params = {}
    params['discardInitial_s'] = 0.1
    params['numChannels'] = aData['numChannels'][0,0]
    params['alignHz'] = aData['alignHz'][0,0]
    params['analyzeHz'] = 100
    params['activityChannel'] = 1
    params['savedr'] = savedr
    params['decayTau_s'] = 0.05

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


    # Save as multi-channel TIFF
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

        data_file = os.path.join(params['savedr'], f'lowres_data_DMD{DMDix+1}.npz')
        
        if os.path.exists(data_file):
            print(f'Loading existing low resolution data from {data_file}')
            data_arrays = np.load(data_file)
            lowResData = data_arrays['lowResData']
            lowResDataCt = data_arrays['lowResDataCt'] 
            lowResMotionR = data_arrays['lowResMotionR']
            lowResMotionC = data_arrays['lowResMotionC']
            lowResMotionZ = data_arrays['lowResMotionZ']
        else:
            with mp.Pool(processes=min(mp.cpu_count(),len(trial_info))) as pool:
                results = list(pool.imap(get_trial_data_partial, trial_info))

            lowResData = np.concatenate([r[0] for r in results], axis=1)
            lowResDataCt = np.concatenate([r[1] for r in results], axis=1)
            # lowResBaseline = np.concatenate([r[2] for r in results], axis=1)
            lowResMotionR = np.concatenate([r[2]['motionDSr'] for r in results], axis=0)
            lowResMotionC = np.concatenate([r[2]['motionDSc'] for r in results], axis=0)
            lowResMotionZ = np.concatenate([r[2]['motionDSz'] for r in results], axis=0)
            # lowResBrightness = np.concatenate([r[2]['brightnessDS'] for r in results], axis=0)

            # Save data arrays
            data_arrays = {
                'lowResData': lowResData,
                'lowResDataCt': lowResDataCt,
                'lowResMotionR': lowResMotionR,
                'lowResMotionC': lowResMotionC,
                'lowResMotionZ': lowResMotionZ
            }
            np.savez(data_file, **data_arrays)
            print(f'Saved low resolution data to {data_file}')

        uniqueMotion, motInds = np.unique(np.round(np.concatenate((lowResMotionR,lowResMotionC,lowResMotionZ),axis=1)),axis=0,return_inverse=True)
        motIndsToKeep = (np.bincount(motInds) > 100).nonzero()[0]

        uniqueMotionYX = np.unique(uniqueMotion[motIndsToKeep,:2],axis=0)

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

        # print(uniqueMotion.shape)
        # print(np.histogram(motInds,bins=np.arange(0,np.max(motInds)+1)))
        # import napari
        # data_viewer = napari.view_image(lowResData)
        # data_viewer.add_image(lowResBaseline * lowResDataCt)


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

        data_for_nmf = lowResData / lowResDataCt
        decayTau_frames = params['decayTau_s']*params['alignHz']
        decay_kernel = np.exp(np.linspace(-np.ceil(decayTau_frames*3),0,int(np.ceil(decayTau_frames*3)+1))/decayTau_frames)
        data_for_nmf = signal.convolve2d(data_for_nmf,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same')
        data_for_nmf = torch.from_numpy(data_for_nmf.astype(np.float32))

        background_spatial_components = torch.zeros((numSuperPixels,len(motIndsToKeep)))
        for i in range(len(motIndsToKeep)):
            motion_idx = motIndsToKeep[i]
            motion_frames = (motInds == motion_idx).nonzero()[0]
            background_spatial_components[:,i] = torch.median(data_for_nmf[:, motion_frames],dim=1)[0]
            background_spatial_components[:,i] = background_spatial_components[:,i] / torch.norm(background_spatial_components[:,i])

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

        residual = torch.zeros_like(data_for_nmf)
        residual[:] = torch.nan
        background = torch.zeros_like(data_for_nmf)
        background[:] = torch.nan
        rho = torch.zeros((lowResData.shape[1],len(selPixIdxs)))
        rho[:] = torch.nan
        for i, motion_idx in enumerate(motIndsToKeep):
            motion_frames = (motInds == motion_idx).nonzero()[0]
            
            # Get background component as contiguous vector
            bg_comp = background_spatial_components[:,i]
            
            # Compute background trace more efficiently
            bg_norm_sq = torch.sum(bg_comp * bg_comp)
            bg_trace = (bg_comp @ data_for_nmf[:,motion_frames]) / bg_norm_sq

            background[:,motion_frames] = torch.outer(bg_comp, bg_trace)
            
            # Compute residual in-place
            residual[:,motion_frames] = data_for_nmf[:,motion_frames] - background[:,motion_frames]
            
            # Use pre-computed HD matrix
            rho[motion_frames,:] = (HD_mots[i].T @ residual[:,motion_frames]).T

        rho = rho.numpy()

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
            # dilated_nan_mask = ndimage.binary_dilation(nan_mask, structure=dilation_struct)
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
            # valid_maxima = (batch_rho > 0) & local_maxima & ~dilated_nan_mask
            valid_maxima = local_maxima & ~dilated_nan_mask
            # Add squared values for valid time points
            start_offset = batch_start - padded_start
            end_offset = start_offset + (batch_end - batch_start)
            # rho_activity += np.sum(np.where(valid_maxima[start_offset:end_offset], 
            #                               batch_rho[start_offset:end_offset]**2, 0), 
            #                      axis=0)
            rho_activity += np.sum(np.where(valid_maxima[start_offset:end_offset], 
                                          batch_rho[start_offset:end_offset]**3, 0), 
                                 axis=0)
        
        

        # lam = 0.05 * np.nanpercentile(rho,99,axis=0)
        # rho_clip = np.clip(rho - lam, 0, None)

        # rho_energy = np.nansum(rho_clip**2,axis=0)

        # energy_im = np.zeros((dmdPixelsPerColumn*dmdPixelsPerRow,))
        # energy_im[selPixIdxs] = rho_energy
        # energy_im = energy_im.reshape((dmdPixelsPerColumn,dmdPixelsPerRow))
        # energy_im[energy_im == 0] = np.nan

        # print("Computing skewness without intermediate tensors...")
        
        # rho_centered = rho - np.nanmean(rho,axis=0,keepdims=True)
        # cubed_values = rho_centered**3
        
        # # Sum along the axis
        # rho_skew = np.nansum(cubed_values, axis=0)
        
        # # Clean up the intermediate tensor
        # del cubed_values

        # print(f"Skewness computed, shape: {rho_skew.shape}")

        # skew_im = np.zeros((dmdPixelsPerColumn*dmdPixelsPerRow,))
        # skew_im[selPixIdxs] = rho_skew
        # skew_im = skew_im.reshape((dmdPixelsPerColumn,dmdPixelsPerRow))
        # skew_im[skew_im == 0] = np.nan

        act_im = rho_activity.copy()
        nan_mask = (act_im == 0) | (np.isnan(act_im))
        act_im[nan_mask] = np.nan

        med_act_im = ndimage.generic_filter(act_im, np.nanmedian, size=(15,15))
        act_im = act_im - med_act_im
        act_im[nan_mask] = np.nan
        # act_im[nan_mask] = np.median(act_im[~nan_mask])

        
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

        # def compute_topological_prominence(act_im, maxima_coords, connectivity=8):
        #     # maxima_coords: either (rows, cols) tuple of arrays or Nx2 array
        #     if isinstance(maxima_coords, tuple):
        #         peak_rows, peak_cols = maxima_coords
        #         peaks_rc = np.column_stack([peak_rows, peak_cols])
        #     else:
        #         peaks_rc = np.asarray(maxima_coords)
        #         if peaks_rc.ndim != 2 or peaks_rc.shape[1] != 2:
        #             raise ValueError("maxima_coords must be (rows, cols) or Nx2 array")

        #     H, W = act_im.shape
        #     N = H * W

        #     img = np.asarray(act_im, dtype=float).copy()
        #     img[~np.isfinite(img)] = -np.inf  # exclude NaNs

        #     # Map each peak to a unique id
        #     peak_ids = -np.ones((H, W), dtype=int)
        #     for pid, (r, c) in enumerate(peaks_rc):
        #         peak_ids[r, c] = pid

        #     # Union-Find structures
        #     parent = np.arange(N, dtype=np.int32)
        #     rank = np.zeros(N, dtype=np.int8)
        #     active = np.zeros(N, dtype=bool)

        #     comp_peak_id = -np.ones(N, dtype=np.int32)        # valid only at roots
        #     comp_peak_height = np.full(N, -np.inf, dtype=float)
        #     prominences = np.full(len(peaks_rc), np.nan, dtype=float)

        #     # Helpers
        #     def idx_rc_to_flat(r, c):
        #         return r * W + c

        #     def find(x):
        #         while parent[x] != x:
        #             parent[x] = parent[parent[x]]
        #             x = parent[x]
        #         return x

        #     def union(a, b, level):
        #         ra, rb = find(a), find(b)
        #         if ra == rb:
        #             return ra

        #         # union by rank
        #         if rank[ra] < rank[rb]:
        #             ra, rb = rb, ra
        #         parent[rb] = ra
        #         if rank[ra] == rank[rb]:
        #             rank[ra] += 1

        #         ida, idb = comp_peak_id[ra], comp_peak_id[rb]
        #         ha, hb = comp_peak_height[ra], comp_peak_height[rb]

        #         if ida != -1 and idb != -1:
        #             # both components have peaks; lower peak dies at saddle "level"
        #             if ha >= hb:
        #                 prominences[idb] = hb - level
        #                 comp_peak_id[ra], comp_peak_height[ra] = ida, ha
        #             else:
        #                 prominences[ida] = ha - level
        #                 comp_peak_id[ra], comp_peak_height[ra] = idb, hb
        #         elif ida != -1:
        #             comp_peak_id[ra], comp_peak_height[ra] = ida, ha
        #         else:
        #             comp_peak_id[ra], comp_peak_height[ra] = idb, hb

        #         return ra

        #     # Neighbor offsets
        #     if connectivity == 4:
        #         nbrs = [(-1,0), (1,0), (0,-1), (0,1)]
        #     else:
        #         nbrs = [(-1,-1),(-1,0),(-1,1),(0,-1),(0,1),(1,-1),(1,0),(1,1)]

        #     vals = img.ravel()
        #     order = np.flatnonzero(np.isfinite(vals))
        #     order = order[np.argsort(vals[order])[::-1]]  # process high to low

        #     # Fast row/col access
        #     rows = np.arange(H)[:, None] * np.ones(W, dtype=int)
        #     cols = np.ones((H, 1), dtype=int) * np.arange(W)
        #     rows = rows.ravel()
        #     cols = cols.ravel()

        #     for p in order:
        #         r, c = rows[p], cols[p]
        #         if vals[p] == -np.inf:
        #             continue

        #         active[p] = True
        #         parent[p] = p
        #         rank[p] = 0

        #         # seed peak if this pixel is a detected maximum
        #         pid = peak_ids[r, c]
        #         if pid != -1:
        #             comp_peak_id[p] = pid
        #             comp_peak_height[p] = vals[p]
        #         else:
        #             comp_peak_id[p] = -1
        #             comp_peak_height[p] = -np.inf

        #         # merge with active neighbors
        #         for dr, dc in nbrs:
        #             rr, cc = r + dr, c + dc
        #             if 0 <= rr < H and 0 <= cc < W:
        #                 q = idx_rc_to_flat(rr, cc)
        #                 if active[q]:
        #                     union(p, q, level=vals[p])

        #     # Peaks that never merged into a higher peak: use base level
        #     base = np.nanmin(act_im)
        #     for pid, (r, c) in enumerate(peaks_rc):
        #         if np.isnan(prominences[pid]):
        #             prominences[pid] = act_im[r, c] - base

        #     return prominences
        
        # def peak_widths_at_half_prominence(act_im, maxima_coords, prominence, connectivity=8):
        #     # maxima_coords: (rows, cols) tuple or Nx2 array
        #     if isinstance(maxima_coords, tuple):
        #         rows, cols = maxima_coords
        #         peaks = np.column_stack([rows, cols])
        #     else:
        #         peaks = np.asarray(maxima_coords)

        #     H, W = act_im.shape
        #     act = np.asarray(act_im, float)
        #     act[~np.isfinite(act)] = np.nan

        #     peak_vals = act[peaks[:,0], peaks[:,1]]
        #     saddle = peak_vals - prominence
        #     half_level = saddle + 0.5 * prominence

        #     # 8-conn or 4-conn structuring element
        #     structure = np.ones((3,3), bool) if connectivity == 8 else np.array([[0,1,0],[1,1,1],[0,1,0]], bool)

        #     widths = np.full(len(peaks), np.nan)
        #     areas = np.full(len(peaks), np.nan)

        #     finite_mask = np.isfinite(act)
        #     for i, (r, c) in enumerate(peaks):
        #         thr = half_level[i]
        #         if not np.isfinite(thr):
        #             continue

        #         mask = (act >= thr) & finite_mask
        #         labels, nlab = ndimage.label(mask, structure=structure)
        #         lab = labels[r, c]
        #         if lab == 0:
        #             continue

        #         comp_area = np.count_nonzero(labels == lab)
        #         areas[i] = comp_area
        #         # Equivalent diameter (2 * sqrt(area/pi)) as a width estimate
        #         widths[i] = 2.0 * np.sqrt(comp_area / np.pi)

        #     return widths, areas, half_level
            
        # # Calculate prominences
        # peak_prominences = compute_topological_prominence(act_im, maxima_coords)
        # widths, areas, _ = peak_widths_at_half_prominence(act_im,maxima_coords,peak_prominences)

        mad = np.median(np.abs(act_im[~np.isnan(act_im)] - np.median(act_im[~np.isnan(act_im)])))
        thresh = 3 * 1.4826 * mad
        # thresh = 1.75 * maxima_values[min(len(maxima_values)-1,int(numSuperPixels*9 * 0.01))]

        # # Create a mask of nan values
        # nan_mask = np.isnan(act_im)

        # # Dilate the nan mask to cover 3 pixel radius
        # dilated_nan_mask = ndimage.binary_dilation(nan_mask, structure=np.ones((7,7)))

        # Calculate threshold for top 5% of values
        # top_5_threshold = np.nanpercentile(act_im, 95)

        # Remove maxima that are within 3 pixels of nan values, unless they are in top 5%
        # valid_maxima = (~dilated_nan_mask[maxima_coords]) | (maxima_values > top_5_threshold)
        valid_maxima = maxima_values > thresh
        # valid_maxima = peak_prominences > thresh
        maxima_coords = (maxima_coords[0][valid_maxima], maxima_coords[1][valid_maxima])
        maxima_values = maxima_values[valid_maxima]
        source_seeds = np.array(maxima_coords).T

        nSources = source_seeds.shape[0]

        source_params = torch.cat([
            torch.tensor(source_seeds, dtype=torch.float32),
            torch.ones(nSources, 2, dtype=torch.float32),
            torch.zeros(nSources, 1, dtype=torch.float32)
        ], dim=1)

        def sel_pix_gaussian_2d(gaussian_params):
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

            normalization_factor = 1 / (2 * np.pi * y_sigmas * x_sigmas * torch.sqrt(1 - corr_coef**2))
            
            return normalization_factor * torch.exp(exponent)  # Shape: [nPixels, nSources]

        A = torch.zeros((nPixels, nSources))
        A[selPixIdxs,:] = sel_pix_gaussian_2d(source_params * torch.tensor([1, 1, 1, 1, 1]))

        # NMF parameters
        als_nmf_iters = 4
        mult_nmf_max_iters = 10
        outer_loop_iters = 3
        nmf_tol = 1e-6
        # Gaussian fitting optimization parameters
        learning_rate = 0.1
        num_epochs = 200
        gd_tol = 1e-8
        
        phi_lowRes = torch.zeros(lowResData.shape[1], nSources+1, dtype=torch.float32)
        phi_lowRes[:] = np.nan

        X_mots = [None] * len(motIndsToKeep)
        
        for outer_loop_iter in range(outer_loop_iters):

            for i, motion_idx in enumerate(motIndsToKeep):
                # Extract data for the most common motion mode
                motion_frames = (motInds == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                X = X / torch.max(X, dim=0, keepdim=True)[0]
                X = torch.clamp(X, min=torch.exp(torch.tensor(-3.0))) - torch.exp(torch.tensor(-3.0))

                # newR = refR+uniqueMotion[motion_mode_idx,0].astype(int)
                # newC = refC+uniqueMotion[motion_mode_idx,1].astype(int)
                # newD = torch.from_numpy(fastZ2RefZ[f'DMD{DMDix+1}'][refD.numpy()].flatten()-1)+uniqueMotion[motion_mode_idx,2].astype(int)

                # add background spatial component
                # background_spatial_component = torch.as_tensor(refStack[f'DMD{DMDix+1}'][0][newD,newR,newC].reshape((-1,1)))
                X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # normalize spatial components
                X = X / torch.norm(X,dim=0,keepdim=True)
                X_support = X > 0

                # initialize phi by projected least squares
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]  # This gives all time points at once
                
                # Add small regularization to ensure stability
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                
                # Solve the system for all time points at once
                # We need to solve (X^T * X) * phi = X^T * data for each column of data
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
                
                # Ensure non-negativity
                phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:], min=0)
            
                # Initialize variables for NMF
                error_values = []
                prev_reconstruction_error = float('inf')

                for iter_idx in tqdm(range(mult_nmf_max_iters),desc=f'Multiplicative NMF {i+1}/{len(motIndsToKeep)}'):
                    # Multiplicative update for NMF
                    # Update phi (temporal components) using multiplicative update rule
                    # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
                    numerator = X.T @ data_for_nmf[:, motion_frames]
                    denominator = (X.T @ X) @ phi_lowRes[motion_frames,:].T + 1e-10
                    phi_lowRes[motion_frames,:] = phi_lowRes[motion_frames,:] * (numerator / denominator).T
                    
                    # Update X (spatial components) using multiplicative update rule
                    # X = X * (data * phi^T) / (X * phi * phi^T + epsilon)
                    numerator = data_for_nmf[:, motion_frames] @ phi_lowRes[motion_frames,:]
                    denominator = X @ (phi_lowRes[motion_frames,:].T @ phi_lowRes[motion_frames,:]) + 1e-10
                    X = X * (numerator / denominator)

                    # Apply spatial support constraint
                    X[~X_support] = 0
                    X[:,-1] = background_spatial_components[:,i]
                    
                    # Normalize X and phi to avoid scaling ambiguity
                    for s in range(nSources):
                        norm = torch.norm(X[:, s])
                        if norm > 0:
                            X[:, s] = X[:, s] / norm
                            phi_lowRes[motion_frames, s] = phi_lowRes[motion_frames, s] * norm
                    
                    # Calculate reconstruction error
                    reconstruction = X @ phi_lowRes[motion_frames,:].T
                    current_error = torch.mean((data_for_nmf[:, motion_frames] - reconstruction)**2).item()
                    error_values.append(current_error)
                    
                    # # Check convergence
                    # if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    #     break
                        
                    prev_reconstruction_error = current_error
                
                X_mots[i] = X

            # # Plot reconstruction error at the end
            # plt.figure(figsize=(10, 5))
            # plt.title('NMF Reconstruction Error')
            # plt.xlabel('Iteration')
            # plt.ylabel('Mean Squared Error')
            # plt.grid(True)
            # plt.plot(range(1, len(error_values) + 1), error_values, 'b-', marker='o')
            # plt.tight_layout()


            # Fit A_final for X=HA where each column of A_final is a Gaussian using gradient descent
            print("Fitting Gaussian spatial components using gradient descent...")
            
            # Make parameters require gradients
            optim_params = source_params.clone().requires_grad_(True)
            
            # Initialize Adam optimizer
            optimizer = torch.optim.Adam([optim_params], lr=learning_rate)
            
            # Optimization loop
            losses = []
            for epoch in tqdm(range(num_epochs),desc='Gaussian fitting'):
                # Zero gradients
                optimizer.zero_grad()
                
                # Compute all Gaussian values at once
                # all_gaussian_values = amplitudes * torch.exp(y_term + x_term)  # Shape: [nPixels, nSources]
                all_gaussian_values = sel_pix_gaussian_2d(optim_params * torch.tensor([1, 1, 1, 1, 1]))

                total_loss = 0
                
                i = epoch % len(motIndsToKeep)
                motion_idx = motIndsToKeep[i]

                # Compute sparse matrix multiplication for all sources at once
                total_fit = torch.sparse.mm(H_mots[i], all_gaussian_values)

                
                total_fit = total_fit / torch.max(total_fit, dim=0, keepdim=True)[0]
                total_fit = torch.clamp(total_fit, min=torch.exp(torch.tensor(-3.0))) - torch.exp(torch.tensor(-3.0))
                norms = torch.norm(total_fit, dim=0, keepdim=True)
                total_fit = torch.where(norms > 0, total_fit / norms, total_fit)

                # total_fit[total_fit < torch.max(total_fit,dim=0,keepdim=True)[0] * torch.exp(torch.tensor(-3.0))] = 0
            
                # Compute loss across all sources
                support = torch.logical_or(X_mots[i][:,:-1] > 0, total_fit > 0)
                loss = torch.mean((total_fit[support] - X_mots[i][:,:-1][support]) ** 2)
                
                loss.backward()
                optimizer.step()
                
                # Apply constraints after optimizer step
                with torch.no_grad():
                    # Ensure amplitude and sigmas are positive
                    # optim_params[:,0].clamp_(min=1e-5)
                    optim_params[:,2].clamp_(min=0.3, max=5)
                    optim_params[:,3].clamp_(min=0.3, max=5)

                if i == len(motIndsToKeep)-1:
                    total_loss = 0
                    for i, motion_idx in enumerate(motIndsToKeep):
                        total_fit = torch.sparse.mm(H_mots[i], all_gaussian_values)

                        total_fit = total_fit / torch.max(total_fit, dim=0, keepdim=True)[0]
                        total_fit = torch.clamp(total_fit, min=torch.exp(torch.tensor(-3.0))) - torch.exp(torch.tensor(-3.0))
                        norms = torch.norm(total_fit, dim=0, keepdim=True)
                        total_fit = torch.where(norms > 0, total_fit / norms, total_fit)

                        # total_fit[total_fit < torch.max(total_fit,dim=0,keepdim=True)[0] * torch.exp(torch.tensor(-3.0))] = 0

                        # Compute loss across all sources
                        support = torch.logical_or(X_mots[i][:,:-1] > 0, total_fit > 0)
                        total_loss += torch.mean((total_fit[support] - X_mots[i][:,:-1][support]) ** 2)

                    # Check for convergence
                    losses.append(total_loss.item())
                    if epoch // len(motIndsToKeep) > 0 and abs(losses[-1] - losses[-2]) < gd_tol:
                        print(f"Converged at epoch {epoch+1} with loss: {loss.item():.8f}")
                        break

                # with torch.no_grad():
                #     optim_params[norms.squeeze() == 0,2] = optim_params[norms.squeeze() == 0,2] * 5
                #     optim_params[norms.squeeze() == 0,3] = optim_params[norms.squeeze() == 0,3] * 2
                
                # Print progress
                # if (epoch + 1) % 10 == 0:
                #     print(f"Epoch {epoch+1}/{num_epochs}, Loss: {loss.item():.8f}")
            
            # Update parameters
            source_params = optim_params.detach()

            A = torch.zeros((nPixels, nSources), dtype=torch.float32)
            # A_final[selPixIdxs, :] = (amplitudes * torch.exp(y_term + x_term))
            A[selPixIdxs, :] = sel_pix_gaussian_2d(source_params * torch.tensor([1, 1, 1, 1, 1]))
            
            print("Gaussian fitting complete")

            # X_update = torch.sparse.mm(H, A[selPixIdxs,:])
            # X_update = torch.clamp(X_update, min=1e-4) - 1e-4
            # norms = torch.norm(X_update, dim=0, keepdim=True)
            # X_update = torch.where(norms > 0, X_update / norms, X_update)

            # # Remove sources with all zeros
            # non_zero_sources = torch.any(X_update > 0, dim=0)
            # X_update = X_update[:, non_zero_sources]
            # source_params = source_params[non_zero_sources]
            # nSources = torch.sum(non_zero_sources).item()

            # X_tensor = torch.cat((X_update, background_spatial_component.unsqueeze(-1)),dim=1)
            # phi_lowRes = torch.cat((phi_lowRes[:,non_zero_sources.nonzero()[:,0]], phi_lowRes[:,-1].unsqueeze(-1)),dim=1)

            for i, motion_idx in enumerate(motIndsToKeep):
                # Extract data for the most common motion mode
                motion_frames = (motInds == motion_idx).nonzero()[0]

                # project image space (A) into superpixel space (X)
                X = torch.sparse.mm(H_mots[i], A[selPixIdxs,:])
                thresh = torch.max(X,dim=0,keepdim=True)[0] * torch.exp(torch.tensor(-3.0))
                X = torch.clamp(X, min=thresh) - thresh
                # X[~X_support] = 0
                # X = X / torch.max(X, dim=0, keepdim=True)[0]
                # X = torch.clamp(X, min=torch.exp(torch.tensor(-3.0))) - torch.exp(torch.tensor(-3.0))

                # newR = refR+uniqueMotion[motion_mode_idx,0].astype(int)
                # newC = refC+uniqueMotion[motion_mode_idx,1].astype(int)
                # newD = torch.from_numpy(fastZ2RefZ[f'DMD{DMDix+1}'][refD.numpy()].flatten()-1)+uniqueMotion[motion_mode_idx,2].astype(int)

                # add background spatial component
                # background_spatial_component = torch.as_tensor(refStack[f'DMD{DMDix+1}'][0][newD,newR,newC].reshape((-1,1)))
                X = torch.concat((X,background_spatial_components[:,i].unsqueeze(-1)),dim=1)

                # normalize spatial components
                norms = torch.norm(X,dim=0,keepdim=True)
                if torch.any(norms == 0):
                    print(f"Warning: {np.nonzero(norms.numpy().squeeze() == 0)[0]} norms are zero for motion {i}")
                # X = torch.where(norms > 0, X / norms, X)
                # X_support = X > 0

                # initialize phi by projected least squares
                XtX = X.T @ X
                Xtd = X.T @ data_for_nmf[:, motion_frames]  # This gives all time points at once
                
                # Add small regularization to ensure stability
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                
                # Solve the system for all time points at once
                # We need to solve (X^T * X) * phi = X^T * data for each column of data
                phi_lowRes[motion_frames,:] = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T

            # phi_lowRes = phi_lowRes.clamp(min=0)
            
            # sort sources by variance
            sortorder = np.argsort(-np.nansum((phi_lowRes[:,:nSources].numpy()-np.nanmean(phi_lowRes[:,:nSources].numpy(),axis=0))**2,axis=0))
            source_params = source_params[sortorder,:]
            A = A[:,sortorder]
            sortorder = np.append(sortorder,max(sortorder)+1)
            phi_lowRes = phi_lowRes[:,sortorder]

            # distances = []
            
            # for i in range(nSources):
            #     for j in range(i+1, nSources):
            #         dist = np.sqrt(np.sum((source_params[i,:2].numpy() - source_params[j,:2].numpy())**2))
            #         distances.append(dist)
            
            # # Plot histogram of distances
            # plt.figure()
            # plt.hist(distances, bins=50)
            # plt.xlabel('Distance between source pairs (pixels)')
            # plt.ylabel('Count')
            # plt.title('Distribution of pairwise distances between sources')
            # plt.show()

            if (outer_loop_iter+1) % 4 == 0:
                residual = data_tensor - X_tensor @ phi_lowRes.T
                SNR = torch.zeros(nSources)
                for i in range(nSources-1,-1,-1):
                    sorted_vals = torch.sort(X_tensor[:,i], descending=True)[0]
                    cumsum_vals = torch.cumsum(sorted_vals, dim=0)
                    total_mass = cumsum_vals[-1]
                    half_mass_idx = torch.searchsorted(cumsum_vals, total_mass * 0.5)
                    thresh = sorted_vals[half_mass_idx]
                    validPixs = (X_tensor[:,i] >= thresh).nonzero()[:,0]
                    varExpSource = torch.sum((X_tensor[validPixs,i].unsqueeze(-1) * phi_lowRes[:,i].unsqueeze(-1).T)**2)
                    varResidual = torch.sum(residual[validPixs]**2)
                    SNR[i] = varExpSource / varResidual
                
                SNR_cut = 0.5
                keepSources = (SNR > SNR_cut).nonzero()[:,0]
                nSources = keepSources.shape[0]
                source_params = source_params[keepSources,:]
                keepSources = np.append(keepSources,phi_lowRes.shape[1]-1)
                phi_lowRes = phi_lowRes[:,keepSources]
                X_tensor = X_tensor[:,keepSources]

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
                                              dr=dr, trialTable=trialTable, A_final=A, psf=psf)

        with mp.Pool(processes=min(mp.cpu_count(),len(trial_info))) as pool:
            results = list(pool.imap(get_high_res_traces_partial, trial_info))
        
        # Save all peak data to a single file
        # source_extraction_data = {
        #     'phi': [r[0] for r in results],
        #     'timePts': [r[1] for r in results],
        #     'source_params': source_params,
        #     'A_final': A_final,
        #     'selPixIdxs': [r[2] for r in results],
        # }
        # output_filename = os.path.join(params['savedr'], f'source_extraction_data_DMD{DMDix+1}.npz')
        # np.savez(output_filename, **source_extraction_data)
        # print(f"Saved source extraction data to {output_filename}")

        # Save data to HDF5 file
        output_h5_filename = os.path.join(params['savedr'], f'experiment_summary_nonmax_supp_updated_DoG.h5')
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
            temporal_group = source_group.create_group('temporal')

            spatial_group.create_dataset('source_params', data=source_params.numpy())
            spatial_group.create_dataset('fp_masks', data=A.numpy().reshape(dmdPixelsPerRow,dmdPixelsPerColumn,-1).transpose(2,0,1))
            spatial_group.create_dataset('fp_coords', data=source_params[:,:2].numpy())

            dF = np.concatenate([r[0] for r in results], axis=0)
            temporal_group.create_dataset('dF', data=dF)

            trial_start_idxs = np.concatenate([[0], np.cumsum([len(r[1]) for r in results])[:-1]])

            frame_group = dmd_group.create_group('frame_info')
            frame_group.create_dataset('trial_start_idxs', data=trial_start_idxs)
            frame_group.create_dataset('discard_frames', data=np.any(np.isnan(dF), axis=1))
            # frame_group.create_dataset('selPixIdxs', data=[r[2] for r in results])

            globalF = np.concatenate([r[3] for r in results], axis=0)
            global_group = dmd_group.create_group('global')
            global_group.create_dataset('F', data=globalF)

            visualizations = dmd_group.create_group('visualizations')
            visualizations.create_dataset('act_im', data=act_im)

        print(f"Added DMD{DMDix+1} data to {output_h5_filename}")

if __name__ == '__main__':
    main()

# for iter_idx in tqdm(range(als_nmf_iters),desc='Alternating least squares NMF'):
#     # Update phi (temporal components) with X fixed
#     # Formulate as least squares: min ||data_t - X * phi_t||^2 for all t
#     XtX = X_tensor.T @ X_tensor
#     Xtd = X_tensor.T @ data_tensor  # This gives all time points at once

#     # Add small regularization to ensure stability
#     regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])

#     # Solve the system for all time points at once
#     # We need to solve (X^T * X) * phi = X^T * data for each column of data
#     phi_lowRes = torch.linalg.solve(
#         regularized_XtX,
#         Xtd
#     ).T

#     # Ensure non-negativity
#     phi_lowRes = torch.clamp(phi_lowRes, min=0)

#     # Update X (spatial components) with phi fixed
#     # Formulate as least squares: min ||data - X * phi^T||^2 for all sources
#     phi_phi_t = phi_lowRes.T @ phi_lowRes  # Shape: [sources, sources]
#     data_phi = data_tensor @ phi_lowRes  # Shape: [pixels, sources]

#     regularized_phi_phi_t = phi_phi_t + 1e-10 * torch.eye(phi_phi_t.shape[0])

#     X_update = torch.linalg.solve(regularized_phi_phi_t, data_phi.T).T

#     # Ensure non-negativity
#     X_tensor = torch.clamp(X_update, min=0)
#     X_tensor[~X_support] = 0
#     X_tensor[:,-1] = background_spatial_component

#     # Normalize X and phi to avoid scaling ambiguity
#     for s in range(nSources):
#         norm = torch.norm(X_tensor[:, s])
#         if norm > 0:
#             X_tensor[:, s] = X_tensor[:, s] / norm
#             phi_lowRes[:, s] = phi_lowRes[:, s] * norm

#     # Calculate reconstruction error
#     reconstruction = X_tensor @ phi_lowRes.T
#     current_error = torch.mean((data_tensor - reconstruction)**2).item()
#     error_values.append(current_error)

#     # Check convergence
#     if abs(prev_reconstruction_error - current_error) < nmf_tol:
#         break
    
#     prev_reconstruction_error = current_error


            
            # # Initialize parameters based on the current X_tensor
            # for s in range(nSources):
            #     # Get the current spatial component
            #     source_weights = X_tensor[:, s].numpy()
            #     source_weights_filtered = convMatrix @ source_weights
                
            #     # Find the maximum value and its position
            #     max_idx = np.argmax(source_weights_filtered)
            #     # max_val = source_weights[max_idx]

            #     frame = np.zeros((800,1280))
            #     frame[:] = np.nan
            #     frame[refR + uniqueMotion[motion_mode_idx,0].astype(int), refC + uniqueMotion[motion_mode_idx,1].astype(int)] = source_weights
                
            #     patch_size = 21  # 5x5 patch around maximum
            #     max_r = refR[max_idx] + uniqueMotion[motion_mode_idx,0].astype(int)
            #     max_c = refC[max_idx] + uniqueMotion[motion_mode_idx,1].astype(int)
                
            #     # Define patch boundaries
            #     r_start = max(0, max_r - patch_size//2)
            #     r_end = min(800, max_r + patch_size//2 + 1)
            #     c_start = max(0, max_c - patch_size//2)
            #     c_end = min(1280, max_c + patch_size//2 + 1)
                
            #     # Extract patch from frame
            #     patch = frame[r_start:r_end, c_start:c_end]
            #     del frame
                
            #     # Calculate centroid if patch has valid values
            #     if not np.all(np.isnan(patch)):
            #         # Create coordinate grids for the patch
            #         r_coords, c_coords = np.meshgrid(
            #             np.arange(r_start, r_end),
            #             np.arange(c_start, c_end),
            #             indexing='ij'
            #         )
                    
            #         # Get valid (non-NaN) values
            #         valid_mask = ~np.isnan(patch)
            #         if np.any(valid_mask):
            #             weights = patch[valid_mask]
            #             weights = np.maximum(weights, 0)  # Ensure non-negative weights
                        
            #             if np.sum(weights) > 0:
            #                 # Calculate weighted centroid
            #                 centroid_r = np.sum(r_coords[valid_mask] * weights) / np.sum(weights)
            #                 centroid_c = np.sum(c_coords[valid_mask] * weights) / np.sum(weights)
            #             else:
            #                 print("No valid values in patch")
            #                 centroid_r = max_r
            #                 centroid_c = max_c
            #         else:
            #             print("No valid values in patch")
            #             centroid_r = max_r
            #             centroid_c = max_c


            #     # Use the position of maximum value as initial mean
            #     # source_params[s, 0] = torch.tensor(max_val) # amplitude
            #     source_params[s, 0] = centroid_r
            #     source_params[s, 1] = centroid_c
            #     # source_params[s,2:] = source_params[s,2:]
            
            # Optimization loop for all sources together
            # print("Fitting all Gaussian sources simultaneously")

            

        # phi_lowRes = torch.zeros(lowResData.shape[1], nSources+1, dtype=torch.float32)
        # phi_lowRes[:] = np.nan
        
        # nonneg = False

        # for motion_idx in motIndsToKeep:

        #     # Extract data for the most common motion mode
        #     motion_frames = (motInds == motion_idx).nonzero()[0]
            
        #     data_tensor = torch.from_numpy(data_for_nmf[:, motion_frames].astype(np.float32))

        #     sparseHIndsShifted = sparseHInds.copy()
        #     sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)
        #     sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        #     sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        #     H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

        #     # project image space (A) into superpixel space (X)
        #     X = torch.sparse.mm(H, A_final[selPixIdxs,:])
        #     X = torch.clamp(X, min=1e-4) - 1e-4

        #     # newR = refR+uniqueMotion[motion_idx,0].astype(int)
        #     # newC = refC+uniqueMotion[motion_idx,1].astype(int)
        #     # newD = torch.from_numpy(fastZ2RefZ[f'DMD{DMDix+1}'][refD.numpy()].flatten()-1)+uniqueMotion[motion_idx,2].astype(int)

        #     # add background spatial component
        #     # background_spatial_component = torch.as_tensor(refStack[f'DMD{DMDix+1}'][0][newD,newR,newC].reshape((-1,1)))
        #     background_spatial_component = torch.median(data_tensor,dim=1)[0]
        #     background_spatial_component = background_spatial_component / torch.norm(background_spatial_component)
        #     X = torch.concat((X,background_spatial_component.unsqueeze(-1)),dim=1)

        #     XtX = X.T @ X
        #     Xtd = X.T @ data_tensor  # This gives all time points at once
            
        #     # Add small regularization to ensure stability
        #     regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
            
        #     # Solve the system for all time points at once
        #     # We need to solve (X^T * X) * phi = X^T * data for each column of data
        #     phi_lowRes[motion_frames,:] = torch.linalg.solve(
        #         regularized_XtX,
        #         Xtd
        #     ).T
            
        #     if nonneg:
        #         # Ensure non-negativity
        #         phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:], min=0)
                    
        #         # # Normalize X and phi to avoid scaling ambiguity
        #         # for s in range(nSources):
        #         #     norm = torch.norm(X_tensor[:, s])
        #         #     if norm > 0:
        #         #         X_tensor[:, s] = X_tensor[:, s] / norm
        #         #         phi_lowRes[motion_frames, s] = phi_lowRes[motion_frames, s] * norm
                
        #         error_values = []
        #         prev_reconstruction_error = float('inf')
        #         for iter_idx in tqdm(range(mult_nmf_max_iters),desc='Multiplicative NMF'):
        #             # Multiplicative update for NMF
        #             # Update phi (temporal components) using multiplicative update rule
        #             # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
        #             numerator = X.T @ data_tensor
        #             denominator = (X.T @ X) @ phi_lowRes[motion_frames,:].T + 1e-10
        #             phi_update = phi_lowRes[motion_frames,:] * (numerator / denominator).T
        #             phi_lowRes[motion_frames,:] = phi_update
                    
        #             # # Normalize X and phi to avoid scaling ambiguity
        #             # for s in range(nSources):
        #             #     norm = torch.norm(X_tensor[:, s])
        #             #     if norm > 0:
        #             #         X_tensor[:, s] = X_tensor[:, s] / norm
        #             #         phi_lowRes[motion_frames, s] = phi_lowRes[motion_frames, s] * norm
                    
        #             # Calculate reconstruction error
        #             reconstruction = X @ phi_lowRes[motion_frames,:].T
        #             current_error = torch.mean((data_tensor - reconstruction)**2).item()
        #             error_values.append(current_error)
                    
        #             # Check convergence
        #             if abs(prev_reconstruction_error - current_error) < nmf_tol:
        #                 break
                        
        #             prev_reconstruction_error = current_error