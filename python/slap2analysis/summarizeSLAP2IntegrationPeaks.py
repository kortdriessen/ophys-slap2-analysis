# import all dependencies
import os.path
import sys
import time
import datetime
from pathlib import Path
from tkinter import filedialog

import numpy as np
import matplotlib.pyplot as plt
import torch
import scipy.io as spio
from scipy import sparse, signal, stats
from scipy.optimize import curve_fit
import h5py
import napari
import importlib

from skimage import io as skimio
from skimage import draw, measure
from sklearn.cluster import DBSCAN, OPTICS, HDBSCAN

from slap2_utils import DataFile

sys.path.append(str(Path(__file__).parent.parent))
import reconstruct

sys.path.append(str(Path(__file__).parent))
import deconvlucy

import multiprocessing as mp
from functools import partial

from tqdm import tqdm

import warnings

from qtpy.QtWidgets import QPushButton, QApplication, QMessageBox

import plotly.graph_objects as go

def get_trial_data(trial_info, DMDix, params, refStack, subsampleMatrixInds, sparseHInds, sparseHVals, 
                allSuperPixelIDs, dr, trialTable):
    trialIx, keepTrial = trial_info
    
    if not keepTrial:
        return [], [], []

    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]

    nPixels = dmdPixelsPerColumn * dmdPixelsPerRow

    aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]

    numCycles = aData['motionDSr'].shape[0]
    uniqueMotion, motInds = np.unique(np.round(np.concatenate((aData['motionDSr'],aData['motionDSc'],aData['motionDSz']+11),axis=1)),axis=0,return_inverse=True)

    Afinal = torch.empty((nPixels,0))
    phiFinal = torch.empty((numCycles,0))

    importlib.reload(reconstruct)

    startTime = time.time()
    baseline = reconstruct.reconstruct(Afinal,phiFinal,refStack[f'DMD{DMDix+1}'][0],torch.from_numpy(aData['brightnessDS']),subsampleMatrixInds.T,sparseHInds,sparseHVals,uniqueMotion,motInds).detach().numpy()
    endTime = time.time()

    # print(f"Elapsed Time: {endTime - startTime} sec")

    # plt.imshow(baseline)
    # plt.colorbar()
    # plt.show()

    source_fn = trialTable['filename'][DMDix,trialIx][0]
    firstLine = trialTable['firstLine'][DMDix,trialIx]
    lastLine = trialTable['lastLine'][DMDix,trialIx]

    hDataFile = DataFile(os.path.join(dr, source_fn))

    dt = 1/params['alignHz']/hDataFile.metaData.linePeriod_s

    DSframes = aData['DSframes'][0]
    nDSframes= len(DSframes)

    data = np.zeros((numSuperPixels,nDSframes))
    dataCt = np.zeros((numSuperPixels,nDSframes))
    
    # Pre-compute time windows for all frames
    timeWindows = [np.arange(max(1,np.floor(DSframes[i]-2*dt)), 
                            min(np.ceil(DSframes[i]+2*dt),hDataFile.numCycles*hDataFile.header['linesPerCycle'])+1) 
                    for i in range(nDSframes)]
    
    # Pre-compute line and cycle indices for all time windows
    lineIndices_all = [(tw - 1) % hDataFile.header['linesPerCycle'] + 1 for tw in timeWindows]
    cycleIndices_all = [np.floor((tw - 1) / hDataFile.header['linesPerCycle']) + 1 for tw in timeWindows]
    
    # Initialize timing variables
    start_time = time.time()
    
    for DSframeIx in range(nDSframes):
        if DSframeIx % 100 == 0:
            avg_time = (time.time() - start_time) / DSframeIx if DSframeIx > 0 else 0
            print(f"{DSframeIx} of {nDSframes}, Average time per frame: {avg_time:.3f} sec")
            
        lineIndices = lineIndices_all[DSframeIx]
        cycleIndices = cycleIndices_all[DSframeIx]
        
        allLineData = hDataFile.getLineData(lineIndices, cycleIndices, params['activityChannel'])
        
        # Vectorize line processing
        valid_lines = np.where([hDataFile.lineDataNumElements[int(li)-1] != 0 for li in lineIndices])[0]
        
        for lineIdx in valid_lines:
            positions = hDataFile.lineSuperPixelIDs[int(lineIndices[lineIdx])-1][0]
            zIdx = hDataFile.lineFastZIdxs[int(lineIndices[lineIdx])-1]
            
            # Compute lookup values and matches in one go
            lookup_values = positions * 100 + zIdx
            matching_mask = np.isin(allSuperPixelIDs[f'DMD{DMDix+1}'], lookup_values)
            matching_indices = np.where(matching_mask)[0]
            
            if len(matching_indices) > 0:
                # Create a mapping from lookup values to their indices
                value_to_pos = dict(zip(lookup_values.astype(np.uint32), range(len(lookup_values))))
                # Get the positions in lineData for each matching index
                matched_positions = [value_to_pos[int(allSuperPixelIDs[f'DMD{DMDix+1}'][idx])] for idx in matching_indices]
                
                data[matching_indices, DSframeIx] += allLineData[lineIdx][matched_positions,0]
                dataCt[matching_indices, DSframeIx] += 1
    
    return data/100, dataCt, baseline, aData

    return trial_peaks_gaussian, trial_original_peaks_gaussian, trial_peak_vals_gaussian

def process_timepoint(residual_frame, neighborMatrix_sparse):
    return np.max(neighborMatrix_sparse.multiply(np.tile(residual_frame,
                    (neighborMatrix_sparse.shape[1],1)).T).toarray(),axis=0)
                    
def gaussian_2d(xy, amplitude, x0, y0, sigma_x, sigma_y):
    x, y = xy
    exponent = -((x - x0)**2 / (2 * sigma_x**2) + (y - y0)**2 / (2 * sigma_y**2))
    return amplitude/(2*np.pi*sigma_x*sigma_y) * np.exp(exponent)

# make a function to get the nearest neighbors of a pixel in a sparse image
def getNeighbors(image, pixel):
    colNeighbors = getColNeighbors(image, pixel)
    rowNeighbors = getRowNeighbors(image, pixel)

    # combine colNeighbors and rowNeighbors
    neighbors = colNeighbors + rowNeighbors

    return neighbors

def getColNeighbors(image, pixel):
    # get the row and column of the pixel
    r = pixel[0]
    c = pixel[1]

    # column neighbors are the first pixels above or below that are not nan
    colNeighbors = np.where(~np.isnan(image[:,c]))[0]
    colNeighbors = colNeighbors[max(0,(colNeighbors < r).sum()-1):len(colNeighbors)-max(0,(colNeighbors > r).sum()-1)]

    colNeighbors = np.array([colNeighbors, np.ones(colNeighbors.shape)*c],dtype=np.int64)
    colNeighbors = colNeighbors[:,(colNeighbors[0] != r) | (colNeighbors[1] != c)]

    colNeighbors = list(zip(colNeighbors[0],colNeighbors[1]))

    return colNeighbors

def getRowNeighbors(image, pixel):
    r = pixel[0]
    c = pixel[1]

    # row neighbors are the first pixels left or right (or slightly diagonal) that are not nan
    rowNeighbors = np.where(~np.isnan(image[max(0,r-2):min(r+3,image.shape[0]),max(0,c-5):min(c+6,image.shape[1])]))

    # sort rowNeighbors by rowNeighbors[1]
    rowNeighbors = np.array(sorted(zip(rowNeighbors[0],rowNeighbors[1]),key=lambda x: x[1]),dtype=np.int64).T

    if (len(rowNeighbors) == 0):
        print(pixel)

    rowNeighbors[0] = rowNeighbors[0] + r-2
    rowNeighbors[1] = rowNeighbors[1] + c-5
    rowNeighbors = rowNeighbors[:,max(0,(rowNeighbors[1] < c).sum()-1):len(rowNeighbors[1])-max(0,(rowNeighbors[1] > c).sum()-1)]

    # remove current pixel from rowNeighbors
    rowNeighbors = rowNeighbors[:,(rowNeighbors[0] != r) | (rowNeighbors[1] != c)]

    rowNeighbors = list(zip(rowNeighbors[0],rowNeighbors[1]))

    return rowNeighbors

def find_trial_peaks(trialData,dmdPixelsPerColumn, dmdPixelsPerRow, refR, refC, params, DMDix, convMatrix, neighborMatrix_sparse):
    
    trial_ix, data, dataCt, baseline, (motionDSr, motionDSc, motionDSz) = trialData

    numCycles = motionDSr.shape[0]
    numSuperPixels = data.shape[0]

    # print(f"{trial_ix}; numSuperPixels: {numSuperPixels}, numCycles: {numCycles} data: {data.shape}, motionDSr: {motionDSr.shape}")

    # get peaks
    residual = data - baseline*dataCt
    residual_norm = residual / np.sqrt(baseline*dataCt)

    residual = residual / dataCt

    decayTau_frames = params['decayTau_s']*params['alignHz']

    decay_kernel = np.exp(np.linspace(-np.ceil(decayTau_frames*3),0,int(np.ceil(decayTau_frames*3)+1))/decayTau_frames)

    residual_filt = signal.convolve2d(convMatrix @ residual_norm,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same')

    # Find temporal peaks using vectorized operations
    temporalPeaks = (residual_filt > np.roll(residual_filt,1,axis=1)) & (residual_filt > np.roll(residual_filt,-1,axis=1))

    # Get max values for each neighbor group at each timepoint using matrix multiplication
    # This avoids the double loop and uses efficient sparse matrix operations
    neighbor_maxes = np.zeros_like(residual_filt)

    neighbor_sums = np.array(neighborMatrix_sparse.sum(axis=1)).flatten() 
    neighbor_mask = neighbor_sums > 0
    
    for t in range(residual_filt.shape[1]):
        neighbor_maxes[:,t] = process_timepoint(residual_filt[:,t], neighborMatrix_sparse)

    # Create mask for valid peaks using broadcasting
    valid_peaks = (neighbor_mask[:,None] & 
                    (residual_filt >= neighbor_maxes))

    # Apply masks to get final peaks
    allPeaks = temporalPeaks & valid_peaks

    # Get peak locations and values efficiently
    allPeakLocs = np.where(allPeaks.ravel())[0]
    peakVals = residual_filt.ravel()[allPeakLocs]

    kde = stats.gaussian_kde(peakVals)
    x_range = np.linspace(peakVals.min(), peakVals.max(), 1000)

    mode_idx = np.argmax(kde(x_range))
    bottom_half = peakVals[peakVals < x_range[mode_idx]]
    # bottom_half_kde = stats.gaussian_kde(bottom_half)
    bottom_half_x_range = np.linspace(bottom_half.min(), bottom_half.max(), 500)

    def truncated_normal(x, mu, sigma, K):
        return K * stats.truncnorm.pdf(x, -np.inf, 0, loc=mu, scale=sigma)

    # print('fitting peak value pdf...')
    popt, _ = curve_fit(lambda x, mu, sigma, K: truncated_normal(x, mu, sigma, K), bottom_half_x_range, kde(bottom_half_x_range), p0=[x_range[mode_idx], np.std(bottom_half), 1])

    fit_mu = popt[0]
    fit_sigma = popt[1]
    fit_K = popt[2]

    # plt.hist(peakVals,bins=100)
    # plt.plot(x_range,truncated_normal(x_range,fit_mu,fit_sigma,fit_K,-np.inf,np.inf),'r-')

    noise_dist = stats.norm.pdf(x_range,fit_mu,fit_sigma)*fit_K*2

    signal_dist = kde(x_range) - noise_dist
    signal_dist[x_range < fit_mu] = 0
    signal_dist[signal_dist < 0] = 0

    # plt.plot(x_range,kde(x_range),'k-')
    # plt.plot(x_range,noise_dist,'r-')
    # plt.plot(x_range,signal_dist,'b-')

    total_kde_integral = np.trapz(kde(x_range), x_range)
    # total_noise_dist = np.trapz(noise_dist, x_range)
    # total_signal_dist = np.trapz(signal_dist, x_range)

    # plt.plot(x_range,kde(x_range) / total_kde_integral,'k-')
    # plt.plot(x_range,noise_dist / total_kde_integral,'r-')
    # plt.plot(x_range,signal_dist / total_kde_integral,'b-')

    valid_x_range = (signal_dist > max(signal_dist) * 0.1) | (noise_dist > max(noise_dist) * 0.1)
    intersection_idx = np.argmin(np.abs(noise_dist[valid_x_range] - signal_dist[valid_x_range]) / (noise_dist[valid_x_range]+1e-8))
    intersection_idx = np.where(valid_x_range)[0][intersection_idx]
    peakVal_thresh = x_range[intersection_idx]
    intersection_y = noise_dist[intersection_idx]

    # # Plot intersection point
    # plt.plot(peakVal_thresh, intersection_y, 'go', label='Noise/Signal Intersection')

    # plt.xlabel('Peak Values')
    # plt.ylabel('Probability Density')
    # plt.title('PDF of Peak Values')
    # plt.legend()
    # plt.grid(True)
    # plt.show()

    # print(f"{trial_ix} plotting peak values pdf")
    plt.figure()
    plt.plot(x_range,kde(x_range) / total_kde_integral,'k-')
    plt.plot(x_range,noise_dist / total_kde_integral,'r-')
    plt.plot(x_range,signal_dist / total_kde_integral,'b-')
    plt.plot(peakVal_thresh, intersection_y / total_kde_integral, 'go')
    plt.xlabel('Peak Values')
    plt.ylabel('Probability Density')
    plt.title('PDF of Peak Values')
    plt.savefig(os.path.join(params['savedr'], f'peak_value_distributions_DMD{DMDix+1}_T{trial_ix+1}.png'))
    plt.close()

    # low = np.percentile(peakVals,1)
    # high = np.percentile(peakVals,50)

    # finalPeakLocs = allPeakLocs[peakVals > (high + 2*(high - low))] # pick top 6 sigma peaks
    # print(f"{trial_ix} choose top peaks")
    finalPeakLocs = allPeakLocs[peakVals > peakVal_thresh]

    finalPeaks = np.zeros((numSuperPixels,numCycles))
    finalPeaks[np.unravel_index(finalPeakLocs,(numSuperPixels,numCycles))] = 1

    # peaks_viewer = napari.view_image(residual_filt)
    # peaks_viewer.add_image(temporalPeaks)
    # peaks_viewer.add_image(allPeaks)
    # peaks_viewer.add_image(finalPeaks)
    # napari.run()

    trial_peaks = []
    trial_original_peaks = []
    trial_peak_vals = []

    for sp, t in zip(*np.where(finalPeaks)):
        tmpR = refR[sp]
        tmpC = refC[sp]

        newR = int(tmpR + np.round(motionDSr[t]))
        newC = int(tmpC + np.round(motionDSc[t]))

        residual_filt_frame = np.zeros((dmdPixelsPerColumn, dmdPixelsPerRow,len(decay_kernel)))
        residual_filt_frame[:] = np.nan

        for ix in range(len(decay_kernel)):
            dt = ix - len(decay_kernel)//2
            if t+dt >= 0 and t+dt < numCycles:
                residual_filt_frame[refR.int() + int(np.round(motionDSr[t+dt])), refC.int() + int(np.round(motionDSc[t+dt])),ix] = residual[:,t+dt]

        # residual_filt_frame[refR.int() + int(np.round(motionDSr[t])), refC.int() + int(np.round(motionDSc[t])),0] = residual_filt[:,t]
        # if t < numCycles-1:
        #     residual_filt_frame[refR.int() + int(np.round(motionDSr[t+1])), refC.int() + int(np.round(motionDSc[t+1])),1] = residual_filt[:,t+1]
        # if t > 0:
        #     residual_filt_frame[refR.int() + int(np.round(motionDSr[t-1])), refC.int() + int(np.round(motionDSc[t-1])),2] = residual_filt[:,t-1]

        with warnings.catch_warnings():
            warnings.filterwarnings("ignore", category=RuntimeWarning)
            residual_filt_frame = np.nanmean(residual_filt_frame,axis=2)

        patchWidth = 21
        peakPatch = residual_filt_frame[np.int64(newR)-patchWidth//2:np.int64(newR)+patchWidth//2+1,np.int64(newC)-patchWidth//2:np.int64(newC)+patchWidth//2+1]
        
        # Create x and y coordinates
        y, x = np.indices(peakPatch.shape)
        xy = np.column_stack((x.ravel(), y.ravel()))
        
        # Initial guess for parameters
        initial_guess = [np.nanmax(peakPatch), patchWidth//2, patchWidth//2, 2, 2]
        # Fit the 2D Gaussian
        try:
            # Ignore NaNs in peakPatch
            valid_mask = ~np.isnan(peakPatch.ravel())
            x_valid = x.ravel()[valid_mask]
            y_valid = y.ravel()[valid_mask]
            peakPatch_valid = peakPatch.ravel()[valid_mask]
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", category=UserWarning)
                popt, _ = curve_fit(gaussian_2d, (x_valid, y_valid), peakPatch_valid, p0=initial_guess)
            amplitude, x0, y0, sigma_x, sigma_y = popt
            
            # Update peak position based on the Gaussian fit
            peakR = np.int64(newR) - patchWidth//2 + y0
            peakC = np.int64(newC) - patchWidth//2 + x0

            if amplitude <= 0 or sigma_x >= 5 or sigma_y >= 15 or sigma_x < 1 or sigma_y < 1 or y0 < 0 or y0 >= patchWidth or x0 < 0 or x0 >= patchWidth:
                # trial_peaks.append((peakR,peakC,t))
                # trial_original_peaks.append((np.int64(newR),np.int64(newC),t))
                # trial_peak_vals.append(residual_filt[sp,t])
                continue

            trial_peaks.append((peakR,peakC,t))
            trial_original_peaks.append((np.int64(newR),np.int64(newC),t))
            trial_peak_vals.append(amplitude)
            # trial_peak_vals_gaussian.append(residual_filt[sp,t])

            # trial_peak_sds.append((sigma_x,sigma_y))
        except:
            continue
            # If fitting fails, use quadratic interpolation
            # rowNeighbors = getRowNeighbors(residual_filt_frame,(np.int64(newR),np.int64(newC)))
            # colNeighbors = getColNeighbors(residual_filt_frame,(np.int64(newR),np.int64(newC)))

            # if len(rowNeighbors) >= 2:
            #     # ratioC = max(1e-6,(R_filt_frame[np.int64(newR),np.int64(newC)]-R_filt_frame[rowNeighbors[0][0],rowNeighbors[0][1]])/(R_filt_frame[np.int64(newR),np.int64(newC)]-R_filt_frame[rowNeighbors[1][0],rowNeighbors[1][1]]))
            #     # dC = (1-ratioC)/(1+ratioC)/2
            #     # peakC = np.int64(newC) - dC

            #     peakC = (np.int64(newC)*residual_filt_frame[np.int64(newR),np.int64(newC)] + rowNeighbors[0][1]*residual_filt_frame[rowNeighbors[0][0],rowNeighbors[0][1]] + rowNeighbors[1][1]*residual_filt_frame[rowNeighbors[1][0],rowNeighbors[1][1]])/(residual_filt_frame[np.int64(newR),np.int64(newC)] + residual_filt_frame[rowNeighbors[0][0],rowNeighbors[0][1]] + residual_filt_frame[rowNeighbors[1][0],rowNeighbors[1][1]])
            # else:
            #     continue
            #     peakC = np.int64(newC)

            # if len(colNeighbors) >= 2:
            #     # ratioR = max(1e-6,(R_filt_frame[np.int64(newR),np.int64(newC)]-R_filt_frame[colNeighbors[0][0],colNeighbors[0][1]])/(R_filt_frame[np.int64(newR),np.int64(newC)]-R_filt_frame[colNeighbors[1][0],colNeighbors[1][1]]))
            #     # dR = (1-ratioR)/(1+ratioR)/2 * 4.5
            #     # peakR = np.int64(newR) - dR

            #     peakR = (np.int64(newR)*residual_filt_frame[np.int64(newR),np.int64(newC)] + colNeighbors[0][0]*residual_filt_frame[colNeighbors[0][0],colNeighbors[0][1]] + colNeighbors[1][0]*residual_filt_frame[colNeighbors[1][0],colNeighbors[1][1]])/(residual_filt_frame[np.int64(newR),np.int64(newC)] + residual_filt_frame[colNeighbors[0][0],colNeighbors[0][1]] + residual_filt_frame[colNeighbors[1][0],colNeighbors[1][1]])
            # else:
            #     continue
            #     peakR = np.int64(newR)
    return trial_peaks, trial_original_peaks, trial_peak_vals

def get_traces(trial_ix, roi_masks, DMDix, trialTable):
    # print(f"Getting trial {trial_ix+1} traces")
    integration_movie = skimio.imread(trialTable['fnRegDSInt'][DMDix,trial_ix][0])
    integration_movie = integration_movie[::2]
    integration_movie[np.isnan(integration_movie)] = 0
    
    traces, _ , _ , _ = np.linalg.lstsq(roi_masks.reshape((roi_masks.shape[0],-1)).T, integration_movie.T.reshape((-1,integration_movie.shape[0])), rcond=None)

    return traces

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
            tiff_fn = os.path.splitext(os.path.basename(trialTable['fnRegDS'][DMDix,trialIx][0]))[0]
            if not os.path.exists(os.path.join(dr, tiff_fn + '.tif')):
                print(f'Missing tiff file: {tiff_fn}')
                keepTrials[DMDix,trialIx] = False
                
            # Check alignment data files
            align_fn = os.path.splitext(os.path.basename(trialTable['fnAdata'][DMDix,trialIx][0]))[0]
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
    params['analyzeHz'] = 200
    params['activityChannel'] = 1
    params['savedr'] = savedr
    params['decayTau_s'] = 0.03

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

    del sparseMaskInds_refs, allSuperPixelIDs_refs

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
            numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
            IM_dil5 = trialTable['refStack'][0,DMDix]['IM'][0,0].transpose(2,1,0)
            IM_dil5 = IM_dil5.reshape(-1, numChannels, IM_dil5.shape[1], IM_dil5.shape[2]).transpose(1,0,2,3)

            tmp = deconvlucy.deconvlucy(refStack[f'DMD{DMDix+1}'][0], IM_dil5[0], num_iter=50)
            tmp = tmp[np.argmax(np.sum(tmp, axis=(1,2)))]

            # Fit elongated 2D Gaussian to the PSF
            def gaussian2d(xy, amplitude, x0, y0, sigma_x, sigma_y):
                x, y = xy
                x_centered = x - x0
                y_centered = y - y0
                
                # Calculate Gaussian without rotation
                exp_term = (x_centered**2)/(2*sigma_x**2) + (y_centered**2)/(2*sigma_y**2)
                return amplitude * np.exp(-exp_term)

            # Create coordinate grid
            y, x = np.indices(tmp.shape)
            xy = (x, y)

            # Initial parameter guesses
            amplitude_guess = np.max(tmp)
            y0_guess, x0_guess = np.unravel_index(np.argmax(tmp), tmp.shape)
            sigma_guess = 3

            initial_guess = [amplitude_guess, x0_guess, y0_guess, sigma_guess, sigma_guess]

            # Flatten data for fitting
            tmp_flat = tmp.ravel()

            # Fit the 2D Gaussian
            popt, _ = curve_fit(lambda xy, *params: gaussian2d(xy, *params).ravel(),
                                xy, tmp_flat, p0=initial_guess)

            # Store fitted parameters
            amplitude, x0, y0, sigma_x, sigma_y = popt
            tmp = tmp[round(y0)-round(2.5*sigma_y):round(y0)+round(2.5*sigma_y)+1, round(x0)-round(2.5*sigma_x):round(x0)+round(2.5*sigma_x)+1]

            psf[f'DMD{DMDix+1}'] = (tmp+tmp[-1::-1,:]+tmp[:,::-1]+tmp[-1::-1,::-1])/4

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
    
    psf = {
        'DMD1': psf_combined[0],
        'DMD2': psf_combined[1]
    }

    del psf_combined

    for DMDix in range(nDMDs-1, -1, -1):
        print(f'Processing DMD{DMDix+1}')

        dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
        dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
        numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
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
        refD = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
        refC = torch.div((refPixs - refD * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
        refR = refPixs % dmdPixelsPerColumn

        openPixs = torch.from_numpy(sparseMaskInds[f'DMD{DMDix+1}'][:,0]-1)
        d = torch.div(openPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
        c = torch.div((openPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
        r = openPixs % dmdPixelsPerColumn

        filterSize = psf[f'DMD{DMDix+1}'].shape[0]*psf[f'DMD{DMDix+1}'].shape[1]
        sparseHInds = np.zeros((2,numSuperPixels*filterSize))
        sparseHVals = np.zeros((numSuperPixels*filterSize,))

        for spIdx in range(subsampleMatrixInds.shape[0]):
            tmpMap = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow))
            tmpMap[:] = np.nan

            tmpMap[refR[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[0]//2+1,refC[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

            sparseHInds[0,spIdx*filterSize:(spIdx+1)*filterSize] = subsampleMatrixInds[spIdx,1] - 1
            sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize] = np.where(~np.isnan(tmpMap.flatten()))[0]

            sparseHVals[spIdx*filterSize:(spIdx+1)*filterSize] = tmpMap.flatten()[sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize].astype(np.uint32)]

        convMatrix = np.zeros((numSuperPixels,numSuperPixels))

        for spIdx in range(subsampleMatrixInds.shape[0]):
            tmpMap = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow))
            tmpMap[refR[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[0]//2+1,refC[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

            convMatrix[spIdx,:] = tmpMap[refR.int(),refC.int()]

        trial_info = [(i, keepTrials[DMDix,i]) for i in range(nTrials)]
        # processes = []

        get_trial_data_partial = partial(
            get_trial_data,
            DMDix=DMDix,
            params=params,
            refStack=refStack,
            subsampleMatrixInds=subsampleMatrixInds,
            sparseHInds=sparseHInds,
            sparseHVals=sparseHVals,
            allSuperPixelIDs=allSuperPixelIDs,
            dr=dr,
            trialTable=trialTable
        )

        # for info in trial_info:
        #     p = mp.Process(target=get_trial_data_partial, args=(info,))
        #     p.start()
        #     processes.append(p)

        # for p in processes:
        #     p.join()

        with mp.Pool(processes=min(mp.cpu_count(),len(trial_info))) as pool:
            results = list(pool.imap(get_trial_data_partial, trial_info))

        # Concatenate trial data
        # data = [r[0] for r in results]
        # dataCt = [r[1] for r in results]
        # baseline = [r[2] for r in results]
        # motionDSr = [r[3]['motionDSr'] for r in results]
        # motionDSc = [r[3]['motionDSc'] for r in results]
        # motionDSz = [r[3]['motionDSz'] for r in results]

        neighborMatrix = np.zeros((numSuperPixels,numSuperPixels))

        refPixIm = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow))
        refPixIm[:] = np.nan
        refPixIm[refR.int(),refC.int()] = 1

        for sp in range(numSuperPixels):
            neighbors = getNeighbors(refPixIm,(np.int64(refR[sp]),np.int64(refC[sp])))
            # find the indices of the neighbors in subsampleMatrixInds
            if len(neighbors) > 1:
                neighborInds = []
                for neigh in neighbors:
                    tmp = np.where((refR.int() == neigh[0]) & (refC.int() == neigh[1]))[0]
                    neighborInds.append(tmp)

                neighborMatrix[sp,neighborInds] = 1

        # Convert sparse neighborMatrix to CSR format for efficient operations
        neighborMatrix_sparse = sparse.csr_matrix(neighborMatrix)

        trialInfo = [(ix, r[0], r[1], r[2], (r[3]['motionDSr'], r[3]['motionDSc'], r[3]['motionDSz'])) for ix, r in enumerate(results)]

        find_trial_peaks_partial = partial(
            find_trial_peaks,
            dmdPixelsPerColumn = dmdPixelsPerColumn,
            dmdPixelsPerRow = dmdPixelsPerRow,
            refR = refR,
            refC = refC,
            params = params,
            DMDix = DMDix,
            convMatrix = convMatrix,
            neighborMatrix_sparse = neighborMatrix_sparse
        )

        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = list(tqdm(pool.imap_unordered(find_trial_peaks_partial, trialInfo), total=len(trialInfo),desc="Fitting Peaks"))      

        peaks = []
        original_peaks = []
        peak_vals = []
        peak_trials = []

        for trial_ix, (peak, original_peak, peak_val) in enumerate(results):
            peaks.extend(peak)
            original_peaks.extend(original_peak) 
            peak_vals.extend(peak_val)
            peak_trials.extend([trial_ix]*len(peak))

        # Save all peak data to a single file
        peak_data = {
            'peaksPix': np.array(peaks),
            'originalPeaks': np.array(original_peaks),
            'finalPeakVals': np.array(peak_vals),
            'peakTrials': np.array(peak_trials)
        }
        output_filename = os.path.join(params['savedr'], f'peak_data_DMD{DMDix+1}.npz')
        np.savez(output_filename, **peak_data)
        print(f"Saved peak data to {output_filename}")

        peaksPix = np.array(peaks)
        originalPeaks = np.array(original_peaks)
        finalPeakVals = np.array(peak_vals)
        peakTrials = np.array(peak_trials)

        # peaks histogram
        peaksHist = np.zeros((dmdPixelsPerColumn//2,dmdPixelsPerRow*1))

        scaled_finalPeakVals = np.clip(finalPeakVals, None, np.percentile(finalPeakVals, 95))
        # scaled_finalPeakVals = scaled_finalPeakVals - np.min(scaled_finalPeakVals)
        # scaled_finalPeakVals = scaled_finalPeakVals / np.max(scaled_finalPeakVals)

        for i in range(len(peaksPix)):
            peaksHist[round(peaksPix[i][0]/2),round(peaksPix[i][1]*1)] += np.sqrt(scaled_finalPeakVals[i])

        peaksHist[peaksHist == 0] = np.nan

        roi_viewer = napari.view_image(refStack[f'DMD{DMDix+1}'][0], colormap='cyan')

        # Create a napari viewer with the correct aspect ratio
        roi_viewer.add_image(
            peaksHist,
            name='Peaks Histogram',
            colormap='red',
            scale=(2, 1)
        )

        roi_viewer.add_points(np.array(peaksPix)[:,:2], size=1, face_color='magenta', border_color='magenta',opacity=0.3)

        roi_layer = roi_viewer.add_shapes(name='ROIs',face_color=[1,1,0,0.1],edge_color='yellow',edge_width=1)
        roi_layer.mode = 'add_polygon'

        # Use a list to store both ROIs and done status
        result = {'roi_verts': None, 'roi_shapes': None, 'done': False}

        # Check if saved ROIs exist
        masks_filename = os.path.join(params['savedr'], f'roi_masks_DMD{DMDix+1}.npy')
        if os.path.exists(masks_filename):
            load_saved = QMessageBox.question(None, 'Question', 'Load saved ROIs?') == QMessageBox.Yes
            if load_saved:
                print('Loading saved ROIs')
                # Load saved masks
                result['roi_verts'], result['roi_shapes'] = np.load(masks_filename,allow_pickle=True)
                
                # # Convert masks back to vertices
                # roi_verts = []
                # roi_shapes = []
                # for mask in masks_array:
                #     # Find contours of the mask
                #     contours = measure.find_contours(mask, 0.5)
                #     if len(contours) > 0:
                #         # Take the largest contour
                #         contour = contours[0]
                #         roi_verts.append(contour)
                #         roi_shapes.append('polygon')
                
                # Add shapes to viewer
                # roi_layer.add_shapes(roi_verts, shape_type=roi_shapes)

                roi_layer.data = result['roi_verts']
                roi_layer.shape_type = result['roi_shapes']

                result['done'] = False
                napari.utils.notifications.show_info('Loaded saved ROIs')

        if not result['done']:
            # Callback for Done button
            def on_done():
                if 'ROIs' in roi_viewer.layers:
                    result['roi_verts'] = roi_layer.data
                    result['roi_shapes'] = roi_layer.shape_type
                    result['done'] = True
                    napari.utils.notifications.show_info('ROIs captured!')

            # Add Done button
            done_button = QPushButton("Done")
            done_button.clicked.connect(on_done)

            roi_viewer.window.add_dock_widget(done_button, area='right')
            napari.utils.notifications.show_info('Draw ROIs and click "Done" when ready.')

            # Keep checking until Done is clicked
            while not result['done']:
                QApplication.processEvents()
                time.sleep(0.1)

        if len(result['roi_verts']) > 0:
            np.save(masks_filename, np.array([result['roi_verts'],result['roi_shapes']], dtype=object))

            # Create masks for each ROI
            masks = []
            for i, roi_coords in enumerate(result['roi_verts']):
                mask = np.zeros((dmdPixelsPerColumn, dmdPixelsPerRow), dtype=bool)
                coords = np.round(roi_coords).astype(int)

                if result['roi_shapes'][i] in ['polygon', 'rectangle']:
                    rr, cc = draw.polygon(coords[:,0], coords[:,1])
                elif result['roi_shapes'][i] == 'ellipse':
                    center_y = (coords[0][0] + coords[2][0]) / 2
                    center_x = (coords[0][1] + coords[2][1]) / 2
                    height = (coords[2][0] - coords[0][0]) / 2  # semi-major axis
                    width = (coords[1][1] - coords[0][1]) / 2   # semi-minor axis
                    rr, cc = draw.ellipse(center_y, center_x, height, width)
                else:
                    raise ValueError(f"Unknown ROI shape: {result['roi_shapes'][i]}")
                
                valid = (rr >= 0) & (rr < dmdPixelsPerColumn) & (cc >= 0) & (cc < dmdPixelsPerRow)
                mask[rr[valid], cc[valid]] = True
                masks.append(mask)

            # Save masks to file
            masks_array = np.array(masks)
            print(f"Saved {len(masks)} ROI masks to {masks_filename}")

            get_traces_partial = partial(
                get_traces,
                roi_masks = masks_array,
                DMDix = DMDix,
                trialTable = trialTable
            )

            with mp.Pool(processes=4) as pool:
                results = list(tqdm(pool.imap(get_traces_partial, range(nTrials)), total=nTrials, desc="Getting Traces"))
            
            np.savez(os.path.join(params['savedr'], f'traces_DMD{DMDix+1}.npz'), *results)

            # def plot_time_series(traces, time_points=None, title='Time Series'):
            #     """
            #     Create interactive plotly figure for multiple time series
                
            #     Parameters:
            #     traces: np.array of shape (n_traces, n_timepoints)
            #     time_points: optional array of time values
            #     title: string for plot title
                
            #     Returns:
            #     None (saves interactive HTML file)
            #     """
                
            #     # Create time points if not provided
            #     if time_points is None:
            #         time_points = np.arange(traces.shape[1])
                
            #     # Create figure
            #     fig = go.Figure()
                
            #     # Add each time series
            #     for i in range(traces.shape[0]):
            #         fig.add_trace(
            #             go.Scatter(
            #                 x=time_points,
            #                 y=traces[i],
            #                 name=f'Trace {i+1}',
            #                 mode='lines',
            #                 hovertemplate='Time: %{x}<br>Value: %{y:.2f}<extra></extra>'
            #             )
            #         )
                
            #     # Update layout
            #     fig.update_layout(
            #         title=title,
            #         xaxis_title='Time Point',
            #         yaxis_title='Value',
            #         showlegend=True,
            #         # Add range slider
            #         xaxis=dict(
            #             rangeslider=dict(visible=True),
            #             type='linear'
            #         ),
            #         # Make it more interactive
            #         hovermode='x unified'
            #     )
                
            #     # Save as HTML file
            #     fig.write_html(os.path.join(params['savedr'],f'traces_plot_DMD{DMDix+1}.html'))
                
            #     # Optionally display in notebook
            #     # fig.show()

            # # Use the function
            # plot_time_series(traces, title='20 Time Series')

        # break

if __name__ == '__main__':
    main()

# integrationMovie = np.zeros((residual_filt.shape[1],800,1280))
# integrationMovie[:] = np.nan

# residual_filt_movie = np.zeros((residual_filt.shape[1],800,1280))
# residual_filt_movie[:] = np.nan

# for t in range(residual_filt.shape[1]):
#     integrationMovie[t,refR.int() + int(np.round(motionDSr[t])), refC.int() + int(np.round(motionDSc[t]))] = data[:,t] / 100 / dataCt[:,t]
#     residual_filt_movie[t,refR.int() + int(np.round(motionDSr[t])), refC.int() + int(np.round(motionDSc[t]))] = residual_filt[:,t]

# peaks3D = []
# for sp,t in zip(*np.where(finalPeaks)):
#     peaks3D.append((t, refR[sp].int() + int(np.round(motionDSr[t])), refC[sp].int() + int(np.round(motionDSc[t]))))
