import sys
from pathlib import Path
import numpy as np
import torch
import scipy.io as spio
import scipy.ndimage as ndimage
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

def get_trial_data(trial_info, DMDix, samp_freq, params, refStack, subsampleMatrixInds, fastZ2RefZ, sparseHInds, sparseHVals, 
                allSuperPixelIDs, dr, trialTable):
    trialIx, keepTrial = trial_info
    
    if not keepTrial:
        return [], [], []

    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[2]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[3]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[1]
    numSuperPixels = allSuperPixelIDs[f'DMD{DMDix+1}'].shape[0]
    numFastZs = fastZ2RefZ[f'DMD{DMDix+1}'].shape[0]

    nPixels = dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs

    aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]

    numCycles = aData['motionDSr'].shape[0]
    uniqueMotion, motInds = np.unique(np.round(np.concatenate((aData['motionDSr'],aData['motionDSc'],aData['motionDSz']),axis=1)),axis=0,return_inverse=True)

    Afinal = torch.empty((nPixels,0))
    phiFinal = torch.empty((numCycles,0))

    importlib.reload(reconstruct)

    startTime = time.time()
    # baseline = getBaseline(Afinal,phiFinal,refStack[f'DMD{DMDix+1}'][0],torch.from_numpy(aData['brightnessDS']),subsampleMatrixInds.T,fastZ2RefZ[f'DMD{DMDix+1}'],sparseHInds,sparseHVals,uniqueMotion,motInds)
    baseline = reconstruct.reconstruct(Afinal,phiFinal,refStack[f'DMD{DMDix+1}'][0],torch.from_numpy(aData['brightnessDS']),subsampleMatrixInds.T,fastZ2RefZ[f'DMD{DMDix+1}'],sparseHInds,sparseHVals,uniqueMotion,motInds).detach().numpy()
    endTime = time.time()

    # print(f"Elapsed Time: {endTime - startTime} sec")

    # plt.imshow(baseline)
    # plt.colorbar()
    # plt.show()

    source_fn = trialTable['filename'][DMDix,trialIx][0]
    firstLine = trialTable['firstLine'][DMDix,trialIx]
    lastLine = trialTable['lastLine'][DMDix,trialIx]

    importlib.reload(slap2_utils)
    hDataFile = slap2_utils.DataFile(os.path.join(dr, source_fn))

    dt = 1/samp_freq/hDataFile.metaData.linePeriod_s

    DSframes = aData['DSframes'][0]
    nDSframes= len(DSframes)
    
    # Pre-compute time windows for all frames
    timeWindows = [np.arange(max(1,np.floor(DSframes[i]-2*dt)), 
                            min(np.ceil(DSframes[i]+2*dt),hDataFile.numCycles*hDataFile.header['linesPerCycle'])+1) 
                    for i in range(nDSframes)]
    
    # Pre-compute line and cycle indices for all time windows
    lineIndices_all = [(tw - 1) % hDataFile.header['linesPerCycle'] + 1 for tw in timeWindows]
    cycleIndices_all = [np.floor((tw - 1) / hDataFile.header['linesPerCycle']) + 1 for tw in timeWindows]
    
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
                
                data[matching_indices, DSframeIx] += line_data[matched_positions, 0]
                dataCt[matching_indices, DSframeIx] += 1
    
    return data/100, dataCt, baseline, aData

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
            samp_freq=params['alignHz'],
            params=params,
            refStack=refStack,
            subsampleMatrixInds=subsampleMatrixInds,
            fastZ2RefZ=fastZ2RefZ,
            sparseHInds=sparseHInds,
            sparseHVals=sparseHVals,
            allSuperPixelIDs=allSuperPixelIDs,
            dr=dr,
            trialTable=trialTable
        )

        # if params.get('parallel', True):
        with mp.Pool(processes=min(mp.cpu_count(),len(trial_info))) as pool:
            results = list(pool.imap(get_trial_data_partial, trial_info))
        # else:
        #     results = [get_trial_data_partial(t) for t in trial_info]

        lowResData = np.concatenate([r[0] for r in results], axis=1)
        lowResDataCt = np.concatenate([r[1] for r in results], axis=1)
        lowResBaseline = np.concatenate([r[2] for r in results], axis=1)
        lowResMotionR = np.concatenate([r[3]['motionDSr'] for r in results], axis=0)
        lowResMotionC = np.concatenate([r[3]['motionDSc'] for r in results], axis=0)
        lowResMotionZ = np.concatenate([r[3]['motionDSz'] for r in results], axis=0)
        lowResBrightness = np.concatenate([r[3]['brightnessDS'] for r in results], axis=0)

        uniqueMotion, motInds = np.unique(np.round(np.concatenate((lowResMotionR,lowResMotionC,lowResMotionZ),axis=1)),axis=0,return_inverse=True)
        motIndsToKeep = (np.bincount(motInds) > 100).nonzero()[0]

        uniqueMotionYX = np.unique(uniqueMotion[:,:2],axis=0)

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


        # read in source seeds
        source_seeds = np.load(os.path.join(params['savedr'], f'roi_seeds_DMD{DMDix+1}.npy'))
        nSources = source_seeds.shape[0]

        source_params = torch.cat([
            torch.ones(nSources, 1, dtype=torch.float32),
            torch.tensor(source_seeds, dtype=torch.float32),
            torch.ones(nSources, 2, dtype=torch.float32) * 2
        ], dim=1)

        A = torch.zeros((nPixels, nSources))

        y, x = pixel_coords_tensor[:, 1], pixel_coords_tensor[:, 2]
        y_means = source_params[:, 1].unsqueeze(0)  # Shape: [1, nSources]
        x_means = source_params[:, 2].unsqueeze(0)  # Shape: [1, nSources]
        y_sigmas = source_params[:, 3].unsqueeze(0)  # Shape: [1, nSources]
        x_sigmas = source_params[:, 4].unsqueeze(0)  # Shape: [1, nSources]
        amplitudes = source_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]

        y_term = -0.5 * ((y.unsqueeze(1) - y_means) / y_sigmas) ** 2  # Shape: [nPixels, nSources]
        x_term = -0.5 * ((x.unsqueeze(1) - x_means) / x_sigmas) ** 2  # Shape: [nPixels, nSources]
        A[selPixIdxs,:] = (amplitudes * torch.exp(y_term + x_term))

        motion_mode_idx = np.bincount(motInds).argmax()

        # Extract data for the most common motion mode
        motion_mode_frames = (motInds == motion_mode_idx).nonzero()[0]
        
        decayTau_frames = params['decayTau_s']*params['alignHz']
        decay_kernel = np.exp(np.linspace(-np.ceil(decayTau_frames*3),0,int(np.ceil(decayTau_frames*3)+1))/decayTau_frames)
        data_for_nmf = lowResData / lowResDataCt # - lowResDataCt[:, motion_mode_frames] * lowResBaseline[:, motion_mode_frames]
        data_tensor = torch.from_numpy(data_for_nmf[:, motion_mode_frames].astype(np.float32))
        # data_for_nmf = signal.convolve2d(data_for_nmf,np.expand_dims(decay_kernel / np.sum(decay_kernel),0),mode='same')
        # data_tensor = torch.from_numpy(data_for_nmf[:, motion_mode_frames].astype(np.float32))

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_mode_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_mode_idx,1].astype(int)

        # H = torch.sparse_coo_tensor(sparseHIndsShifted,sparseHVals,(numSuperPixels,nPixels),dtype=torch.float32)

        sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

        # project image space (A) into superpixel space (X)
        X = torch.sparse.mm(H, A[selPixIdxs,:])
        X = torch.clamp(X, min=1e-4) - 1e-4

        newR = refR+uniqueMotion[motion_mode_idx,0].astype(int)
        newC = refC+uniqueMotion[motion_mode_idx,1].astype(int)
        newD = torch.from_numpy(fastZ2RefZ[f'DMD{DMDix+1}'][refD.numpy()].flatten()-1)+uniqueMotion[motion_mode_idx,2].astype(int)

        # add background spatial component
        # background_spatial_component = torch.as_tensor(refStack[f'DMD{DMDix+1}'][0][newD,newR,newC].reshape((-1,1)))
        background_spatial_component = torch.median(data_tensor,dim=1)[0]
        background_spatial_component = background_spatial_component / torch.norm(background_spatial_component)
        X = torch.concat((X,background_spatial_component.unsqueeze(-1)),dim=1)

        # Initialize phi_lowRes for NMF
        # phi_lowRes = torch.rand(len(motion_mode_frames), nSources, dtype=torch.float32) * 0.1
        
        # initialize spatial components
        X_tensor = X.clone().detach() / torch.norm(X,dim=0,keepdim=True)

        # NMF parameters
        als_nmf_iters = 4
        mult_nmf_max_iters = 200
        outer_loop_iters = 3
        nmf_tol = 1e-6
        # Gaussian fitting optimization parameters
        learning_rate = 0.01
        num_epochs = 100
        gd_tol = 1e-8

        for outer_loop_iter in range(outer_loop_iters):
            
            # Initialize variables for NMF
            error_values = []
            prev_reconstruction_error = float('inf')

            X_tensor = X_tensor.clone().detach() / torch.norm(X_tensor,dim=0,keepdim=True)
            X_support = X_tensor > 0
            
            for iter_idx in tqdm(range(als_nmf_iters),desc='Alternating least squares NMF'):
                # Update phi (temporal components) with X fixed
                # Formulate as least squares: min ||data_t - X * phi_t||^2 for all t
                XtX = X_tensor.T @ X_tensor
                Xtd = X_tensor.T @ data_tensor  # This gives all time points at once
                
                # Add small regularization to ensure stability
                regularized_XtX = XtX + 1e-10 * torch.eye(XtX.shape[0])
                
                # Solve the system for all time points at once
                # We need to solve (X^T * X) * phi = X^T * data for each column of data
                phi_lowRes = torch.linalg.solve(
                    regularized_XtX,
                    Xtd
                ).T
                
                # Ensure non-negativity
                phi_lowRes = torch.clamp(phi_lowRes, min=0)
                
                # Update X (spatial components) with phi fixed
                # Formulate as least squares: min ||data - X * phi^T||^2 for all sources
                phi_phi_t = phi_lowRes.T @ phi_lowRes  # Shape: [sources, sources]
                data_phi = data_tensor @ phi_lowRes  # Shape: [pixels, sources]

                regularized_phi_phi_t = phi_phi_t + 1e-10 * torch.eye(phi_phi_t.shape[0])
                
                X_update = torch.linalg.solve(regularized_phi_phi_t, data_phi.T).T
                
                # Ensure non-negativity
                X_tensor = torch.clamp(X_update, min=0)
                X_tensor[~X_support] = 0
                X_tensor[:,-1] = background_spatial_component
                
                # Normalize X and phi to avoid scaling ambiguity
                for s in range(nSources):
                    norm = torch.norm(X_tensor[:, s])
                    if norm > 0:
                        X_tensor[:, s] = X_tensor[:, s] / norm
                        phi_lowRes[:, s] = phi_lowRes[:, s] * norm
                
                # Calculate reconstruction error
                reconstruction = X_tensor @ phi_lowRes.T
                current_error = torch.mean((data_tensor - reconstruction)**2).item()
                error_values.append(current_error)
                
                # Check convergence
                if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    break
                    
                prev_reconstruction_error = current_error
            for iter_idx in tqdm(range(mult_nmf_max_iters),desc='Multiplicative NMF'):
                # Multiplicative update for NMF
                # Update phi (temporal components) using multiplicative update rule
                # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
                numerator = X_tensor.T @ data_tensor
                denominator = (X_tensor.T @ X_tensor) @ phi_lowRes.T + 1e-10
                phi_update = phi_lowRes * (numerator / denominator).T
                phi_lowRes = phi_update
                
                # Update X (spatial components) using multiplicative update rule
                # X = X * (data * phi^T) / (X * phi * phi^T + epsilon)
                numerator = data_tensor @ phi_lowRes
                denominator = X_tensor @ (phi_lowRes.T @ phi_lowRes) + 1e-10
                X_update = X_tensor * (numerator / denominator)

                
                # Apply spatial support constraint
                X_tensor = X_update
                X_tensor[~X_support] = 0
                X_tensor[:,-1] = background_spatial_component
                
                # Normalize X and phi to avoid scaling ambiguity
                for s in range(nSources):
                    norm = torch.norm(X_tensor[:, s])
                    if norm > 0:
                        X_tensor[:, s] = X_tensor[:, s] / norm
                        phi_lowRes[:, s] = phi_lowRes[:, s] * norm
                
                # Calculate reconstruction error
                reconstruction = X_tensor @ phi_lowRes.T
                current_error = torch.mean((data_tensor - reconstruction)**2).item()
                error_values.append(current_error)
                
                # Check convergence
                if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    break
                    
                prev_reconstruction_error = current_error

            sortorder = np.argsort(-np.sum((phi_lowRes[:,:nSources].numpy()-np.mean(phi_lowRes[:,:nSources].numpy(),axis=0))**2,axis=0))
            source_params = source_params[sortorder,:]
            sortorder = np.append(sortorder,max(sortorder)+1)
            phi_lowRes = phi_lowRes[:,sortorder]
            X_tensor = X_tensor[:,sortorder]

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
            
            # Initialize parameters based on the current X_tensor
            for s in range(nSources):
                # Get the current spatial component
                source_weights = X_tensor[:, s].numpy()
                source_weights_filtered = convMatrix @ source_weights
                
                # Find the maximum value and its position
                max_idx = np.argmax(source_weights_filtered)
                max_val = source_weights[max_idx]

                frame = np.zeros((800,1280))
                frame[:] = np.nan
                frame[refR + uniqueMotion[motion_mode_idx,0].astype(int), refC + uniqueMotion[motion_mode_idx,1].astype(int)] = source_weights
                
                patch_size = 21  # 5x5 patch around maximum
                max_r = refR[max_idx] + uniqueMotion[motion_mode_idx,0].astype(int)
                max_c = refC[max_idx] + uniqueMotion[motion_mode_idx,1].astype(int)
                
                # Define patch boundaries
                r_start = max(0, max_r - patch_size//2)
                r_end = min(800, max_r + patch_size//2 + 1)
                c_start = max(0, max_c - patch_size//2)
                c_end = min(1280, max_c + patch_size//2 + 1)
                
                # Extract patch from frame
                patch = frame[r_start:r_end, c_start:c_end]
                del frame
                
                # Calculate centroid if patch has valid values
                if not np.all(np.isnan(patch)):
                    # Create coordinate grids for the patch
                    r_coords, c_coords = np.meshgrid(
                        np.arange(r_start, r_end),
                        np.arange(c_start, c_end),
                        indexing='ij'
                    )
                    
                    # Get valid (non-NaN) values
                    valid_mask = ~np.isnan(patch)
                    if np.any(valid_mask):
                        weights = patch[valid_mask]
                        weights = np.maximum(weights, 0)  # Ensure non-negative weights
                        
                        if np.sum(weights) > 0:
                            # Calculate weighted centroid
                            centroid_r = np.sum(r_coords[valid_mask] * weights) / np.sum(weights)
                            centroid_c = np.sum(c_coords[valid_mask] * weights) / np.sum(weights)
                        else:
                            print("No valid values in patch")
                            centroid_r = max_r
                            centroid_c = max_c
                    else:
                        print("No valid values in patch")
                        centroid_r = max_r
                        centroid_c = max_c


                # Use the position of maximum value as initial mean
                # source_params[s, 0] = torch.tensor(max_val) # amplitude
                source_params[s, 1] = centroid_r / 10
                source_params[s, 2] = centroid_c / 10
                source_params[s,3:] = source_params[s,3:] / 10
            
            # Optimization loop for all sources together
            print("Fitting all Gaussian sources simultaneously")
            
            # Make parameters require gradients
            optim_params = source_params[:,1:].clone().requires_grad_(True)
            
            # Initialize Adam optimizer
            optimizer = torch.optim.Adam([optim_params], lr=learning_rate)
            
            # Optimization loop
            losses = []
            for epoch in tqdm(range(num_epochs),desc='Gaussian fitting'):
                # Zero gradients
                optimizer.zero_grad()
                
                # Forward pass: compute Gaussian values for all sources
                # Compute all Gaussian values at once using vectorized operations
                y, x = pixel_coords_tensor[:, 1], pixel_coords_tensor[:, 2]
                y_means = optim_params[:, 0].unsqueeze(0) * 10  # Shape: [1, nSources]
                x_means = optim_params[:, 1].unsqueeze(0) * 10  # Shape: [1, nSources] 
                y_sigmas = optim_params[:, 2].unsqueeze(0) * 10  # Shape: [1, nSources]
                x_sigmas = optim_params[:, 3].unsqueeze(0) * 10 # Shape: [1, nSources]
                # amplitudes = optim_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]
                
                # Compute terms for all sources at once
                y_term = -0.5 * ((y.unsqueeze(1) - y_means) / y_sigmas) ** 2  # Shape: [nPixels, nSources]
                x_term = -0.5 * ((x.unsqueeze(1) - x_means) / x_sigmas) ** 2  # Shape: [nPixels, nSources]
                
                # Compute all Gaussian values at once
                # all_gaussian_values = amplitudes * torch.exp(y_term + x_term)  # Shape: [nPixels, nSources]
                all_gaussian_values = torch.exp(y_term + x_term)  # Shape: [nPixels, nSources]
                
                # Compute sparse matrix multiplication for all sources at once
                total_fit = torch.mm(H.to_dense(), all_gaussian_values)

                total_fit = torch.clamp(total_fit, min=1e-4) - 1e-4
                norms = torch.norm(total_fit, dim=0, keepdim=True)
                total_fit = torch.where(norms > 0, total_fit / norms, total_fit)
                
                # Compute loss across all sources
                support = torch.logical_or(X_tensor[:,:-1] > 0, total_fit > 0)
                loss = torch.mean((total_fit[support] - X_tensor[:,:-1][support]) ** 2)
                losses.append(loss.item())
                
                # Backward pass: compute gradients
                loss.backward()
                
                # Update parameters using Adam optimizer
                optimizer.step()
                
                # Apply constraints after optimizer step
                with torch.no_grad():
                    # Ensure amplitude and sigmas are positive
                    # optim_params[:,0].clamp_(min=1e-5)
                    optim_params[:,2].clamp_(min=0.01, max=0.5)
                    optim_params[:,3].clamp_(min=0.01)

                # Check for convergence
                if epoch > 0 and abs(losses[-1] - losses[-2]) < gd_tol:
                    print(f"Converged at epoch {epoch+1} with loss: {loss.item():.8f}")
                    break
                else:
                    with torch.no_grad():
                        optim_params[norms.squeeze() == 0,2] = optim_params[norms.squeeze() == 0,2] * 5
                        # optim_params[norms.squeeze() == 0,3] = optim_params[norms.squeeze() == 0,3] * 2
                
                # Print progress
                # if (epoch + 1) % 10 == 0:
                #     print(f"Epoch {epoch+1}/{num_epochs}, Loss: {loss.item():.8f}")
            
            # Update parameters
            source_params[:,1:] = optim_params.detach()
            
            # Update A_final with the fitted Gaussians using vectorized operations
            y, x = pixel_coords_tensor[:, 1], pixel_coords_tensor[:, 2]
            y_means = source_params[:, 1].unsqueeze(0) * 10  # Shape: [1, nSources]
            x_means = source_params[:, 2].unsqueeze(0) * 10  # Shape: [1, nSources]
            y_sigmas = source_params[:, 3].unsqueeze(0) * 10  # Shape: [1, nSources]
            x_sigmas = source_params[:, 4].unsqueeze(0) * 10  # Shape: [1, nSources]
            # amplitudes = source_params[:, 0].unsqueeze(0)  # Shape: [1, nSources]

            # Compute terms for all sources at once
            y_term = -0.5 * ((y.unsqueeze(1) - y_means) / y_sigmas) ** 2  # Shape: [nPixels, nSources]
            x_term = -0.5 * ((x.unsqueeze(1) - x_means) / x_sigmas) ** 2  # Shape: [nPixels, nSources]

            # Compute all Gaussian values at once
            A_final = torch.zeros((nPixels, nSources), dtype=torch.float32)
            # A_final[selPixIdxs, :] = (amplitudes * torch.exp(y_term + x_term))
            A_final[selPixIdxs, :] = torch.exp(y_term + x_term)
            
            print("Gaussian fitting complete")

            X_update = torch.sparse.mm(H, A_final[selPixIdxs,:])
            X_update = torch.clamp(X_update, min=1e-4) - 1e-4
            norms = torch.norm(X_update, dim=0, keepdim=True)
            X_update = torch.where(norms > 0, X_update / norms, X_update)

            # Remove sources with all zeros
            non_zero_sources = torch.any(X_update > 0, dim=0)
            X_update = X_update[:, non_zero_sources]
            source_params = source_params[non_zero_sources]
            nSources = torch.sum(non_zero_sources).item()

            X_tensor = torch.cat((X_update, background_spatial_component.unsqueeze(-1)),dim=1)
            phi_lowRes = torch.cat((phi_lowRes[:,non_zero_sources.nonzero()[:,0]], phi_lowRes[:,-1].unsqueeze(-1)),dim=1)


        phi_lowRes = torch.zeros(lowResData.shape[1], nSources+1, dtype=torch.float32)
        phi_lowRes[:] = np.nan
        
        for motion_idx in motIndsToKeep:

            # Extract data for the most common motion mode
            motion_frames = (motInds == motion_idx).nonzero()[0]
            
            data_tensor = torch.from_numpy(data_for_nmf[:, motion_frames].astype(np.float32))

            sparseHIndsShifted = sparseHInds.copy()
            sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[motion_idx,0].astype(int) * dmdPixelsPerRow + uniqueMotion[motion_idx,1].astype(int)
            sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
            sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
            H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)

            # project image space (A) into superpixel space (X)
            X = torch.sparse.mm(H, A_final[selPixIdxs,:])
            X = torch.clamp(X, min=1e-4) - 1e-4

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
            phi_lowRes[motion_frames,:] = torch.linalg.solve(
                regularized_XtX,
                Xtd
            ).T
                
            # Ensure non-negativity
            phi_lowRes[motion_frames,:] = torch.clamp(phi_lowRes[motion_frames,:], min=0)
                
            # # Normalize X and phi to avoid scaling ambiguity
            # for s in range(nSources):
            #     norm = torch.norm(X_tensor[:, s])
            #     if norm > 0:
            #         X_tensor[:, s] = X_tensor[:, s] / norm
            #         phi_lowRes[motion_frames, s] = phi_lowRes[motion_frames, s] * norm
            
            error_values = []
            prev_reconstruction_error = float('inf')
            for iter_idx in tqdm(range(mult_nmf_max_iters),desc='Multiplicative NMF'):
                # Multiplicative update for NMF
                # Update phi (temporal components) using multiplicative update rule
                # phi = phi * (X^T * data) / (X^T * X * phi + epsilon)
                numerator = X.T @ data_tensor
                denominator = (X.T @ X) @ phi_lowRes[motion_frames,:].T + 1e-10
                phi_update = phi_lowRes[motion_frames,:] * (numerator / denominator).T
                phi_lowRes[motion_frames,:] = phi_update
                
                # # Normalize X and phi to avoid scaling ambiguity
                # for s in range(nSources):
                #     norm = torch.norm(X_tensor[:, s])
                #     if norm > 0:
                #         X_tensor[:, s] = X_tensor[:, s] / norm
                #         phi_lowRes[motion_frames, s] = phi_lowRes[motion_frames, s] * norm
                
                # Calculate reconstruction error
                reconstruction = X @ phi_lowRes[motion_frames,:].T
                current_error = torch.mean((data_tensor - reconstruction)**2).item()
                error_values.append(current_error)
                
                # Check convergence
                if abs(prev_reconstruction_error - current_error) < nmf_tol:
                    break
                    
                prev_reconstruction_error = current_error

        # Save all peak data to a single file
        source_extraction_data = {
            'phi_lowRes': phi_lowRes,
            'source_params': source_params,
            'A_final': A_final,
            'selPixIdxs': selPixIdxs,
        }
        output_filename = os.path.join(params['savedr'], f'source_extraction_data_DMD{DMDix+1}.npz')
        np.savez(output_filename, **source_extraction_data)
        print(f"Saved source extraction data to {output_filename}")














        # # Update A matrix based on the NMF results
        
        # # Convert sparse matrix X back to dense A without converting H to dense
        # # and without using a loop

        # sparseHIndsShiftedSelPix = sparseHIndsShifted.copy()
        # sparseHIndsShiftedSelPix[1] = np.searchsorted(selPixIdxs,sparseHIndsShifted[1])
        # H = torch.sparse_coo_tensor(sparseHIndsShiftedSelPix,sparseHVals,(numSuperPixels,selPixIdxs.shape[0]),dtype=torch.float32)
        
        # # Use sparse matrix operations to solve the system H*A = X_tensor
        # # First compute H^T * X_tensor for all sources at once
        # HtX = torch.sparse.mm(H.t(), X_tensor)
        
        # # Then compute H^T * H (still sparse)
        # HtH = torch.sparse.mm(H.t(), H)

        # # Add small regularization to ensure stability
        # HtH_reg = 1 * torch.eye(HtH.shape[0]) + HtH

        # # Solve the normal equations (H^T * H) * A = H^T * X for all sources at once
        # # using a direct solver that works with dense matrices
        # A_final[selPixIdxs,:] = torch.linalg.solve(HtH_reg, HtX)
        
        # # Ensure non-negativity
        # A_final = torch.clamp(A_final, min=0)

        # X_error_values = []
        # for iter_idx in tqdm(range(100)):
        #     numerator = torch.sparse.mm(H.t(), X_tensor)
        #     denominator = torch.sparse.mm(H.t(), H) @ A_final[selPixIdxs,:] + 1e-10
        #     A_update = A_final[selPixIdxs,:] * (numerator / denominator)
        #     A_final[selPixIdxs,:] = A_update

        #     reconstruction = torch.sparse.mm(H, A_final[selPixIdxs,:])
        #     X_error_values.append(torch.mean((X_tensor - reconstruction)**2).item())





        # refPixs = torch.from_numpy(subsampleMatrixInds[0])

        # d = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor').int()
        # c = torch.div((refPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor').int()
        # r = (refPixs % dmdPixelsPerColumn).int()
        # d = torch.from_numpy(fastZ2RefZ[d.numpy()].flatten()-1)

        # mIdxs = motInds[framesToUse]

        # for m in np.unique(mIdxs):
        #     i = np.where(mIdxs == m)[0]
        #     t = framesToUse[i]

        #     newR = r+uniqueMotion[m,0].astype(int)
        #     newC = c+uniqueMotion[m,1].astype(int)
        #     newD = d+uniqueMotion[m,2].astype(int)

        #     sparseHIndsShifted = sparseHInds.copy()
        #     sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[m,0].astype(int) * dmdPixelsPerRow + uniqueMotion[m,1].astype(int)

        #     H = torch.sparse_coo_tensor(sparseHIndsShifted,sparseHVals,(numSuperPixels,nPixels),dtype=torch.float32)

        #     X = torch.sparse.mm(H, A)



        # phi_lowRes = torch.rand(lowResData.shape[1], nSources) * 0.1  # Initialize random traces
        # bg_brightness = torch.from_numpy(lowResBrightness)

        # # Create parameters to optimize
        # source_params.requires_grad = True
        # phi_lowRes.requires_grad = True
        # bg_brightness.requires_grad = True
        
        # # Set up optimizer with all parameters to be optimized
        # lr = 0.01  # Learning rate
        # optimizer = torch.optim.Adam([source_params, phi_lowRes, bg_brightness], lr=lr)
        
        # prev_error = float('inf')

        # max_iter = 1000
        # tol = 1e-5
        # for i in range(max_iter):
        #     # Create A matrix based on source parameters (x, y, x_sigma, y_sigma)
        #     A = torch.zeros((nPixels, nSources))
        #     for i in range(nSources):
        #         # Extract parameters for each source
        #         y_pos, x_pos = source_params[i, 0], source_params[i, 1]  # x,y positions
        #         y_sigma, x_sigma = source_params[i, 2], source_params[i, 3]  # x,y sigmas
        #         amplitude = source_params[i, 4]
                
        #         # Define patch bounds (3 sigma radius)
        #         patch_radius = 3
        #         r_min = max(0, int(y_pos.detach().item() - patch_radius * y_sigma.detach().item()))
        #         r_max = min(dmdPixelsPerColumn, int(y_pos.detach().item() + patch_radius * y_sigma.detach().item()) + 1)
        #         c_min = max(0, int(x_pos.detach().item() - patch_radius * x_sigma.detach().item()))
        #         c_max = min(dmdPixelsPerRow, int(x_pos.detach().item() + patch_radius * x_sigma.detach().item()) + 1)
                
        #         # Create grid for the Gaussian using PyTorch
        #         y = torch.arange(r_min, r_max, dtype=torch.float32)
        #         x = torch.arange(c_min, c_max, dtype=torch.float32)
        #         Y, X = torch.meshgrid(y, x, indexing='ij')

        #         # Calculate Gaussian values using PyTorch operations
        #         x_diff = X - x_pos
        #         y_diff = Y - y_pos

        #         # Compute the Gaussian using PyTorch operations
        #         gaussian_patch = torch.exp(
        #             -0.5 * (
        #                 (x_diff ** 2) / (x_sigma ** 2) +
        #                 (y_diff ** 2) / (y_sigma ** 2)
        #             )
        #         )
                
        #         # Normalize the Gaussian
        #         if gaussian_patch.sum() > 0:
        #             gaussian_patch = gaussian_patch / gaussian_patch.sum() * amplitude
                
        #         # Place the Gaussian in the A matrix
        #         for z in range(numFastZs):
        #             for r_idx, r in enumerate(range(r_min, r_max)):
        #                 for c_idx, c in enumerate(range(c_min, c_max)):
        #                     pixel_idx = z * dmdPixelsPerColumn * dmdPixelsPerRow + r * dmdPixelsPerRow + c
        #                     if 0 <= pixel_idx < nPixels:
        #                         A[pixel_idx, i] = gaussian_patch[r_idx, c_idx]

        #     # Reconstruct the low-resolution data
        #     dataEst_lowRes = reconstruct.reconstruct(A, phi_lowRes, refStack[f'DMD{DMDix+1}'][0], bg_brightness,subsampleMatrixInds.T,fastZ2RefZ[f'DMD{DMDix+1}'],sparseHInds,sparseHVals,uniqueMotion,motInds)

        #     loss = torch.mean((dataEst_lowRes - torch.from_numpy(lowResData))**2)

        #     optimizer.zero_grad()
        #     loss.backward()
        #     optimizer.step()

        #     error = loss.item()
        #     if abs(prev_error - error) < tol:
        #         break

        #     prev_error = error

        #     if i % 10 == 0:
        #         print(f"Iteration {i}, Error: {error:.6f}")

        # print(f"Converged in {i} iterations")
        

if __name__ == '__main__':
    main()