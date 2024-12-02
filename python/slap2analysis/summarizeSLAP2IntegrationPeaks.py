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
from scipy import sparse, signal
import h5py
import napari
import importlib

from skimage import io as skimio, restoration
from sklearn.cluster import DBSCAN, OPTICS, HDBSCAN

sys.path.append('C:\\Users\\michael.xie\\Documents\\SLAP2_Utils\\')
from slap2_utils.datafile import DataFile
from slap2_utils.subclasses.metadata import MetaData
from slap2_utils.utils.trace import Trace, TracePixel

sys.path.append('C:\\Users\\michael.xie\\Documents\\ophys-slap2-analysis\\python')
import reconstruct
import deconvlucy

# load SLAP2 data folder
# dr = filedialog.askdirectory(initialdir = 'Z:\\scratch\\ophys\\Michael', \
#                                     title = "Select data directory")
dr = 'Z:\\scratch\\ophys\\Michael\\slap2_integration+raster\\slap2_760268_2024-11-05_12-35-49\\fov1\\experiment1'
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
params['analyzeHz'] = 1/aData['frametime'][0,0]  # analyze conventional recordings at acquisition framerate

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
    
    # allSuperPixelIDs
    allSuperPixelIDs_refs = lt['allSuperPixelIDs'][:].flat
    allSuperPixelIDs = {'DMD1': refs[allSuperPixelIDs_refs[0]][:].T,
                        'DMD2': refs[allSuperPixelIDs_refs[1]][:].T}
    
    # sparseMaskInds
    sparseMaskInds_refs = lt['sparseMaskInds'][:].flat
    sparseMaskInds = {'DMD1': refs[sparseMaskInds_refs[0]][:].T,
                      'DMD2': refs[sparseMaskInds_refs[1]][:].T}

del sparseMaskInds_refs, allSuperPixelIDs_refs

# Print shapes to verify
print("Shapes:")
for DMDix in range(2):
    print(f"allSuperPixelIDs DMD{DMDix+1}: {allSuperPixelIDs['DMD'+str(DMDix+1)].shape}")
    print(f"sparseMaskInds DMD{DMDix+1}: {sparseMaskInds['DMD'+str(DMDix+1)].shape}")

refStack = {}
for DMDix in range(2):
    pattern = f"**/*DMD{DMDix+1}_CONFIG2-REFERENCE*"
    matching_files = list(Path(dr).glob(pattern))
    if matching_files:
        first_file = str(matching_files[0])
        print(f"DMD{DMDix+1} ref stack file: {first_file}")
        numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
        refStackTmp = skimio.imread(first_file) / 100
        refStackTmp = refStackTmp.reshape(-1, numChannels, refStackTmp.shape[1], refStackTmp.shape[2]).transpose(1,0,2,3)
        refStack['DMD'+str(DMDix+1)] = refStackTmp
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
    for DMDix in range(2):
        print(f"Calculating PSF for DMD{DMDix+1}")
        numChannels = len(trialTable['refStack'][0,DMDix]['channels'][0,0].T)
        IM_dil5 = trialTable['refStack'][0,DMDix]['IM'][0,0].transpose(2,1,0)
        IM_dil5 = IM_dil5.reshape(-1, numChannels, IM_dil5.shape[1], IM_dil5.shape[2]).transpose(1,0,2,3)

        tmp = deconvlucy.deconvlucy(refStack[f'DMD{DMDix+1}'][0], IM_dil5[0], num_iter=100)
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
        from scipy.optimize import curve_fit
        popt, _ = curve_fit(lambda xy, *params: gaussian2d(xy, *params).ravel(),
                            xy, tmp_flat, p0=initial_guess)

        # Store fitted parameters
        amplitude, x0, y0, sigma_x, sigma_y = popt
        tmp = tmp[round(y0)-round(2.5*sigma_y):round(y0)+round(2.5*sigma_y)+1, round(x0)-round(2.5*sigma_x):round(x0)+round(2.5*sigma_x)+1]

        psf[f'DMD{DMDix+1}'] = (tmp+tmp[-1::-1,:]+tmp[:,::-1]+tmp[-1::-1,::-1])/4

    # save as two channel tiff
    # Create 2-channel array for PSF
    psf_combined = np.zeros((2, psf['DMD1'].shape[0], psf['DMD1'].shape[1]))
    psf_combined[0] = psf['DMD1']
    psf_combined[1] = psf['DMD2']
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

    dmdPixelsPerColumn = refStack[f'DMD{DMDix+1}'].shape[1]
    dmdPixelsPerRow = refStack[f'DMD{DMDix+1}'].shape[2]
    numRefStackZs = refStack[f'DMD{DMDix+1}'].shape[0]
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
        tmpMap[refR[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[0]//2+1,refC[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

        sparseHInds[0,spIdx*filterSize:(spIdx+1)*filterSize] = subsampleMatrixInds[spIdx,1] - 1
        sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize] = np.where(tmpMap.flatten() > 0)[0]

        sparseHVals[spIdx*filterSize:(spIdx+1)*filterSize] = tmpMap.flatten()[sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize].astype(np.uint32)]

    convMatrix = np.zeros((numSuperPixels,numSuperPixels))

    for spIdx in range(subsampleMatrixInds.shape[0]):
        tmpMap = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow))
        tmpMap[refR[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[0]//2+1,refC[spIdx].int()-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx].int()+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

        convMatrix[spIdx,:] = tmpMap[refR.int(),refC.int()]


    for trialIx in range(nTrials):
        if keepTrials[DMDix,trialIx]:
            aData = spio.loadmat(trialTable['fnAdataInt'][DMDix,trialIx][0])['aData'][0,0]

            numCycles = aData['motionDSr'].shape[0]
            uniqueMotion, motInds = np.unique(np.round(np.concatenate((aData['motionDSr'],aData['motionDSc'],aData['motionDSz']),axis=1)),axis=0,return_inverse=True)

            Afinal = torch.empty((nPixels,0))
            phiFinal = torch.empty((numCycles,0))

            importlib.reload(reconstruct)

            startTime = time.time()
            baseline = reconstruct.reconstruct(Afinal,phiFinal,refStack[f'DMD{DMDix+1}'][0],torch.from_numpy(aData['brightnessDS']),subsampleMatrixInds.T,sparseHInds,sparseHVals,uniqueMotion,motInds).detach().numpy()
            endTime = time.time()

            print(f"Elapsed Time: {endTime - startTime} sec")
