import os
import numpy as np
import torch
import h5py
import scipy.io as spio
import datetime
from pathlib import Path
import skimage.io as skimio
from scipy import stats
from scipy.optimize import curve_fit
from scipy.ndimage import maximum_filter, generate_binary_structure
from scipy import signal
import napari
from PyQt5.QtWidgets import QApplication, QMessageBox, QPushButton
from tkinter import filedialog
from tqdm import tqdm
import time

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
    params['synDensity'] = 2 # synapses per micron

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
            psf[f'DMD{DMDix+1}'] = skimio.imread(os.path.join(str(Path(__file__).parent), 'psfs', 'dil-17.tif'))
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

        convMatrix = np.zeros((numSuperPixels,numSuperPixels))

        for spIdx in range(subsampleMatrixInds.shape[0]):
            tmpMap = np.zeros((numFastZs,dmdPixelsPerColumn,dmdPixelsPerRow))
            tmpMap[refD[spIdx],
                    refR[spIdx]-psf[f'DMD{DMDix+1}'].shape[0]//2:refR[spIdx]+psf[f'DMD{DMDix+1}'].shape[0]//2+1,
                    refC[spIdx]-psf[f'DMD{DMDix+1}'].shape[1]//2:refC[spIdx]+psf[f'DMD{DMDix+1}'].shape[1]//2+1] = torch.from_numpy(psf[f'DMD{DMDix+1}'])

            convMatrix[spIdx,:] = tmpMap[refD,refR,refC]

        output_filename = os.path.join(params['savedr'], f'peak_data_DMD{DMDix+1}.npz')
        peaks_data = np.load(output_filename)

        peaksPix = peaks_data['peaksPix']
        peaksSigmas = peaks_data['peaksSigmas']
        peakVals = peaks_data['finalPeakVals']
        # peakTrials = peaks_data['peakTrials']

        kde = stats.gaussian_kde(peakVals)
        x_range = np.linspace(peakVals.min(), peakVals.max(), 1000)

        mode_idx = np.argmax(kde(x_range))
        bottom_half = peakVals[peakVals < x_range[mode_idx]]
        bottom_half_x_range = np.linspace(bottom_half.min(), bottom_half.max(), 500)

        def truncated_normal(x, mu, sigma, K):
            return K * stats.truncnorm.pdf(x, -np.inf, 0, loc=mu, scale=sigma)

        # print('fitting peak value pdf...')
        popt, _ = curve_fit(lambda x, mu, sigma, K: truncated_normal(x, mu, sigma, K), bottom_half_x_range, kde(bottom_half_x_range), p0=[x_range[mode_idx], np.std(bottom_half), 1])

        fit_mu = popt[0]
        fit_sigma = popt[1]
        fit_K = popt[2]

        noise_dist = stats.norm.pdf(x_range,fit_mu,fit_sigma)*fit_K*2

        signal_dist = kde(x_range) - noise_dist
        signal_dist[x_range < fit_mu] = 0
        signal_dist[signal_dist < 0] = 0

        valid_x_range = (signal_dist > max(signal_dist) * 0.1) | (noise_dist > max(noise_dist) * 0.1)
        intersection_idx = np.argmin(np.abs(noise_dist[valid_x_range] - signal_dist[valid_x_range]) / (noise_dist[valid_x_range]+1e-8))
        intersection_idx = np.where(valid_x_range)[0][intersection_idx]
        peakVal_thresh = x_range[intersection_idx]

        peaksHist = np.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))

        patch_radius = 3  # How many sigmas to include in the patch

        for i in tqdm(range(len(peaksPix))):
            if peakVals[i] > peakVal_thresh:
                mu_r, mu_c = peaksPix[i][1], peaksPix[i][2]
                sigma_r, sigma_c = max(1e-3, peaksSigmas[i][0]), max(1e-3, peaksSigmas[i][1])
                if sigma_r > 1e2 or sigma_c > 1e2:
                    continue

                # Define patch bounds
                r_min = max(0, int(mu_r - patch_radius * sigma_r))
                r_max = min(dmdPixelsPerColumn, int(mu_r + patch_radius * sigma_r) + 1)
                c_min = max(0, int(mu_c - patch_radius * sigma_c))
                c_max = min(dmdPixelsPerRow, int(mu_c + patch_radius * sigma_c) + 1)

                # Create local grid
                x, y = np.mgrid[r_min:r_max, c_min:c_max]
                pos = np.dstack((x, y))
                rv = stats.multivariate_normal([mu_r, mu_c], [[sigma_r, 0], [0, sigma_c]])
                peaksHist[r_min:r_max, c_min:c_max] += rv.pdf(pos).reshape((r_max-r_min,c_max-c_min))

        # Assume peaksHist is your 2D array
        neighborhood = generate_binary_structure(2, 2)  # 8-connectivity
        local_max = (peaksHist == maximum_filter(peaksHist, footprint=neighborhood))
        # Optionally, mask out zeros or NaNs if needed:
        local_max &= (peaksHist > 0)  # or ~np.isnan(peaksHist)

        # Get coordinates of local maxima
        maxima_coords = np.argwhere(local_max)
        histPeakVals = peaksHist[maxima_coords[:, 0], maxima_coords[:, 1]]

        histPeakVals_thresh = np.percentile(histPeakVals, 100*(1-(params['synDensity']*numSuperPixels/5/4/len(histPeakVals))))

        roi_viewer = napari.view_image(trialTable['refStack'][0,DMDix]['IM'][0,0].T, colormap='cyan')

        # Create a napari viewer with the correct aspect ratio
        roi_viewer.add_image(
            peaksHist,
            name='Peaks Histogram',
            colormap='red',
            opacity=0.5
        )
        
        # todo: test and make this compatible with fastZ imaging
        # roi_viewer.add_points(np.array([[int(fastZ2RefZ[f'DMD{DMDix+1}'][int(x[0])])-1,x[1],x[2]] for x in peaksPix]), size=1, face_color='magenta', edge_color='magenta', opacity=0.3)
        roi_viewer.add_points(np.array([[x[1],x[2]] for x in peaksPix]), size=1, face_color='magenta', edge_color='magenta', opacity=0.3)

        roi_layer = roi_viewer.add_points(maxima_coords[histPeakVals > histPeakVals_thresh], name='Sources',size=3,symbol='cross',face_color='yellow',edge_color='yellow')

        # Use a list to store both ROIs and done status
        result = {'roi_seeds': None, 'done': False}

        # Check if saved ROIs exist
        masks_filename = os.path.join(params['savedr'], f'roi_seeds_DMD{DMDix+1}.npy')
        if os.path.exists(masks_filename):
            load_saved = QMessageBox.question(None, 'Question', 'Load saved ROIs?') == QMessageBox.Yes
            if load_saved:
                print('Loading saved ROIs')
                # Load saved masks
                result['roi_seeds'] = np.load(masks_filename)

                roi_layer.data = result['roi_seeds']

                result['done'] = False
                napari.utils.notifications.show_info('Loaded saved ROIs')

        if not result['done']:
            # Callback for Done button
            def on_done():
                if 'Sources' in roi_viewer.layers:
                    result['roi_seeds'] = roi_layer.data
                    result['done'] = True
                    napari.utils.notifications.show_info('Sources captured!')

            # Add Done button
            done_button = QPushButton("Done")
            done_button.clicked.connect(on_done)

            roi_viewer.window.add_dock_widget(done_button, area='right')
            napari.utils.notifications.show_info('Draw ROIs and click "Done" when ready.')

            # Keep checking until Done is clicked
            while not result['done']:
                QApplication.processEvents()
                time.sleep(0.1)

        if len(result['roi_seeds']) > 0:
            np.save(masks_filename, result['roi_seeds'])

            print(f"Saved {len(result['roi_seeds'])} sources to {masks_filename}")


if __name__ == '__main__':
    main()