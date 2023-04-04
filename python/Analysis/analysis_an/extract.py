from __future__ import division
from read_intan_CLP_file.read_intan_CLP_file import read_intan_CLP_file
import numpy as np
from scipy import signal as scipysig
import os, warnings, json
from lab.classes.pklExperiment import pklExperiment
import process as proc
import cPickle as pkl
import pandas as pd
pd.set_option("display.max_colwidth", 10000)
from collections import OrderedDict, Iterable

def get_immediate_subdirectories(a_dir):
    return [name for name in os.listdir(a_dir)
            if os.path.isdir(os.path.join(a_dir, name))]

def get_intan_AI_during_DI_transitions(intan_clp_pkl_fpathname, DI_ch_idx = 1):
    """
    Returns raw Intan AI data between the first and last DI transition.
    
    Parameters
    ----------
    intan_clp_pkl_fpathname : str
        Path to extracted intan clamp .pkl signals file.
        
    DI_ch_idx : int
        1-based Intan digital input channel index used to define the start and end of the trace.
        
    Returns
    -------
    out : dict with keys
        fs : float
            Raw Intan sampling frequency in [Hz].
            
        sig : numpy.ndarray
            Sampled AI signal.
    """
    out = {}
    # try to import previous .pkl AI waveforms
    with open(intan_clp_pkl_fpathname, 'rb') as f:    
        intan_data = pkl.load(f)

    intan_AI = intan_data['sig']['ai']
    transition_idx = intan_data['sig']['di_transitions'][DI_ch_idx]*intan_data['fs']
    
    out['fs'] = intan_data['fs']
    out['sig'] = intan_AI[int(transition_idx[0]):int(transition_idx[-1])]
    
    return out 
    
def get_intan_AI_at_DI_transition(intan_pkl_fpathname, DI_ch_idx, periodic_sampling = True, known_transition_frequency = 0,
    timing_backup_DI_ch_idx = None, timing_backup_DI_ch_offset = 0, crop = (0, None)):

    """
    Samples Intan clamp and RHD2000 AI channel when a low-to-high transition occurs on the selected DI channel. Both AI and DI channels are assumed to be digitized by the same DAQ at the same time.
    
    Parameters
    ----------
    intan_pkl_fpathname : str
        Path to extracted intan clamp or RHD .pkl signals file.
        
    DI_ch_idx : int
        Intan digital input channel index. For clamp system use a 1-based indexing, while for RHD2000 system use a 0-based indexing.
            
    periodic_sampling : bool
        If True, flag indicates that resampling is done periodically by the DI channel.
        
    known_transition_frequency : int, float
        If known, frequency of the digital input transitions that can be used in case transitions are missing, e.g. if sampling on the DI channel was too slow.
    
    timing_backup_DI_ch_idx : int
        Backup digital input channel from which a zero time reference can be checked and estimated for DI_ch_idx used to sample the AI channel
        in case transitions are missed.
    
    timing_backup_DI_ch_offset : float
        Zero-time reference offset in [s] between backup and AI sampling DI channels. This is equipment dependent and should be reproducible between measurements.
        While the known_transition_frequency can be used to fill-in missing transitions, it would not be possible to determine if the first transition was missed.
        Convention here is that a positive value means that the backup channel DI transitions occur delayed w.r.t. the AI sampling DI_ch_idx. 
    
    crop : tuple
        Index of DI triggers to consider when sampling AI. Defaults to all received triggers. This is useful to ignore sampling the AI channel beyond an expected number of DI triggers.
        
    Returns
    -------
    out : dict 
        "fs" : None or float
            Sampling frequency determined from the first two transitions in [Hz] if periodic_sampling = True.

        "fs_raw" : float
            Original AI sampling frequency.
            
        "sig" : dict of numpy.ndarray
            Sampled AI signal at DI transitions.

        "sig_raw" : dict of numpy.ndarray
            Raw AI sampled at its original frequency, but timed to the start of the first DI transition.

        "timing_backup_DI_ch_offset" : float    
            Forwards provided parameter. See above.

        "known_transition_frequency" : float
            Forwards provided parameter. See above.

        "info" : dict with keys
            Amplifier acquisition info:
            "bandwidth": tuple
                (f_low, f_high) Low and high pass filter cutoff frequencies in [Hz].
            "dsp_filter": dict
                Digital signal processing offset removal filter.
                "cutoff_freq": float
                    DSP cutoff frequency in [Hz].
                "enabled": bool
            "ai_unit": str
                AI unit.


    """
    out = {}
    out['sig'] = {} # AI signal sampled at DI transitions
    out['sig_raw'] = {} # raw AI signal sampled at its original rate but timed to the start of the first DI transition
    out['timing_backup_DI_ch_offset'] = timing_backup_DI_ch_offset
    out['known_transition_frequency'] = known_transition_frequency

    with open(intan_pkl_fpathname, 'rb') as f:    
        intan_data = pkl.load(f)

    if type(intan_data['sig']['ai']) is dict:
        chan_labels = intan_data['sig']['ai'].keys() # multiple AI channels.
        sig = intan_data['sig']['ai']
    else:
        chan_labels = [0] # intan system has only one AI channel.
        sig = {0: intan_data['sig']['ai']}

    # ensure that all AI channels have the same number of samples
    n_ai = None
    for clabel in chan_labels:
        if n_ai == None:
            n_ai = len(sig[clabel])
        else:
            if n_ai!=len(sig[clabel]):
                raise Exception('Number of AI samples in channel ({}) is not the same as in the other channels ({}).'.format(len(sig[clabel]), n_ai))

    # check if there are DI transitions to use for sampling AI signals or if a backup timing DI channel was provided with a known transition frequency
    if not len(intan_data['sig']['di_transitions'][DI_ch_idx][crop[0]:crop[1]]):
        if not(timing_backup_DI_ch_idx != None and known_transition_frequency):
            raise Exception('Intan digital input channel {} does not have any transitions to use for sampling AI channels'.format(DI_ch_idx))

    reconstruct_transitions = False
    if len(intan_data['sig']['di_transitions'][DI_ch_idx][crop[0]:crop[1]]):
        transitions_diff = np.diff(intan_data['sig']['di_transitions'][DI_ch_idx][crop[0]:crop[1]])
        # check if indeed periodic sampling in two ways:
        # 1) if DI transition reference frequency is not known, just make sure that transitions are periodic within sampling resolution
        # 2) if DI transition reference frequency was provided, ensure that measured transition frequency is within the provided reference frequency given sampling resolution 
        if known_transition_frequency:
            uneven_samples = np.nonzero((abs(transitions_diff-1/known_transition_frequency)*intan_data['fs'])>1.1)[0] 
        else:
            uneven_samples = np.nonzero((abs(transitions_diff-transitions_diff[0])*intan_data['fs'])>1.1)[0]
            
        if len(uneven_samples) and periodic_sampling:
            # reconstruct DI transitions if backup DI transition info was given
            reconstruct_transitions = True
        else:
            out['DI_transition_idx'] = (intan_data['sig']['di_transitions'][DI_ch_idx][crop[0]:crop[1]]*intan_data['fs']).astype('int')
    else:
        # there are no DI transitions, need to reconstruct transitions
        reconstruct_transitions = True
    
    # reconstruct transitions
    if reconstruct_transitions:
        if timing_backup_DI_ch_idx != None and known_transition_frequency:
            out['DI_transition_idx'] = ((np.arange(0, int((n_ai/intan_data['fs']-(intan_data['sig']['di_transitions'][timing_backup_DI_ch_idx][0]-timing_backup_DI_ch_offset))* \
                              known_transition_frequency))/known_transition_frequency + intan_data['sig']['di_transitions'][timing_backup_DI_ch_idx][0] - \
                              timing_backup_DI_ch_offset)[crop[0]:crop[1]]*intan_data['fs']).astype(int)
        else:
            raise Exception('No backup timing DI channel was given and intan digital input {} transitions are not equally timed while periodic sampling was requested. {} uneven DI transitions found at times (s): {}\n'.\
                        format(DI_ch_idx, len(uneven_samples), uneven_samples/intan_data['fs']))
    
    if periodic_sampling and len(out['DI_transition_idx'])>=2 and not len(uneven_samples):
        out['fs'] = 1/np.mean(transitions_diff)
    elif known_transition_frequency:
        out['fs'] = known_transition_frequency  
    else:
        out['fs'] = None
          
    # downsampling anti-aliasing filter
    if out['fs'] is not None:
        # using an IIR filter is very slow, need a faster FIR filter but not sure yet what to choose
        # note: cannot use scipy.decimate because the downsampling factor must be an integer.
        # note: cannot use a single pass because of phase distortions
        cheby1_sos = scipysig.cheby1(4, 0.05, 0.8/(intan_data['fs']/out['fs']), btype = 'low', output = 'sos', analog = False)
    else:
        raise Exception('Anti-aliasing filter cannot be used for reducing sampling rate.')

    for clabel in chan_labels:
        out['sig'][clabel] = np.array([scipysig.sosfiltfilt(cheby1_sos, sig[clabel])[idx] for idx in out['DI_transition_idx']])
        out['sig_raw'][clabel] = sig[clabel][out['DI_transition_idx'][0]:out['DI_transition_idx'][-1]]

    out['fs_raw'] = intan_data['fs']
    out['info'] = {
        "bandwidth": intan_data['bandwidth'],
        "dsp_filter": intan_data['dsp_filter'],
        "ai_unit": intan_data['ai_unit']
    }
    
    return out
    
def get_imaging_ROIs(imaging_folder_path, origin = 'sima'):
    """
    Extracts region-of-interest (ROI) fluorescence data from a folder containing files following the naming convention 'ROI_<ROI name>_<channel name>'.
    
    Parameters
    ----------
    imaging_folder_path : str
        Folder path where imaging data is stored.
        
    origin : str
        Flag used for different ways in which ROI signals are structured.
    
    Returns
    -------
    dict with keys:
        'ROIs': dict of numpy.ndarray
            Fluorescence data grouped by ROI group, ROI name and by channel name.
        'fs': float
            Frame rate.
        'nframes': int
            Number of frames.
    """
    ROI_data = {}
    if origin == 'ImageJ':
        
        raise Exception('Code outdated.\n')
        
        rois_folder_path = os.path.join(imaging_folder_path, 'moco/ROI traces')
        roi_files = [f for f in os.listdir(rois_folder_path) if os.path.isfile(os.path.join(rois_folder_path, f)) and f.startswith('ROI')]
        for f_name in roi_files:
            split_f_name = f_name.split('_')
            roi_name = split_f_name[1]
            chan_name = split_f_name[2].split('.')[0]
            if f_name.endswith('.csv'):
                # comma separated txt file
                df = pd.read_csv(os.path.join(rois_folder_path, f_name))            
            elif f_name.endswith('.txt'):
                # tab separated txt file
                df = pd.read_csv(os.path.join(rois_folder_path, f_name), sep='\t')
            else:
                raise Exception('File extension import not implemented.')
            df = df.rename(columns = {'X0': 'frame', 'Y0': 'sig'})
            if roi_name not in ROI_data:
                ROI_data[roi_name] = {}
            ROI_data[roi_name][chan_name] = np.array(df.sig)
            
    elif origin == 'sima':
        
        imaging_data = {}
        imaging_data['moco_dxdy'] = {}
        
        def check_df_over_f(ROIs):
            """
            Ensure fluorescence was extracted raw instead of dF/F because this is more meaningful for later analysis
            such as background subtraction.
            """
            for roi_list_name in ROIs:
                if not 'df_over_f' in ROIs[roi_list_name]:
                    raise Exception('Cannot load fluorescence data without knowing if it is raw fluorescence or relative change dF/F. Re-extract signals_<n>.pkl')
                if ROIs[roi_list_name]['df_over_f']:
                    raise Exception('Fluorescence signals are relative, i.e. dF/F. Re-extract using raw fluorescence flag -rawf.')    
        
        ROI_data = {}
        sima_folders = [f for f in os.listdir(imaging_folder_path) if f.endswith('.sima')]

        # check consistency of imaging data source
        if not len(sima_folders):
            return {} # no ROI data 
        
        # cycle over .sima folders and use their "_<ID>.sima" suffixes to label ROI groups if more than one sima folder
        for sf in sima_folders:
            if len(sima_folders) > 1:
                sf_suffix = '_'+os.path.splitext(sf)[0].split('_')[-1]
            else:
                sf_suffix = ''

            # extract frame rate, scan start timestamp from .json file if it exists
            json_nframes = None
            json_imaging_info = [f for f in os.listdir(imaging_folder_path) if f.endswith('.json')]
            if len(json_imaging_info) == 1:
                with open(os.path.join(imaging_folder_path, json_imaging_info[0])) as f:
                    json_data = json.load(f)
                    try:
                        imaging_data['fs'] = 1/json_data['hRoiManager']['scanFramePeriod']
                    except KeyError:
                        imaging_data['fs'] = json_data['frame_rate']

                    try:
                        imaging_data['scan_start_tstamp'] = json_data['ScanStartTime']
                    except KeyError:
                        imaging_data['scan_start_tstamp'] = json_data['scan_start_tstamp']

                    try:
                        json_nframes = json_data['hStackManager']['framesPerSlice']
                    except KeyError:
                        json_nframes = None

                    try:
                        # PMT gains, 1-indexed for available channels
                        imaging_data['pmt_gains'] = {pmt_idx+1: pmt_g for pmt_idx, pmt_g in enumerate(json_data['hPmts']['gains'])}
                    except KeyError:
                        imaging_data['pmt_gains'] = json_data['pmt_gains'] # must be a dict with 'Ch<channel number>' keys and PMT gain values

                    # fluorescence digitization rate in [Hz] used to calculate photon flux correctly
                    try:
                        imaging_data['sample_rate'] = json_data['hScan2D']['sampleRate']
                    except KeyError:
                        imaging_data['sample_rate'] = json_data['pixel_sample_rate']
            
            elif len(json_imaging_info) > 1:
                raise Exception("There must be a single '.json' file in the imaging folder.")
            else:
                imaging_data['fs'] = None
                imaging_data['scan_start_tstamp'] = None
                imaging_data['pmt_gains'] = {}
                imaging_data['sample_rate'] = None
                    
            rois_folder_path = os.path.join(imaging_folder_path, sf)
            sig_files = [f for f in os.listdir(rois_folder_path) if f.startswith('signals_')]
            if not len(sig_files):
                warnings.warn("There are no imaging signal files at '{}'. Extract signals first.\n".format(rois_folder_path)) 

            # read motion correction pixel displacements
            with open(os.path.join(rois_folder_path, "sequences.pkl"), 'rb') as fh:
                seq = pkl.load(fh)
            dx = seq[0]['base']['displacements'][:,0,0]
            dy = seq[0]['base']['displacements'][:,0,1]
            # store motion correction per sima folder
            # if there is a single sima folder, sf_suffix[1:]  = '', otherwise sf_suffix[1:] = '<ID>'
            imaging_data['moco_dxdy'][sf_suffix[1:]] = np.hstack([dx[:,np.newaxis],dy[:,np.newaxis]])

            for fname in sig_files:
                # 0 - indexed channel number
                ch_num = int(os.path.splitext(fname)[0].split('_')[-1])
                ch_name = 'Ch{}'.format(ch_num+1)
                with open(os.path.join(rois_folder_path, fname), 'rb') as fh:
                    ROIs = pkl.load(fh)
                    check_df_over_f(ROIs)
                    
                # change nan pixel numbers to 0 (makes combining signals from multiple ROI patches of same name easier)
                for roi_list_name in ROIs:
                    ROIs[roi_list_name]['total_pixel_samples'][0][~np.isfinite(ROIs[roi_list_name]['total_pixel_samples'][0])] = 0

                # calculate background signal from all ROI patches that belong to either a 'background' or 'bg' labelled ROI group
                # avg_bkgrnd.shape = (nframes,)
                avg_bkgrnd = 0
                for roi_list_name in ROIs:
                    if roi_list_name in ['background', 'bg']:
                        # calculate average background level
                        avg_bkgrnd = np.nansum(ROIs[roi_list_name]['total'][0], axis = 0)/np.sum(ROIs[roi_list_name]['total_pixel_samples'][0], axis = 0)
                        # if background cannot be determined (i.e. for some frames there are no available background pixels for any), assume 0, which is less drastic compared to assuming nan
                        avg_bkgrnd[~np.isfinite(avg_bkgrnd)] = 0
                        break
                
                # calculate ROI patch signals
                for roi_list_name in ROIs:
                    roi_group_name = roi_list_name+sf_suffix
                    if roi_group_name not in ROI_data:
                        # roi_group_name = ROI group name that takes into account multiple .sima folders in a single scan (useful for independent motion correction)
                        ROI_data[roi_group_name] = {}

                    for roi_idx, roi in enumerate(ROIs[roi_list_name]['rois']):
                        # ensure that even if there is no user provided ROI label or if it was deleted, there is a valid label, which will be its ROI id
                        if roi['label']:
                            roi_name = roi['label']
                        else:
                            roi_name = roi['id']

                        if roi_name not in ROI_data[roi_group_name]:
                            ROI_data[roi_group_name][roi_name] = {}
                        
                        if roi_list_name in ['background', 'bg']:
                            # as is
                            patch_signal = ROIs[roi_list_name]['total'][0][roi_idx,:]
                        else:
                            # subtract background
                            patch_signal = ROIs[roi_list_name]['total'][0][roi_idx,:] - avg_bkgrnd*ROIs[roi_list_name]['total_pixel_samples'][0][roi_idx,:]

                        if ch_name not in ROI_data[roi_group_name][roi_name]:
                            ROI_data[roi_group_name][roi_name][ch_name] = {
                                'total_sig': patch_signal,
                                'npix_samples_sig': ROIs[roi_list_name]['total_pixel_samples'][0][roi_idx,:]
                            }
                        else:
                            ROI_data[roi_group_name][roi_name][ch_name]['total_sig'] = np.nansum([ROI_data[roi_group_name][roi_name][ch_name]['total_sig'], patch_signal], axis = 0)
                            ROI_data[roi_group_name][roi_name][ch_name]['npix_samples_sig'] += ROIs[roi_list_name]['total_pixel_samples'][0][roi_idx,:]
                        
                    # calculate average signals and tidy up dict
                    for roi_name in ROI_data[roi_group_name]:
                        ROI_data[roi_group_name][roi_name][ch_name]['total_sig'][ROI_data[roi_group_name][roi_name][ch_name]['npix_samples_sig'] == 0] = np.nan
                        ROI_data[roi_group_name][roi_name][ch_name]['avg_sig'] = ROI_data[roi_group_name][roi_name][ch_name]['total_sig']/ROI_data[roi_group_name][roi_name][ch_name]['npix_samples_sig']
                        # The total number of photons in a patch is calculated by combining the average number of photons in a DAQ sample times the median number of DAQ samples in the patch.
                        # This approach is more robust to jitter, similar to calculating the average, but here the total number of photons allows for an estimate of the patch shot noise SNR
                        ROI_data[roi_group_name][roi_name][ch_name+'_total'] = {'sig': ROI_data[roi_group_name][roi_name][ch_name]['avg_sig']*np.median(ROI_data[roi_group_name][roi_name][ch_name]['npix_samples_sig'])}
                        ROI_data[roi_group_name][roi_name][ch_name+'_avg'] = {'sig': ROI_data[roi_group_name][roi_name][ch_name]['avg_sig']}
                        del ROI_data[roi_group_name][roi_name][ch_name]

            # check if number of frames in all ROIs and channels is the same
            nframes = None
            for roi_group in ROI_data:
                for roi_name in ROI_data[roi_group]:
                    for ch_name in ROI_data[roi_group][roi_name]:
                        if nframes is None:
                            nframes = len(ROI_data[roi_group][roi_name][ch_name]['sig'])
                        elif nframes != len(ROI_data[roi_group][roi_name][ch_name]['sig']):
                            raise Exception('ROI {} from group {} for channel {} failed consistency check for equal number of frames. Found {} frames, but other ROIs have {}.'.format(
                            roi_name, roi_group, ch_name, len(ROI_data[roi_group][roi_name][ch_name]['sig']), n_frames))  
                      
            if nframes == None:
                imaging_data['nframes'] = 0
            else:
                imaging_data['nframes'] = nframes
                
            # check if there are any frames, that the number of frames is the same as in the .json file if it exists
            if json_nframes != None and json_nframes != imaging_data['nframes'] and imaging_data['nframes']:
                raise Exception('Number of frames in .json file ({}) is not the same as the number of extracted frames ({}).\n'.format(json_nframes, imaging_data['nframes']))
            
            imaging_data['ROIs'] = ROI_data  
    else:
        raise Exception('ROI format origin not implemented.')

    return imaging_data
    
def get_treadmill(behavior_filename):
    """
    Extracts raw treadmill position data.
    
    Parameters
    ----------
    behavior_filename : str
        Path to .pkl file containing behavor data. Use /code/analysis/automaticScripts/tdmlPickler.py to convert .tdml file to .pkl  
    
    Returns
    -------
    dict of 1D numpy.ndarray
        Obtains normalized belt position data as dict with labels:
        'pos_time' - Time in [s] from trial start when position is measured. Key is missing if mouse did not move
        'img_sync_time' - Time in [s] when sync TTL pulse transitions from low to high; note both position and image sync times are measured on the same board.
        'pos' - Belt position normalized to belt length. Key is missing if mouse did not move.
        'blength' - Belt length in [mm].
        'dist' - Travelled distance in [mm]. Key is missing if mouse did not move.
        'dur' - Recoding duration in [s].
    """
    bd = pklExperiment(behavior_filename).behaviorData()
    bdata = {}
    # process belt position if mouse moves or if it is present
    if 'treadmillPosition' in bd:
        # fix multiple initial 0 timestamps by truncating data
        start_idx = -1
        for idx in range(len(bd['treadmillPosition'][:, 0])):
            if bd['treadmillPosition'][idx, 0] != 0:
                break
            else:
                start_idx += 1
        assert start_idx >= 0   
    
        bdata['pos_time'] = bd['treadmillPosition'][start_idx:, 0]
        bdata['pos'] = bd['treadmillPosition'][start_idx:, 1]
        # obtain unwrapped belt distance
        # note: here using a numpy shortcut, not the most efficient but it works with minimum headache
        bdata['dist'] = np.unwrap(bdata['pos']*2*np.pi)/(2*np.pi)*bd['trackLength']
        
    bdata['blength'] = bd['trackLength']
    bdata['img_sync_time'] = bd['image_sync_pin'][:, 0]
    bdata['dur'] = bd['recordingDuration']
    
    return bdata
 
# warning: outdates as of 6/24/20
def get_intan_clp_combined_signals(folder_path, roi_label, ratio = ''):
    """
    Obtains synchronized imaging, behavior and electrophysiology data using intan timing info. Synchronization is done with respect to the imaging data.
    
    Hardware settings
    -----------------
    - Frame sync signal must be connected to Digital In 1.
    - Behavior controller heartbeat must be connected to Digital In 2.
    
    Parameters
    ----------
    folder_path : str
        Path to folder containing imaging, behavior and electrophysiology data.
        
    roi_label : str
        ROI group label.
        
    ratio : str of list of str
        If non-empty string or list of strings of the form '<Ch 1 name>/<Ch 2 name>', computes also ratio of debleached fluorescence signals in the corresponding channels.
        
    Once signals are aligned, they are saved in /analysis/aligned_signals.pkl file that will be reloaded next time. If extraction should be done again, the file should be deleted manually.
    
    Returns
    -------
    signals : dict
       Synchronized behavior treadmill position, frames and electrophysiology signals with dict having keys:
       'fs' : float
           Sampling frequency in [Hz].
       'e-phys' : dict of np.array in [mV].
           Electrophysiology signals e.g. 'LFP' for local field potential. 
    """
    
    # Intan measured signals
    n_intan_clp_files = 0
    for fname in os.listdir(folder_path):
        if fname.startswith('intan_clp'):
            intan_clamp_fpath = os.path.join(folder_path, fname)
            n_intan_clp_files += 1
    assert n_intan_clp_files == 1
    
    # extracted ROI signals folder path
    rois_folder_path = folder_path
    # behavior file
    behavior_fpath = os.path.join(folder_path, 'behavior.pkl')
    # aligned signals that will be reloaded again
    aligned_sig_fpath = os.path.join(folder_path, 'analysis/aligned_signals.pkl')
    
    # if previous aligned signals exist, reload, otherwise align again and save
    signals = {}
    if os.path.isfile(aligned_sig_fpath):
        with open(aligned_sig_fpath, 'rb') as f:
            signals = pkl.load(f)
    else:
        # convert to list of str to iterate over
        if isinstance(ratio, basestring):
            if ratio:
                ratio = [ratio]
            else:
                ratio = []
        # get imaging frame rate and local field potential signal synced to each frame
        print("Syncing LFP to frames...")
        sampled_AI = get_intan_AI_at_DI_transition(intan_clamp_fpath, 1)
        signals['fs'] = sampled_AI['fs']
        signals['LFP'] = sampled_AI['sig']
        print("done.\n")
        
        # get raw fluorescence from ROIs
        print("Getting ROI fluorescence...")
        signals['ROIs'] = get_imaging_ROIs(rois_folder_path)[roi_label]
        print("done.\n")
        
        # calculate ratiometric fluorescence signals if needed
        for r_name in ratio:
            ch_A_name, ch_B_name = r_name.split('/')
            for roi_name in signals['ROIs']:
                # add ratiometric signals to dict
                signals['ROIs'][roi_name][r_name] = proc.debleach(signals['ROIs'][roi_name][ch_A_name])/ \
                proc.debleach(signals['ROIs'][roi_name][ch_B_name])
        
        # get treadmill data
        if os.path.isfile(behavior_fpath):
            print("Getting behavior data...")
            bd = get_treadmill(behavior_fpath, smooth_velocity = 10)
            print("done.\n")
            
            # get intan DI transition indices
            print("Getting Intan digital sync data...")
            DI_transition_idx = get_intan_DI_transitions(intan_clp_AUX_filename, {'frame': 1, 'behavior sync': 2}, False)['transitions']
            print("done.\n")
            
            # ==== interpolate position data at frame times ====
            # convert frame times to behavior controller time reference
            # warning: this does not work if values fall outside of available behavior sync time stamps
            # warning: enforce intan recorded behavior sync times to be as many as times recorded by the behavior controller
            print("Interpolating belt position data...")
            # get frame times relative to the arduino time stamp
            frame_time_ref_to_behavior = np.interp(DI_transition_idx['frame'], DI_transition_idx['behavior sync'][:len(bd['img_sync_time'])], bd['img_sync_time'])
            # get traveled distance at frame time stamps 
            signals['dist'] = np.interp(frame_time_ref_to_behavior, bd['pos_time'], bd['dist'])
            
            # wrap distance again using the belt length
            signals['pos'] = proc.wrap_dist_on_belt(signals['dist'], bd['blength'])
            print("done.\n")
          
        # TESTING
        # saving signals does not always work
        # ====== save signals ======
        # create dict with all signals suitable for converting to a pandas data frame
        save_signals = {}
        # add fluorescence traces
        for ROI_name in signals['ROIs']:
            for chan_name in signals['ROIs'][ROI_name]:
                save_signals['ROI_'+ROI_name+'_'+chan_name] = signals['ROIs'][ROI_name][chan_name]
        # add remaining signals and data
        save_signals['fs'] = signals['fs']
        save_signals['LFP'] = signals['LFP']
        if 'dist' in signals:
          save_signals['dist'] = signals['dist']
        if 'pos' in signals:
          save_signals['pos'] = signals['pos']
          
        df = pd.DataFrame(save_signals)
        signals_folder_path = os.path.join(folder_path, 'analysis')
        if not os.path.exists(signals_folder_path):
            os.mkdir(signals_folder_path)
        #df.to_csv(extracted_sig_fname)
        
    return signals    
    
