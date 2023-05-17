from . import util
from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import butter, sosfiltfilt
import bottleneck as bn, numpy as np, pandas as pd
from oasis.functions import estimate_time_constant as oasis_estimate_time_constant, oasisAR1
#from oasis import nan_oasisAR1

# ================================================================================================
# CLASSES
# ================================================================================================


# ================================================================================================
# FUNCTIONS
# ================================================================================================

def fl_baseline(sig, med_filt_ns = 700, gauss_filt_sigma = 3):
    """
    Obtains fluorescence baseline. Fluorescence may have events and bleach.
    """
  
    nans = np.isnan(sig)
    # interpolate nans
    no_nans_sig = pd.Series(sig).interpolate(limit_direction = "both").to_numpy()
    filt_sig = median_filter(no_nans_sig, size = med_filt_ns, mode = "nearest") # width should be at least twice the width of an event
    filt_sig = gaussian_filter1d(filt_sig, sigma = gauss_filt_sigma) # reduce the magnitude of the minima, while being narrow enough to not contaminate the baseline with events
    
    # put back nans
    filt_sig[nans] = np.nan

    return filt_sig

def mov_avg(sig, n, no_delay = False):
    """
    Applies a moving average filter.
    
    Parameters
    ----------
    sig : array_like
        Input signal.
        
    n : int
        n-sample moving average.
        
    no_delay : bool
        If True, the delay introduced by the filter is compensated by shifting the result. This shift will introduce 'nan' elements.
        For the shift to be accurate, use odd n.

    Returns
    -------
    1D numpy array
        Returned array has same length as input array with np.nan padded to the beginning and end of the signal.
    """
    if n > 2:
        if no_delay:
            return np.append(np.insert((np.convolve(sig, np.ones((n,))/n, mode = 'valid')[:-int((n-1)/2)]).astype(float), [0]*int((n-1)/2), np.nan), [np.nan]*(n-1))
        else:
            return np.append(np.convolve(sig, np.ones((n,))/n, mode = 'valid'), [np.nan]*(n-1))
    elif n == 2:
        return np.append(np.convolve(sig, np.ones((n,))/n, mode = 'valid'), [np.nan]*(n-1))
    else:
        return np.copy(sig) #copying keeps the output consistent and there is no risk of altering the original data
    
def integrate(sig, n, mode = 'sum'):
    """
    Downsamples signal by summing neighboring samples.

    Parameters
    ----------
    sig : 1D numpy.ndarray
        Signal to integrate.

    n : int
        Number of samples to integrate.

    mode : str
        Integrator mode. Choose between sum, 'sum' and average, 'avg'.

    Returns
    -------
    1D numpy.ndarray
    """
    assert n>0
    if n > 1:
        n_resampled = int(len(sig)/n)
        out = np.empty((n_resampled,))
        for i in range(n_resampled):
            out[i] = bn.nansum(sig[n*i:n*(i+1)])
        if mode == 'sum':
            return out
        elif mode == 'avg':
            return out/n
        else:
            raise ValueError()
    else:
        return sig

def filt_butter_hp(sig, fs, fc, order):
    """
    High-pass butterworth filter.

    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : float
        Sampling frequency in [Hz].

    fc : float
        Cutoff frequency in [Hz].

    order : int
        Single-pass filter order. Actual order is double because of zero-phase filtering.
    """
    # bring mean to 0 and set nans to 0
    adjusted_sig = sig - bn.nanmean(sig)
    nans = np.isnan(adjusted_sig)
    adjusted_sig[nans] = 0
            
    nyq = 0.5 * fs
    # have to set padlen = 0 because end of signal is distorted otherwise
    out = sosfiltfilt(butter(order, fc / nyq, btype = 'highpass', output = 'sos'), adjusted_sig, padlen = 0)
    # put back nans
    out[nans] = np.nan
    
    return out

def bleach_fit_fl(sig, fs, exclude = [], fit_par = {}):
    """
    Fits a double exponential or natural cubic spline function to fluorescence data under bleaching conditions.
    The exponential fit is of the form:
        amp_fast*exp(t*ln(0.5)/tau_half_fast)+amp_slow*exp(t*ln(0.5)/tau_half_slow)+steady
    If bleaching is not well described by a bi-exponential function, a natural cubic spline fit is more appropriate (default).
    
    Parameters
    ----------
    sig : 1D numpy.nd array
        Fluorescence data. 
        
    fs : float
        Sampling rate in [Hz].
        
    exclude : list of 2 element iterable
        Start and end times in [s] of intervals to exclude from fitting.
        
    fit_par : dict
        Fitting parameters, dict with keys:
            'type': optional, str, default 'nat-cubic-spline'
                Fit method. Choose between 'exp' and 'nat-cubic-spline' for bi-exponential or natural cubic splines fit.
            'exp_tau_half_min': optional, float, default 0.1 [s]
                Shortest fluorescence half-life bleaching time constant expected in [s] controlling the bi-exponential fit.
            'spline_tstart': optional, float, default 0.1 [s]
                Start time in [s] to consider for geometrically placing knot points for natural cubic spline fitting. Cannot be 0.
            'spline_nknots': optional, int, default 10
                Number of cubic spline knot points to control smoothess. For typical fluorescence bleaching use 5-15 points.
            'outlier_perc" : optional, float, default 1
                Lower and upper 100-perc percentile cutoff for extreme events.
            'npass' : optional, int, default 3
                Number of times to exclude extreme events and refit. Typically 3x the function is robust to extreme events.
        
    Returns
    -------
    tuple (bleach_fit, exp_fit_par) where:
        bleach_fit : 1D numpy array
            Fit curve.
        exp_fit_par : dict
            If fit_par['type'] == 'exp', dict with keys:
                'amp_fast': fast bleaching component amplitude.
                'tau_half_fast': fast component half-life in [s].
                'amp_slow': slow bleaching component amplitude.
                'tau_half_slow': slow component half-life in [s].
                'steady': state state component amplitude.
            If fit_par['type'] == 'nat-cubic-spline' dict is empty.
    """
    # default bleach fit parameters
    util.set_default_keys({
        'type': 'nat-cubic-spline',
        'exp_tau_half_min': 0.1,
        'spline_tstart': 0.1,
        'spline_nknots': 10,
        'outlier_perc': 1,
        'npass': 3,
        }, fit_par)

    tmp_sig = np.copy(sig).astype(float)
    # exclude frames
    nan_intervals(tmp_sig, fs, exclude)
        
    for pass_idx in range(fit_par['npass']):
        if fit_par['type'] == 'exp':
            bleach_fit, exp_fit_par = _double_exp_bleach_fit(tmp_sig, fs, tau_half_min = fit_par['exp_tau_half_min'])
        elif fit_par['type'] == 'nat-cubic-spline':
            exp_fit_par = {}
            bleach_fit = _natural_cubic_spline_bleach_fit(tmp_sig, fs, tstart = fit_par['spline_tstart'], nknots = fit_par['spline_nknots'])
        else:
            raise ValueError("fit_par can be 'exp' or 'nat-cubic-spline'; it is '%s'."%fit_par['type'])

        # standardize signal to unit variance and zero mean, regardless of bleaching
        # the second step where a division by the variance is needed is because
        # the fluorescence signal mean is related to the poisson distribution mean by a multiplicative factor
        sig_z_score = (tmp_sig-bleach_fit)/bleach_fit**0.5
        sig_z_score /= bn.nanvar((sig_z_score))**0.5
        # get lower and upper percentiles and use them to exclude extremes to improve exponential fit
        sig_perc = np.nanpercentile(sig_z_score, [fit_par['outlier_perc'],100-fit_par['outlier_perc']])
        tmp_sig[np.less(sig_z_score,sig_perc[0]) + np.greater(sig_z_score,sig_perc[1])] = np.nan
    
    # exclude frames from fit as well (fit can become quickly innacurate beyond regions where no data is available)
    for ex in exclude:
        if ex:
            if ex[1] is not None:
                bleach_fit[int(ex[0]*fs):int(ex[1]*fs)] = np.nan
            else:
                bleach_fit[int(ex[0]*fs):] = np.nan
        
    return bleach_fit, exp_fit_par

def nan_intervals(sig, fs, intervals):
    """
    Sets given intervals to nan in the given array.

    Parameters
    ----------
    sig : 1D numpy array
        Signal.

    fs : float
        Sampling frequency in [Hz].

    intervals : list of 2 element tuples
        (start, end) intervals in [s] to be set to nan.

    Returns
    -------
    None
    """
    for interval in intervals:
        if interval:
            if interval[1] is not None:
                sig[max(0,int(interval[0]*fs)):min(len(sig),int(interval[1]*fs))] = np.nan
            else:
                sig[max(0,int(interval[0]*fs)):] = np.nan

def denoise_fl_events(sig, fs, tau_decay = None):
    """
    Denoises fluorescence traces that contain single decay time constant events. If not provided,
    event time constant is automatically estimated. The provided fluorescence trace should have a stationary baseline.

    Parameters
    ----------
    sig : 1D np.array
        Fluorescence signal.
    fs : float
        Sampling rate in [Hz].
    tau_decay : None or float
        Event decay time constant in [s].

    Returns
    -------
    tuple
        Denoised fluorescence trace and estimated or provided decay time constant (denoised_fl, tau_decay).
    """
    nan_mask = ~np.isnan(sig)
    no_nan_sig = sig[nan_mask]
    if tau_decay is None:
        g = oasis_estimate_time_constant(no_nan_sig, p = 1)[0]
    else:
        g = np.exp(-1/(tau_decay*fs))

    tau_decay = -1/fs/(np.log(g))

    denoised_fl_,_ = oasisAR1(no_nan_sig, g = g)
    denoised_fl = np.zeros_like(sig)
    denoised_fl[nan_mask] = denoised_fl_
    return denoised_fl, tau_decay
