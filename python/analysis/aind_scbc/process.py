from scipy.ndimage import gaussian_filter1d, median_filter
from scipy.signal import butter, sosfiltfilt

def bleach_fit1(sig, fs):
    smoothed_sig = gaussian_filter1d(sig, sigma = 3)
    med_filt_sig = median_filter(smoothed_sig, size = 1085, mode = "nearest")
   
    for _ in range(3):
        fc = 0.1 # cutoff frequency in [Hz]
        nyq = 0.5 * fs
        med_filt_sig = sosfiltfilt(butter(2, fc / nyq, btype = 'lowpass', output = 'sos'), np.minimum(smoothed_sig, med_filt_sig), padlen = 0)

    return med_filt_sig

def bleach_fit2(sig, fs, exclude):
    bleach_fit_par = {
        "type": "nat-cubic-spline",
        "spline_tstart": 0.1,
        "spline_nknots": 10
    }
    # bleaching fit function
    blfit, _ = an_proc.bleach_fit_fl(sig = sig, fs = fs, exclude = exclude, fit_par = bleach_fit_par)
    return blfit

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
