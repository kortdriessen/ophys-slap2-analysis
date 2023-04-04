from __future__ import division
import os, numpy as np
from scipy import nanmean, signal as scipysig, stats
from scipy.interpolate import interp1d as scipy_interp1d
import pandas as pd
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
import mne
import sklearn
import sklearn.linear_model
from sklearn.pipeline import Pipeline as sklearn_Pipeline
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
from statsmodels import robust
import pywt
from scipy import ndimage
from lmfit import Model as lmfitModel, conf_interval as lmfit_conf_int # Non-Linear Least-Squares Minimization and Curve-Fitting
from collections import Iterable
from . import util
from scipy.io import savemat
from future.utils import iteritems
import bottleneck as bn # speeds up certain operations on numpy arrays with nan

# ================================================================================================
# LAMBDA FUNCTIONS
# ================================================================================================

# Forward and inverse Anscombe transformations for variance-stabilizing transformation that transforms 
# a random variable with a Poisson distribution into one with an approximately standard Gaussian distribution.
# Notes:
# - the approximation holds well for poisson distributions with mean >= 4
# - closed-form unbiased inverse according to Makitalo, M., & Foi, A. (2011). "A closed-form approximation of the exact unbiased inverse 
#   of the Anscombe variance-stabilizing transformation". IEEE transactions on image processing, 20(9), 2697-2698.   
fwd_anscombe = lambda x: 2*(x+3/8.)**0.5
inv_anscombe = lambda y: y**2/4.-1/8.+(3/2.)**0.5/(4.*y)-11/(8.*y**2)+(3/2.)**0.5*5/(8.*y**3)

# Alternative variance stabilization for Poisson distribution with mean x
# the transformation will have variance 1+(3/8)*x+Order(1/x^2)
# The advantage of this transform is that ratiometric imaging can be used to cancel out modulation of signals due to motion artifacts
simple_var_stab = lambda x: 2*x**0.5

# ================================================================================================
# CLASSES
# ================================================================================================

class _AbstractSpline(sklearn.base.BaseEstimator, sklearn.base.TransformerMixin):
    """Base class for all spline basis expansions."""

    def __init__(self, max=None, min=None, n_knots=None, n_params=None, knots=None):
        if knots is None:
            if not n_knots:
                n_knots = self._compute_n_knots(n_params)
            knots = np.linspace(min, max, num=(n_knots + 2))[1:-1]
            max, min = np.max(knots), np.min(knots)
        self.knots = np.asarray(knots)

    @property
    def n_knots(self):
        return len(self.knots)

    def fit(self, *args, **kwargs):
        return self

class _NaturalCubicSpline(_AbstractSpline):
    """Apply a natural cubic basis expansion to an array.
    The features created with this basis expansion can be used to fit a
    piecewise cubic function under the constraint that the fitted curve is
    linear *outside* the range of the knots..  The fitted curve is continuously
    differentiable to the second order at all of the knots.
    This transformer can be created in two ways:
      - By specifying the maximum, minimum, and number of knots.
      - By specifying the cutpoints directly.  

    If the knots are not directly specified, the resulting knots are equally
    space within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.
    Parameters
    ----------
    min: float 
        Minimum of interval containing the knots.
    max: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.
    """

    def _compute_n_knots(self, n_params):
        return n_params

    @property
    def n_params(self):
        return self.n_knots - 1

    def transform(self, X, **transform_params):
        X_spl = self._transform_array(X)
        if isinstance(X, pd.Series):
            col_names = self._make_names(X)
            X_spl = pd.DataFrame(X_spl, columns=col_names, index=X.index)
        return X_spl

    def _make_names(self, X):
        first_name = "{}_spline_linear".format(X.name)
        rest_names = ["{}_spline_{}".format(X.name, idx)
                      for idx in range(self.n_knots - 2)]
        return [first_name] + rest_names

    def _transform_array(self, X, **transform_params):
        X = X.squeeze()
        try:
            X_spl = np.zeros((X.shape[0], self.n_knots - 1))
        except IndexError: # For arrays with only one element
            X_spl = np.zeros((1, self.n_knots - 1))
        X_spl[:, 0] = X.squeeze()

        def d(knot_idx, x):
            def ppart(t): return np.maximum(0, t)

            def cube(t): return t*t*t
            numerator = (cube(ppart(x - self.knots[knot_idx]))
                         - cube(ppart(x - self.knots[self.n_knots - 1])))
            denominator = self.knots[self.n_knots - 1] - self.knots[knot_idx]
            return numerator / denominator

        for i in range(0, self.n_knots - 2):
            X_spl[:, i+1] = (d(i, X) - d(self.n_knots - 2, X)).squeeze()
        return X_spl

# ================================================================================================
# FUNCTIONS
# ================================================================================================

def _get_natural_cubic_spline_model(x, y, minval = None, maxval = None, n_knots = None, knots = None):
    """
    Get a natural cubic spline model for the data.

    For the knots, give (a) `knots` (as an array) or (b) minval, maxval and n_knots.

    If the knots are not directly specified, the resulting knots are equally
    spaced within the *interior* of (max, min).  That is, the endpoints are
    *not* included as knots.

    Parameters
    ----------
    x: np.array of float
        The input data
    y: np.array of float
        The outpur data
    minval: float 
        Minimum of interval containing the knots.
    maxval: float 
        Maximum of the interval containing the knots.
    n_knots: positive integer 
        The number of knots to create.
    knots: array or list of floats 
        The knots.

    Returns
    --------
    model: a model object
        The returned model will have following method:
        - predict(x):
            x is a numpy array. This will return the predicted y-values.
    """

    if knots is not None:
        spline = _NaturalCubicSpline(knots = knots)
    else:
        spline = _NaturalCubicSpline(max = maxval, min = minval, n_knots = n_knots)

    p = sklearn_Pipeline([
        ('nat_cubic', spline),
        ('regression', sklearn.linear_model.LinearRegression(fit_intercept = True))
    ])

    p.fit(x, y)

    return p

def export_aligned_sig(file_name, sig, keypaths, path_mrkr = '/', ignore_missing_paths = False):
    """
    Export selected aligned signals to different file types.

    Parameters
    ----------    
    fname : str
        File path name and extension. File data type is inferred from the extension. Allowed types:
            '.mat' - Matlab
            '.csv' - Comma separated text file

    sig : dict
        Aligned signals file.

    keypaths : str, iterable of str or dict
        Key paths containing a path marker, e.g. '/' used to look up the value of an item in the nested source dict.
        If dict use as:
        {
            <old keypath>: '<new name>'
        }
    
    path_mrkr : str
        Path marker for key names in the target dict.

    ignore_missing_paths : bool
        If True, broken or missing paths are not retrieved, otherwise, if False (default) an exception is generated.
    """

    selection = util.get_dpath(source = sig, keypaths = keypaths, path_mrkr = path_mrkr, ignore_broken_paths = ignore_missing_paths)
    if not selection:
        return

    fname, fext = os.path.splitext(file_name)
    
    if fext == '.mat':
        # convert '/' to '_' and make compatible with matlab
        adjusted_selection = {k.replace('/','_'):v for k,v in iteritems(selection)}
        savemat(file_name, adjusted_selection, do_compression = False)
    elif fext == '.csv':
        df = pd.DataFrame.from_dict(selection)
        df.to_csv(file_name)
    else:
        raise ValueError("Trying to convert to '%s' file type. Not implemented/valid."%fext)

def wavelet_psd(sig, fs, wt_freq, wt_ncycles = 5):
    """
    Calculates the wavelet PSD of a given signal.

    sig : np.array
        Signal to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    wt_freq : tuple
        Wavelet transform over [lower, upper] frequency interval and number of frequency bins as (f_lower, f_upper, nbins):
            f_lower, f_upper : float
            nbins : int
        Bin size is (upper-lower)/nbins
        
    wt_ncycles : int
        Number of cycles in a wavelet at all frequencies. This is the tradeoff between temporal and frequency estimation.
    """
    
    # wavelet frequency bins
    wt_freqs = np.linspace(*wt_freq)
    # generate wavelets
    w = mne.time_frequency.tfr.morlet(sfreq = fs, freqs = wt_freqs, n_cycles = wt_ncycles, sigma = None, zero_mean = False)
    # perform continuous wavelet transform with set of wavelets
    wt = np.reshape(mne.time_frequency.tfr.cwt(np.reshape(sig,(1, -1)), w, use_fft = True, mode = 'same', decim = 1), (len(wt_freqs), -1))
    # calculate power spectral density
    wt_psd = np.abs(wt)**2

    return wt_psd
    
def calc_inst_phase_amp_freq(sig, fs, phase_offset = 0):
    """
    Calculates instantaneous amplitude, phase and frequency of a given signal using the analytical hilbert transform.

    Parameters
    ----------
    sig : 1D numpy ndarray
        Input signal.

    fs : float
        Sampling frequency in [Hz].

    Returns
    -------
    tuple of (amp, phase, freq)
        where:
        - peak to peak amplitude of signal is 2*amp
        - phase is measured in [deg] and -180+phase_offset deg corresponds to the oscillation through.
    """
    z = scipysig.hilbert(sig)
    amp = np.abs(z)
    phase = np.angle(z)
    freq = np.diff(phase)/(2*np.pi)*fs

    return (amp, phase*180/np.pi+phase_offset, freq)
 
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
        
def wzfl_threshold(ch1, ch2, zlevel = -2.5, wlevel = 3):
    """
    Z-score fluorescence thresholding with denoising wavelet filter.
    """
    # denoise fluorescence signals using discrete wavelet transform and initialize debleached signals
    ch1_db = wfl_denoise(ch1, level = wlevel)
    ch2_db = wfl_denoise(ch2, level = wlevel)

    ch1_db = debleach(ch1_db)
    ch2_db = debleach(ch2_db)
    r = ch1_db/ch2_db
    thresh = (r-bn.nanmean(r))/bn.nanstd(r) < zlevel
    #ch1_db[thresh] = np.nan
    #ch2_db[thresh] = np.nan
    return thresh

def wfl_denoise(sig, wbasis = 'db4', level = 3):
    """
    Wavelet denoising of poisson noise limited fluorescence signals. 
    """
    coeffs = pywt.wavedec(sig, wbasis, level = level)
    coeff_arr, coeff_slices = pywt.coeffs_to_array(coeffs)
    # calculate universal threshold similiar to wavethresh R-package, using method after Donoho and Johnstone
    retained_coeff = coeff_arr[coeff_slices[0][0].stop:]
    nans, x = nan_helper(retained_coeff)
    nan_idx = x(nans)
    retained_coeff_nonans = np.delete(retained_coeff, nan_idx)

    uthresh = (2*np.log(len(sig)))**0.5*robust.mad(retained_coeff_nonans)**2
    # threshold coefficients
    for i in range(1, len(coeffs)):
        coeffs[i] = pywt.threshold(coeffs[i], uthresh, mode = 'soft')
    # reconstruct signal and apply inverse anscombe transform
    return pywt.waverec(coeffs, wbasis)

def filt_wavelet_hp(sig, fs, settings):
    """
    Wavelet based high-pass filter.
    
    Parameters
    ----------
    sig : 1D numpy array
        Signal to be filtered.
        
    fs : float
        Sampling frequency in [Hz].
        
    fc : float
        Cutoff target frequency in [Hz]. Actual cutoff will vary as level decomposition will be chosen to be as close as possible.
    
    wbasis : str
        Wavelet basis function. See PyWavelets for more options.

    "settings": dict
        Filter settings
        {
            "target_fc": float, mandatory
                Target high-pass cutoff frequency in [Hz], which will be adjusted to closes available frequency.
            "adjusted_fc": float, reserved, do not pass, added back to dict
                Adjusted high-pass cutoff frequency in [Hz].
            'wbasis' : str, optional
                Wavelet basis, default 'db4'.
        }
    """
    # set default basis
    if 'wbasis' in settings:
        wbasis = settings['wbasis']
    else:
        wbasis = 'db4'

    level = int(round(np.log2(fs/settings['target_fc'])-1))
    settings['adjusted_fc'] = fs/2**(level+1)
    coeffs = pywt.wavedec(sig, wbasis, level = level)
    coeffs[0][:] = 0 # dump detail coefficients
    return pywt.waverec(coeffs, wbasis)

def zthresh_mask(sig, zlevel, siglevel = 1, mask_dilation = 0, direction = 'falling', fs = 0, window = 0):
    """
    Generates a binary mask for crossing a z-score and optionally signal level threshold. Signal may contain nan values.
    
    sig : 1D numpy.ndarray
        Signal.
    zlevel : float
        Z-score level.
    siglevel : float
        Mean-normalized signal level to be crossed, [0,1].
    mask_dilation : int
        Numer of iterations to use for dilating the threshold-crossing mask. This is handy to blank a larger duration of a brief event.
    direction : str
        Crossing direction. Choose from 'falling' or 'rising'.
    fs : float
        Signal sampling frequency in [Hz].
    window : int or float
        Moving z-score window size. If fs is non-zero, it is duration in [s], otherwise it is the number of samples. For the last samples within
        the window size, the same z-score is kept as for the last valid windowed calculation.
    """
    if fs:
        nw = int(window*fs)
    else:
        nw = int(window)

    if not window or window and len(sig)<=nw: # if signal is shorter or equal than the window, do not use a moving window
        sig_mean = bn.nanmean(sig)
        if direction == 'falling':
            thresh = ((sig-sig_mean)/bn.nanstd(sig) < zlevel) * (sig < sig_mean*siglevel)
        elif direction == 'rising':
            thresh = ((sig-sig_mean)/bn.nanstd(sig) > zlevel) * (sig > sig_mean*siglevel)
        else:
            raise ValueError
    else:
        thresh = np.full((len(sig),),False)
        for i in range(len(sig)-nw):
            sig_mean = bn.nanmean(sig[i:i+nw])
            sig_std = bn.nanstd(sig[i:i+nw])
            if direction == 'falling':
                thresh[i] = (sig[i]-sig_mean)/sig_std < zlevel and sig[i] < sig_mean*siglevel
            elif direction == 'rising':
                thresh[i] = (sig[i]-sig_mean)/sig_std > zlevel and sig[i] > sig_mean*siglevel
            else:    
                raise ValueError
        for i in range(len(sig)-nw, len(sig)):
            if direction == 'falling':
                thresh[i] = (sig[i]-sig_mean)/sig_std < zlevel and sig[i] < sig_mean*siglevel
            elif direction == 'rising':
                thresh[i] = (sig[i]-sig_mean)/sig_std > zlevel and sig[i] > sig_mean*siglevel
            else:    
                raise ValueError

    if mask_dilation:
        return ndimage.morphology.binary_dilation(thresh, iterations = mask_dilation)
    else:
        return thresh

class FIRwin_filter:

    def __init__(self, low, high, fs, ntaps = 1000, pass_zero = False, window = 'hamming'):
        """
        Constructs a windowed finite impulse response bandpass filter.
        
        Parameters
        ----------
        low, high: float
            Low and high bandpass cuttoff frequencies in [Hz].
            
        fs : float
            Sampling frequency in [Hz].
            
        ntaps : int
            Number of taps to use.
            
        pass_zero : bool
            If True, pass DC component.
            
        window : str
            Filter window to use, see scipy.signal for more info.     
        """ 
        self.taps = scipysig.firwin(ntaps, [low, high], fs = fs, pass_zero = pass_zero, \
                    window = window, scale = False)
        self.info = {'low': low, 'high': high, 'fs': fs, 'ntaps': ntaps, 'window': window}
                
    def apply(self, sig):
        """
        Applies the designed filter to a given signal. The filter is applied in a double pass to obtain a zero-phase filter.
        
        Parameters
        ----------
        sig : array_like
        
        Returns
        -------
        array_like
            Filtered signal.
        """
        return scipysig.filtfilt(self.taps, 1, sig, method = 'gust')
        
    def plot_sine_response(self, test_freq):
        """
        Applies the filter to a pure sine wave and plots the time-domain response.

        Parameters
        ----------
        test_freq : float
            Sine wave test frequency in [Hz].
        """

        t = np.arange(-100/test_freq, 100/test_freq, 1/self.info['fs'])

        # original signal
        sig = np.sin(2*np.pi*test_freq*t)
        # windowed FIR filtered signal
        firwin_sig = self.apply(sig)

        # plots
        plt.figure(figsize = (5,3))
        # original signal
        plt.plot(t*test_freq*360, sig, '-b')
        # FIR window filtered signal
        plt.plot(t*test_freq*360, firwin_sig, '-r')
        plt.xlabel('Phase (deg)')
        plt.ylabel('Amplitude')

    def plot_response(self):
        """
        Plots single pass filter magnitude, phase and step impulse responses.
        Note that this filter is applied twice, once forward and once backward to get a zero phase filter with twice the attenuation.
        """
        plot_mfreqz(self.taps, 1, self.info['fs'])
        plot_impz(self.taps, 1, self.info['fs'])
        
class PCAFilter:
    """
    Filters signals by removing common mode artifacts e.g. due to motion using PCA. For this method to work, ensure that the signals of interest are sparse in time
    and for better results, exponential bleaching is removed prior to using this method. 
    """  
    def __init__(self, sig):
        """
        Initialized PCA filter by first computing PCA on supplied signals.
        
        Parameters
        ----------
        sig: iterable of numpy 1D ndarray
        """
        self.ssig = np.column_stack(sig)
        self.pca = sklearn.decomposition.PCA()
        self.pca.fit(self.ssig)
    
    def filt(self, n):
        """
        Filter signals by throwing out the first n PCA components. This assumes that the common mode artifact is dominant over the uncorrelated activity e.g. sparse synaptic
        activation.
        """
        return np.dot(self.pca.transform(self.ssig)[:,n:], self.pca.components_[n:,:]) + np.mean(self.ssig, axis = 0)
        

#def pca_filter(signals, ):

def tslice(sig, fs, t_slice = (0, None), return_intervals = False):
    """
    Returns a time slice of a signal.
    
    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : None or float
        If given, sampling frequency in [Hz].
        
    t_slice : tuple of int, float or None
        Signal slice start and end array indices or times if fs given as (t_start, t_end). If t_end is None, the signal starts at 0 and ends at len(sig) or t = len(sig)/fs

    return_intervals : bool
        Output format, see Returns.
    
    Returns
    -------
    sig_slice if return_intervals is False, otherwise tuple (sig_slice, (t_start, t_end))
        where sig_slice as array_like and standardized start and end array indices or times in [s] as tuple of int or float.
    """
    if len(sig):
        t_start = t_slice[0]
        if t_slice[1] is None:
            if fs is None:
                t_end = len(sig)-1
                sliced_array = sig[t_start:]
            else:
                t_end = (len(sig)-1)/fs
                sliced_array = sig[int(round(t_start*fs)):]
        else:
            if fs is None:
                t_end = min(t_slice[1], len(sig)-1)
                sliced_array = sig[t_start:t_slice[1]+1]
            else:
                t_end = min(t_slice[1],(len(sig)-1)/fs)
                sliced_array = sig[int(round(t_start*fs)):int(round(t_end*fs))+1]
    else:
        sliced_array = np.array([])
        t_start = 0
        t_end = 0

    assert t_start <= t_end

    if return_intervals:
        return sliced_array, (t_start, t_end)
    else:
        return sliced_array
    
def filt_FIR_bp(sig, fs, f_low, f_high, passband_ripple = 1e-3, stopband_attn = 1e-3, transition_bw = 1, pass_zero = False):
    """
    Double pass zero-phase FIR band-pass filter.
    
    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : float
        Sampling frequency in [Hz].
    
    f_low, f_high : float
        Band-pass low and high cut-off frequencies in [Hz].
        
    passband_ripple : float 
        Amplitude variation within passband.
        
    stopband_attn : float
        Stop band amplitude attenuation.
        
    transition_bw : float
        Transition bandwidth in [Hz], i.e. difference between end of pass-band and start of stop-band.
        
    pass_zero : bool
        Whether to pass DC component or not. 
    """
    # determine number of taps needed empirically
    ntaps = int(round(2/3. *np.log10(1/(10*passband_ripple*stopband_attn))*fs/transition_bw))
    
    if ntaps > 1000:
        print('Warning, number of taps {} is very large, this could take too long. Consider changing parameters.'.format(ntaps))
    filt = FIRwin_filter(f_low, f_high, fs, ntaps = ntaps, pass_zero = pass_zero)
    return filt.apply(sig), filt
    
def filt_FIR_hp(sig, fs, f_stop, f_pass, order = 101):
    """
    Zero-phase finite impulse response high-pass filter.

    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : float
        Sampling frequency in [Hz].

    f_stop, f_pass : float
        Transition band frequency limits in [Hz].

    order : int
        Number of taps, must be odd number.
    """
    # design high-pass filter
    f_nyquist = fs / 2.
    desired = (0, 0, 1, 1)
    bands = (0, f_stop, f_pass, f_nyquist)
    filter_coefs = scipysig.firls(order, bands, desired, nyq = f_nyquist)

    # apply high-pass filter
    return scipysig.filtfilt(filter_coefs, [1], sig)

def filt_theta_band(sig, fs, theta_band = (5, 10)):
    """
    Filters a signal to be contained within the theta band.
    
    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : float
        Sampling frequency in [Hz].
        
    theta_band : tuple of float
        Theta-band frequencies as (lower, upper).
                
    Returns
    -------
    array_like
        Filtered theta band signal.
    """
    return filt_FIR_bp(sig, fs, theta_band[0], theta_band[1], transition_bw = 1)[0]
    
def filt_ripple_band(sig, fs, passband = [80, 250]):
    """
    Filters ripple band signal using a zero-phase 2nd-order butterworth filter.

    Parameters
    ----------
    sig : array_like
        Signal.
        
    fs : float
        Sampling frequency in [Hz].

    passband : tuple of float
        Ripple pass-band frequencies as (lower, upper) in [Hz].
    """
    b, a = scipysig.butter(1, [passband[0]/(fs/2.), passband[1]/(fs/2.)], 'bandpass', analog = False) # filtfilt will double the order
    return scipysig.filtfilt(b, a, sig)

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
    out = scipysig.sosfiltfilt(scipysig.butter(order, fc / nyq, btype = 'highpass', output = 'sos'), adjusted_sig, padlen = 0)
    # put back nans
    out[nans] = np.nan
    
    return out

def sos_butter_bp(lowcut, highcut, fs, order = 5):
    """
    Designs a Butterworth band-pass filter using second-order sections.

    Parameters
    ----------
    lowcut, highcut : float
        -3dB cutoff frequencies.

    fs : float
        Sampling frequency in [Hz].

    order : int
        Filter order.

    Returns
    -------
    ndarray
        sos filter coefficients
    """
    nyq = 0.5 * fs
    low = lowcut / nyq
    high = highcut / nyq

    return scipysig.butter(order, [low, high], btype = 'band', output = 'sos')

def filt_butter_bp(sig, fs, lowcut, highcut, order = 2):
    """
    Applies a zero-phase Butterworth filter.

    Parameters
    ----------
    sig : 1D numpy array
        Signal.

    fs : float
        Sampling frequency in [Hz].

    lowcut, highcut : float
        -3dB cutoff frequencies.

    order : int
        Single-pass filter order. For this double-pass zero-phase filter, the order will be double.

    Returns
    -------
    1D numpy array
    """
    return scipysig.sosfiltfilt(sos_butter_bp(lowcut, highcut, fs, order), sig)

def filt_savgol(sig, fs, settings):
    """
    Applies a low-pass zero-phase savgol filter.

    Parameters
    ----------
    sig : 1D numpy array
        Signal.

    fs :
        Sampling frequency in [Hz].

    settings : dict
        Filter settings, dict with keys:
            "target_window": float, mandatory
                 Polynomial fitting window duration in [ms] (actual window duration is adjusted to get 
                 an odd number of samples within the window for a minimum of 3 samples).
            "adjusted_window": float, reserved, do not pass, added back to dict
                Adjusted moving average window duration in [ms].
            "order": int, mandatory
                Savgol filter polynomial fit order
            "cval" : float, optional
                If known, value to use to fill the fitting window at the edges. If not specified, edge windows are interpolated.

    Returns
    -------
    1D numpy array
        Filtered signal.
    """
    # round to closest odd number of samples >= 3
    n_savgol_window = util.round_up_to_odd_int(settings['target_window']*1e-3*fs)
    if n_savgol_window < 3:
        n_savgol_window = 3
        
    # add back adjusted window duration
    settings['adjusted_window'] = n_savgol_window*1e3/fs
    if 'cval' in settings: 
        return scipysig.savgol_filter(sig, n_savgol_window, settings['order'], mode = 'constant', cval = settings['cval'])
    else:
        return scipysig.savgol_filter(sig, n_savgol_window, settings['order'], mode = 'interp')

def filt_nat_cubic_spline(sig, fs, settings):
    """
    Applies a high- or low-pass filter by fitting a natural cubic spline.

    Parameters
    ----------
    sig : 1D numpy.nd array
        Fluorescence data.     
    fs : float
        Sampling rate in [Hz].
    settings : dict
        Filter settings, dict with keys:
        'knot_interval' : float
            Time interval between spline knots in [ms].
        'type' : str
            Filter type, choose 'hp' or 'lp' for high- or low-pass respectively. If high-pass, signal is normalized to spline fit.
        'outlier_perc" : optional, float, default 0
            Lower and upper 100-perc percentile cutoff for extreme events. Extreme events should be brief for this filter to work best.
        'npass' : optional, int, default 2
            Number of times to exclude extreme events and refit. Typically 3x the function is robust to extreme events.
            If 'outlier_perc' == 0 this is overriden and set to 1.
    debug : bool
        If True, return additional information to monitor filter performance.

    Returns
    -------
    If debug is False:
        1D numpy array
            Filtered signal.
    """
    nan_helper = lambda x: (np.isnan(x), lambda z: z.nonzero()[0])
    # default bleach fit parameters
    util.set_default_keys({
        'outlier_perc': 0,
        'npass': 2,
        }, settings)

    assert settings['npass'] > 0
    # override number of passes if there are no outliers
    if settings['outlier_perc'] == 0:
        settings['npass'] = 1

    tmp_sig = np.copy(sig).astype(float)
    outliers = np.full((len(sig),), False)
    updated_outliers = outliers
    time = np.arange(len(sig))/fs
    nknots = int(round(len(tmp_sig)/(settings['knot_interval']*1e-3*fs)))
    knots = list(np.linspace(0, time[-1], nknots))
    sig_nans = None
    # take out knots that fall in regions where the signal is not defined
    nan_knot_idx = []
    for k_idx, k in enumerate(knots):
        if np.isnan(tmp_sig[int(k*fs)]):
            nan_knot_idx.append(k_idx)
    util.del_list_idxs(knots, nan_knot_idx)
    for pass_idx in range(settings['npass']):
        outliers = updated_outliers
        # exclude NaN values in signal to be able to calculate a fit
        nans, x = nan_helper(tmp_sig)
        # retain original signal nans
        if not pass_idx:
            sig_nans = nans
        nan_idx = x(nans)
        sig_nonans = np.delete(tmp_sig, nan_idx)
        time_nonans = np.delete(time, nan_idx)
        
        model = _get_natural_cubic_spline_model(time_nonans, sig_nonans, knots = knots)
        spline_fit = model.predict(time)
        # exclude outliers
        # =======================
        if settings['outlier_perc'] != 0:
            # standardize signal to unit variance and zero mean, regardless of bleaching
            # the second step where a division by the variance is needed is because
            # the fluorescence signal mean is related to the poisson distribution mean by a multiplicative factor
            sig_z_score = (tmp_sig-spline_fit)/spline_fit**0.5
            sig_z_score /= bn.nanvar((sig_z_score))**0.5
            # get lower and upper percentiles and use them to exclude extremes to improve exponential fit
            sig_perc = np.nanpercentile(sig_z_score, [settings['outlier_perc'],100-settings['outlier_perc']])
            updated_outliers = np.logical_or(np.less(sig_z_score,sig_perc[0]), np.greater(sig_z_score,sig_perc[1]))
            tmp_sig[updated_outliers] = np.nan

    if settings['type'] == 'hp':
        return sig/spline_fit
    elif settings['type'] == 'lp':
        spline_fit[sig_nans] = np.nan
        return spline_fit
    else:
        raise ValueError("Filter type can be only 'hp' or 'lp'.")

def filter_bank(sig, fs, filters, hp_dc = 0):
    """
    Applies a series of filters to a given signal.

    Parameters
    ----------
    sig : numpy.array
        Waveform to process.
        
    fs : float
        Sampling frequency in [Hz].

    filters: list of dict
        Filters, their order to be applied, and their settings. Dict has keys:
        'type' : str
            Name of filter to apply.
        'settings' : dict with various filter-specific settings as keys and their values. The following filters are implemented:

        1. moving average, LP filter
        ----------------------------
        Moving average filter.
        {
            "type": "mov_avg",
            "settings":
            {
                "target_window": float, mandatory
                    Requested moving average window duration in [ms]. Actual window duration is rounded up to closest value to have 
                    an odd number of samples in the window and have 0 lag filtering.
                "adjusted_window": float, reserved, do not pass, added back to dict
                    Adjusted moving average window duration in [ms].
            }
        }

        2. integrator, LP filter
        ------------------------
        Integrator filter. Applies an integrator to a given signal of an integer number of samples and upsamples back to obtain
        the same number of samples as the original signal. Signal can contan NaNs.
        {
            "type": "integrator",
            "settings":
            {
                "int_fs" : float
                    Integrator sampling rate in [Hz]. Number of integrator samples is round(fs/settings['int_fs']), where fs is the sampling
                    rate of the original signal.
                'mode' : str
                    Integrator mode, choose between using an average 'avg' or sum, 'sum'.
            }
        }

        3. savgol, LP filter
        --------------------
        Savitzky-Golay low-pass filter.
        {
            "type": "savgol",
            "settings":
            {
                "target_window": float
                     Polynomial fitting window duration in [ms] (actual window duration is adjusted to get 
                     an odd number of samples within the window for a minimum of 3 samples).
                "adjusted_window": float, reserved, do not pass, added back to dict
                    Adjusted moving average window duration in [ms].
                "order": int
                    Savgol filter polynomial fit order
                "cval" : float, optional
                    If known, value to use to fill the fitting window at the edges. If not specified, edge windows are interpolated.
            }
        }

        4. nat_cubic_spline, LP filter
        ------------------------------
        Natural cubic spline filter. Fits and normalizes signal with an equally spaced knot natural cubic spline with given time interval spacing between knots.
        {
            "type": "nat_cubic_spline"
            "settings":
            {
                "type" : str
                    Filter type, choose 'hp' or 'lp' for high- or low-pass respectively. If high-pass, signal is normalized to spline fit.
                "knot_interval": float
                    Time interval between spline knots in [ms].
                'outlier_perc" : optional, float, default 1
                    Lower and upper 100-perc percentile cutoff for extreme events.
                'npass' : optional, int, default 2
                    Number of times to exclude extreme events and refit. Typically 3x the function is robust to extreme events.
                    If 'outlier_perc' == 0 this is overriden and set to 1.
            }
        }

        5. wavelet_hp, HP filter
        ------------------------
        Wavelet-based high-pass filter.
        {
            "type": "wavelet_hp",
            "settings":
            {
                "target_fc": float, mandatory
                    Target high-pass cutoff frequency in [Hz], which will be adjusted to closes available frequency.
                "adjusted_fc": float, reserved, do not pass, added back to dict
                    Adjusted high-pass cutoff frequency in [Hz].
                'wbasis' : str, optional
                    Wavelet basis, default 'db4'.
            }
        }

        6. butter_hp, HP filter
        -----------------------
        Butterworth high-pass filter. Does not work if signal has nan values, replace or cut out.
        {
            "type": "butter_hp",
            "settings":
            {
                "fc": float
                    Cutoff frequency in [Hz].
                "order": int
                    Single-pass filter order, actual order will be double (zero-phase filtering).
            }
        }

    hp_dc : float
        Adds DC offset back after applying a given HP filter.
    """
    filt_sig = sig
    for f in filters:
        if f['type'] == 'savgol':
            filt_sig = filt_savgol(sig = filt_sig, fs = fs, settings = f['settings'])
        
        elif f['type'] == 'mov_avg':
            filt_sig = filt_mov_avg(sig = filt_sig, fs = fs, settings = f['settings'])

        elif f['type'] == 'integrator':
            filt_sig = filt_integrator(sig = filt_sig, fs = fs, settings = f['settings'])

        elif f['type'] == 'nat_cubic_spline':
            filt_sig = filt_nat_cubic_spline(sig = filt_sig, fs = fs, settings = f['settings'])

        elif f['type'] == 'wavelet_hp':
            filt_sig = filt_wavelet_hp(sig = filt_sig, fs = fs, settings = f['settings'])
            filt_sig += hp_dc

        elif f['type'] == 'butter_hp':
            filt_sig = filt_butter_hp(sig = filt_sig, fs = fs, fc = f['settings']['fc'], order = f['settings']['order'])
            filt_sig += hp_dc

        else:
            raise ValueError

    return filt_sig

def filt_mov_avg(sig, fs, settings):
    """
    Zero lag moving average filter.

    Parameters
    ----------
    sig : numpy.array
        Waveform.
        
    fs : float
        Sampling frequency in [Hz].

    settings : dict
        Filter settings. Dict with keys:
            'target_window' : float
                Target moving average window. Actual window is adjusted to round number of samples to
                closest odd number >= 3.
            'adjusted_window' : float
                Actual averaging window.

    Returns
    -------
    1D numpy.ndarray
    """
    # round to closest odd number of samples >= 3
    n_mov_avg = util.round_up_to_odd_int(settings['target_window']*1e-3*fs)
    if n_mov_avg < 3:
        n_mov_avg = 3
    filt_sig = mov_avg(sig = sig, n = n_mov_avg, no_delay = True)
    # store adjusted filter settings
    settings['adjusted_window'] = n_mov_avg*1e3/fs

    return filt_sig

def filt_integrator(sig, fs, settings):
    """
    Applies an integrator to a given signal of an integer number of samples and upsamples back to obtain
    the same number of samples as the original signal. Signal can contan NaNs.

    Parameters
    ----------
    sig : numpy.array
        Waveform.
        
    fs : float
        Sampling frequency in [Hz].

    settings : dict
        Filter settings. Dict with keys:
            'int_fs' : float
                Integrator sampling rate in [Hz]. Number of integrator samples is round(fs/settings['int_fs']).
            'mode' : str
                Integrator mode, choose between using an average 'avg' or sum, 'sum'.

    Returns
    -------
    1D numpy.ndarray
    """
    int_nsamp = int(round(fs/settings['int_fs']))
    if int_nsamp < 1:
        raise Exception('Integrator target frequency {} is higher or equal to an integer multiple >= 2 of the original sampling frequency {}.'.format(
            settings['int_fs'], fs))
    elif int_nsamp == 1:
        return sig

    int_sig = scipysig.resample(integrate(sig = sig, n = int_nsamp, mode = settings['mode']), len(sig))
    int_sig[np.isnan(sig)] = np.nan
    return int_sig

def phase_bin_avg(sig, fs, phase_reset_times, nbins = 12):
    """
    Performs phase-bin averaging of a signal according to irregular defined through times.   Through times are referenced to 0 phase and through-to-through phase is 360 deg.
    
    Parameters
    ----------
    sig : numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    phase_reset_times : np.array
        Time in [s] at which phase is set to 0 deg.
        
    nbins : int
        Number of phase bins to use for averaging.
    """
    # phase bin average
    phase_bin_av = np.zeros((nbins,))
    # number of samples to average per phase bin
    phase_bin_navg = np.zeros((nbins,))
    # sort phases
    phase_reset_idx = (phase_reset_times*fs).astype('int')
    reset_idx_ctr = 0
    for sig_idx in range(phase_reset_idx[0], phase_reset_idx[-2]):
        if sig_idx >= phase_reset_idx[reset_idx_ctr+1]:
            reset_idx_ctr += 1
        bin_idx = int(nbins*(sig_idx-phase_reset_idx[reset_idx_ctr])/(phase_reset_idx[reset_idx_ctr+1] - phase_reset_idx[reset_idx_ctr]))
        phase_bin_av[bin_idx] += sig[sig_idx]
        phase_bin_navg[bin_idx] += 1
        
    return phase_bin_av/phase_bin_navg

# another way to do phase bin averaging using piecewise interpolation
# For some reason this function is an older duplicate               
"""
def phase_bin_avg2(sig, fs, phase_reset_times, nbins = 20):
    
    Parameters
    ----------
    sig : numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    phase_reset_times : np.array
        Time in [s] at which phase is set to 0 deg.
        
    nbins : int
        Number of phase bins to use for averaging.     
    
    # phase bin average
    phase_bin_av = np.zeros((nbins,))
    # number of samples to average per phase bin
    phase_bin_navg = np.zeros((nbins,))
    # calculate normalized phase vector i.e. a full cycle is mapped to [0, 1) 
    time = np.arange(len(sig))/fs
    current_phase_reset_idx = 0
    for t_idx, t in enumerate(time):
        if t < phase_reset_times[0]:
            # calculate phase assuming that the leftmost incomplete cycle has the same frequency as the first complete cycle
            phase = 1-(phase_reset_times[0]-t)/(phase_reset_times[1]-phase_reset_times[0])
        elif t >= phase_reset_times[current_phase_reset_idx]:
            if current_phase_reset_idx < len(phase_reset_times)-1:
                if t >= phase_reset_times[current_phase_reset_idx+1]:
                    current_phase_reset_idx += 1
                if current_phase_reset_idx < len(phase_reset_times)-1:
                    phase = (t-phase_reset_times[current_phase_reset_idx])/(phase_reset_times[current_phase_reset_idx+1]-phase_reset_times[current_phase_reset_idx])
                else:
                    # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                    phase = (t-phase_reset_times[-1])/(phase_reset_times[-1]-phase_reset_times[-2])
            else:
                # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                phase = (t-phase_reset_times[-1])/(phase_reset_times[-1]-phase_reset_times[-2])
        bin_idx = int(phase*nbins)
        phase_bin_av[bin_idx] += sig[t_idx]
        phase_bin_navg[bin_idx] += 1        
            
    return phase_bin_av/phase_bin_navg
"""

def assign_phase_v2(sigs, fband, filt_order = 2):
    """
    Assigns phase values to a fluorescence signal using a hilbert transform of a
    band-passed local field potential. A phase of 0 corresponds to LFP throughs.
    note: nan values are not included.

    Parameters
    ----------
    sigs : list of dict
        Signals to analyze. List of dict with keys and values:
            "fluo" : 1D numpy.array
                Fluorescence waveform to analyze.
            "lfp" : 1D numpy.array
                Local-field potential.
            "gating" : optional, 1D numpy.array
                Array marking portions of the fluorescence signal "fluo" that are included in the phase analysis. 
                True marks samples for inclusion.
            "fs" : float
                Sampling frequency in [Hz].

    fband : tuple
        Minimum and maximum frequencies of the Butterworth LFP bandpass filter.

    filt_order : int
        One way Butterworth filter order. Total order will be double as a zero-phase forward-backward filtering is used.

    Returns
    -------
    (phase, fluo, total_valid_time) : tuple

        phase : 1D np.ndarray
            Concatenated phase in [deg] assigned to each fluorescence sample in the input signals. Zero phase corresponds to LFP throughs.

        fluo : 1D np.ndarray
            Concatenated fluorescence with nan and gating values removed.

        total_valid_time : float
            Total valid time in [s] used for phase estimation, which includes valid non-nan fluorescence samples and a True-valued gate.
    """
    total_valid_time = 0
    phase = []
    fluo = []
    for sig in sigs:
        # band pass LFP signal
        lfp_bp_sig = filt_butter_bp(sig = sig["lfp"], fs = sig["fs"], lowcut = fband[0], highcut = fband[1], order = filt_order)
        _, _phase, _ = calc_inst_phase_amp_freq(sig = lfp_bp_sig, fs = sig["fs"], phase_offset = 180)
        # pick out nan and gated samples
        if "gating" in sig and len(sig["gating"]):
            valid_samples_mask = np.logical_and(np.isfinite(sig["fluo"]), sig["gating"])
        else:
            valid_samples_mask = np.isfinite(sig["fluo"])

        total_valid_time += np.sum(valid_samples_mask)/sig["fs"]
        phase.append(_phase[valid_samples_mask])
        fluo.append(sig["fluo"][valid_samples_mask])

    return np.concatenate(phase), np.concatenate(fluo), total_valid_time

def assign_phase(phase_reset_times, fs, sig, phase_offsets = [], gating = [], drop_nan = True):
    """
    Assigns phase values to fluorescence waveform using a phase reset time vector.
    note: nan values are not included.
    
    Parameters
    ----------
    phase_reset_times : list of np.array 
        Time in [s] at which phase jumps.

    fs : list of float
        Sampling frequency in [Hz].

    sig : list of numpy.array
        Waveforms to analyze.
        
    phase_offsets : list of float
        Phase offsets in [deg] to subtract from the measured phase.

    gating : list of numpy.array
        Optional list of 1D bool numpy array marking portions of the sig that are included in the phase analysis. 
        True marks samples for inclusion.

    drop_nan : bool
        If True, drop phase estimates for nan samples.

    Returns
    -------
    (phase, sig, total_valid_time) : tuple

        phase : 1D np.ndarray
            Phase in [deg] assigned to each sample in the input signals. Phase range is [0, 360) deg.

        sig : 1D np.ndarray
            Concatenated input signal with nan and gating values removed.

        total_valid_time : float
            Total valid time in [s] used for phase estimation, which includes valid non-nan fluorescence samples and a True-valued gate.
    """
    total_valid_time = 0
    sig_out = []
    phase_out = []
    assert len(phase_reset_times) == len(fs)
    if sig:
        assert len(sig) == len(phase_reset_times)
    if gating:
        assert len(gating) == len(phase_reset_times)
    if phase_offsets:
        assert len(phase_offsets) == len(phase_reset_times)
            
    for sig_idx in range(len(phase_reset_times)):
        valid_samples = 0
        time = np.arange(len(sig[sig_idx]))/fs[sig_idx]
        current_phase_reset_idx = 0
        if phase_offsets:
            phase_offset = phase_offsets[sig_idx]
        else:
            phase_offset = 0
        for t_idx, t in enumerate(time):
            if t < phase_reset_times[sig_idx][0]:
                # calculate phase assuming that the leftmost incomplete cycle has the same frequency as the first complete cycle
                phase = 1-(phase_reset_times[sig_idx][0]-t)/(phase_reset_times[sig_idx][1]-phase_reset_times[sig_idx][0])
            elif t >= phase_reset_times[sig_idx][current_phase_reset_idx]:
                if current_phase_reset_idx < len(phase_reset_times[sig_idx])-1:
                    if t >= phase_reset_times[sig_idx][current_phase_reset_idx+1]:
                        current_phase_reset_idx += 1
                    if current_phase_reset_idx < len(phase_reset_times[sig_idx])-1:
                        phase = (t-phase_reset_times[sig_idx][current_phase_reset_idx])/(phase_reset_times[sig_idx][current_phase_reset_idx+1]-phase_reset_times[sig_idx][current_phase_reset_idx])
                    else:
                        # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                        phase = (t-phase_reset_times[sig_idx][-1])/(phase_reset_times[sig_idx][-1]-phase_reset_times[sig_idx][-2])
                else:
                    # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                    phase = (t-phase_reset_times[sig_idx][-1])/(phase_reset_times[sig_idx][-1]-phase_reset_times[sig_idx][-2])
            if gating:
                if drop_nan:
                    if len(gating[sig_idx]) and np.isfinite(sig[sig_idx][t_idx]) and gating[sig_idx][t_idx]:
                        sig_out.append(sig[sig_idx][t_idx])        
                        phase_out.append(((phase%1)*360-phase_offset)%360)
                        valid_samples += 1
                else:
                    sig_out.append(sig[sig_idx][t_idx])        
                    phase_out.append(((phase%1)*360-phase_offset)%360)
                    valid_samples += 1
            else:
                if drop_nan:
                    if np.isfinite(sig[sig_idx][t_idx]):
                        sig_out.append(sig[sig_idx][t_idx])        
                        phase_out.append(((phase%1)*360-phase_offset)%360)
                        valid_samples += 1
                else:
                    sig_out.append(sig[sig_idx][t_idx])        
                    phase_out.append(((phase%1)*360-phase_offset)%360)
                    valid_samples += 1

        total_valid_time += valid_samples/fs[sig_idx]

    return np.array(phase_out), np.array(sig_out), total_valid_time

# another way to do phase bin averaging using piecewise interpolation               
def phase_bin_avg2(sig, fs, phase_reset_times, nbins = 20, alpha = 0.05, bsiter = 5000, n_threads = 20):
    """
    Calculates the phase-bin average of multiple signals and their associated phase reset time to generate an average waveform aligned
    to e.g. each theta-band LFP oscillation cycle.
    
    Parameters
    ----------
    sig : list of numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    phase_reset_times : list of np.array 
        Time in [s] at which phase is set to 0 deg.
        
    nbins : int
        Number of phase bins to use for averaging. 
        
    alpha : float
        Alpha value representing the confidence interval. Defaults to 0.05, i.e., 95th-CI.
    
    bsiter : int
        Number of bootstrap iterations to run. The higher, the more accurate confidence intervals can be estimated,
        but runs also slower.
    
    n_threads : int
        Number of threads to use to speed up boostrap estimation.
        
    """
    # each signal should have a corresponding phase reset vector
    assert len(sig) == len(phase_reset_times)
    
    # phase bin average
    phase_bins = [[] for i in range(nbins)]
    # number of samples to average per phase bin
    # calculate normalized phase vector i.e. a full cycle is mapped to [0, 1) 
    for sig_idx in range(len(sig)):
        time = np.arange(len(sig[sig_idx]))/fs
        current_phase_reset_idx = 0
        for t_idx, t in enumerate(time):
            if t < phase_reset_times[sig_idx][0]:
                # calculate phase assuming that the leftmost incomplete cycle has the same frequency as the first complete cycle
                phase = 1-(phase_reset_times[sig_idx][0]-t)/(phase_reset_times[sig_idx][1]-phase_reset_times[sig_idx][0])
            elif t >= phase_reset_times[sig_idx][current_phase_reset_idx]:
                if current_phase_reset_idx < len(phase_reset_times[sig_idx])-1:
                    if t >= phase_reset_times[sig_idx][current_phase_reset_idx+1]:
                        current_phase_reset_idx += 1
                    if current_phase_reset_idx < len(phase_reset_times[sig_idx])-1:
                        phase = (t-phase_reset_times[sig_idx][current_phase_reset_idx])/(phase_reset_times[sig_idx][current_phase_reset_idx+1]-phase_reset_times[sig_idx][current_phase_reset_idx])
                    else:
                        # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                        phase = (t-phase_reset_times[sig_idx][-1])/(phase_reset_times[sig_idx][-1]-phase_reset_times[sig_idx][-2])
                else:
                    # calculate phase assuming that the rightmost incomplete cycle has the same frequency as the last complete cycle
                    phase = (t-phase_reset_times[sig_idx][-1])/(phase_reset_times[sig_idx][-1]-phase_reset_times[sig_idx][-2])
                    
            bin_idx = int((phase%1)*nbins)
            if not np.isnan(sig[sig_idx][t_idx]):
                phase_bins[bin_idx].append(sig[sig_idx][t_idx])
                
    result = []
    for bin_idx in range(nbins):
        bs_result = bs.bootstrap(np.array(phase_bins[bin_idx]), stat_func = bs_stats.mean, num_iterations = bsiter, num_threads = n_threads, alpha = alpha)
        result.append((bs_result.value, bs_result.lower_bound, bs_result.upper_bound))
        
    return result

def get_oscillation_extrema(sig, fs, passband, bp_order = 2, extrema_order = 1, method = "trough"):
    """
    Obtains theta-oscillation band peak or trough times.
    
    Parameters
    ----------
    sig : numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    passband : tuple
        Lower and upper theta-band frequency limits (f_lower, f_upper).

    bp_order : int
        Butterworth filter order. Since zero-phase filtering is used, total order will be double.
        
    extrema_order : int
        How many points on each side to use for the comparison to consider. 
        
    Returns
    -------
    np.array
        Peak or through times in [s].
    """
    if method == "trough":
        return scipysig.argrelextrema(filt_butter_bp(sig = sig, fs = fs, lowcut = passband[0], highcut = passband[1], order = bp_order), np.less, order = extrema_order)[0]/fs
    elif method == "peak":
        return scipysig.argrelextrema(filt_butter_bp(sig = sig, fs = fs, lowcut = passband[0], highcut = passband[1], order = bp_order), np.greater, order = extrema_order)[0]/fs
    
# note: DOES NOT SEEM TO WORK WELL
def to_phase_domain(sig, fs, phase_reset_times, nbins = 20):
    """
    Converts a waveform from time-domain to phase domain.

    Parameters
    ----------
    sig : numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    phase_reset_times : np.array
        Time in [s] at which phase is set to 0 deg.

    nbins : int
        Number of phase bins dividing a cycle.

    Returns
    -------
    np.array
        Signal is mapped to a regular phase grid

    """
    # convert from time-samples to phase-samples w.r.t phase reset times
    phase_samples = np.interp(np.arange(len(sig))/fs, phase_reset_times, np.arange(len(phase_reset_times))*360, left = np.nan, right = np.nan)
    # obtain phase-regularized signal amplitude
    crop_left = 0
    crop_right = len(phase_samples)
    for i in range(len(phase_samples)):
        if np.isnan(phase_samples[i]):
            crop_left += 1
        else:
            break
    for i in reversed(range(len(phase_samples))):
        if np.isnan(phase_samples[i]):
            crop_right -= 1
        else:
            break
    out = np.interp(np.arange(0, len(phase_reset_times)*360, 360/nbins), phase_samples[crop_left:crop_right], sig[crop_left:crop_right], left = np.nan, right = np.nan)
    # remove nan
    return out[~np.isnan(out)] 

# 04-10-19: AN: this method does not work too well, a great proportion of throughs are detected but there is a non-negligible number that is missed.
def get_wave_throughs_old(sig, fs, est_freq_band = (5, 10), wt_freq = (1, 50, 0.5), wt_ncycles = 5):
    """
    Extracts times at which a given waveform is at a through.
    
    Parameters
    ----------
    sig : numpy.array
        Waveform to analyze.
        
    fs : float
        Sampling frequency in [Hz].
        
    est_freq_band : tuple
        Lower and upper frequency limits for through estimation as (f_lower, f_upper).
        
    wt_freq : tuple
        Wavelet transform range and resolution (bin size) as (f_lower, f_upper, f_bin) in [Hz].
        
    wt_ncycles : int
        Number of cycles in a wavelet at all frequencies. This is the tradeoff between temporal and frequency estimation.
        
    Returns
    -------
    through_times : np.array
        Waveform through times in [s].
    """
    # wavelet frequency bins
    wt_freqs = np.arange(*wt_freq)
    # generate wavelets
    w = mne.time_frequency.tfr.morlet(sfreq = fs, freqs = wt_freqs, n_cycles = wt_ncycles, sigma = None, zero_mean = False)
    # perform continuous wavelet transform with set of wavelets
    wt = np.reshape(mne.time_frequency.tfr.cwt(np.reshape(sig,(1, -1)), w, use_fft = True, mode = 'same', decim = 1), (len(wt_freqs), -1))
    # calculate power spectral density
    wt_psd = np.abs(wt)**2
    # get frequency bin index with highest power spectral density
    max_psd_freq_idx = wt_psd[int((est_freq_band[0]-wt_freq[0])/wt_freq[2]):int((est_freq_band[1]-wt_freq[0])/wt_freq[2])+1,:].argmax(axis=0) + \
                       int(est_freq_band[0]/wt_freq[2])
    # get frequency in the range of est_freq_band with highest power spectral density
    max_psd_freq = max_psd_freq_idx*wt_freq[2]
    # get phase of frequency component with highest power spectral density in the range est_freq_band
    phase_max_psd = np.empty((np.shape(wt)[1],), dtype = complex)
    for idx in range(0, np.shape(wt)[1]):
        phase_max_psd[idx] = np.angle(wt[max_psd_freq_idx[idx], idx])          
    # identify indices where phase wraps around and sign changes
    throughs_idx = np.zeros((np.shape(wt)[1],),dtype=int)
    through_times = []
    for idx in range(0, np.shape(wt)[1]-1):
        if phase_max_psd[idx] > phase_max_psd[idx+1] and phase_max_psd[idx]*phase_max_psd[idx+1] < 0:
            throughs_idx[idx] = 1
            through_times.append(idx/fs)
            
    return np.array(through_times)
    
def plot_mfreqz(b, a, fs):
    """
    Plot FIR filter magnitude and phase.
    
    Parameters
    ----------
    b, a : array_like
        Filter coefficients.
        
    fs : float
        Sampling frequency in [Hz].
    """
    w,h = scipysig.freqz(b,a)
    h_dB = 20 * np.log10 (abs(h))
    plt.figure(figsize = (15,5))
    plt.subplot(211)
    plt.plot(w/max(w)*fs/2,h_dB)
    plt.ylim(-150, 5)
    plt.ylabel('Magnitude (dB)')
    plt.xlabel(r'Frequency (Hz)')
    plt.title(r'Frequency response')
    plt.subplot(212)
    h_Phase = np.unwrap(np.arctan2(np.imag(h), np.real(h)))
    plt.plot(w/max(w)*fs/2, h_Phase)
    plt.ylabel('Phase (radians)')
    plt.xlabel(r'Frequency (Hz)')
    plt.title(r'Phase response')
    plt.subplots_adjust(hspace=0.5)

def plot_impz(b, a, fs):
    """
    Plot filter step and impulse response.
    """
    l = len(b)
    impulse = np.repeat(0.,l); impulse[0] =1.
    x = np.arange(0,l)
    response = scipysig.lfilter(b,a,impulse)
    plt.figure(figsize = (15,5))
    plt.subplot(211)
    plt.stem(x, response)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(r'Impulse response')
    plt.subplot(212)
    step = np.cumsum(response)
    plt.stem(x, step)
    plt.ylabel('Amplitude')
    plt.xlabel(r'n (samples)')
    plt.title(r'Step response')
    plt.subplots_adjust(hspace=0.5)
    
def downsample(sig, k):
    """
    Downsamples a 1D numpy.ndarray a factor of k.
    
    Parameters
    ----------
    sig : 1D numpy.ndarray
        Signal.
        
    k : int
        Downsampling factor.
        
    Returns
    -------
    1D numpy.ndarray
        Downsampled signal.
    """
    pad_size = int(np.ceil(float(sig.size)/k)*k) - sig.size
    sig_padded = np.append(sig, np.zeros(pad_size)*np.NaN)
    return bn.nanmean(sig_padded.reshape(-1,k), axis = 1)
    
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
    
def nan_helper(y):
    """Helper to handle indices and logical indices of NaNs.

    Input:
        - y, 1d numpy array with possible NaNs
    Output:
        - nans, logical indices of NaNs
        - index, a function, with signature indices= index(logical_indices),
          to convert logical indices of NaNs to 'equivalent' indices
    Example:
        >>> # linear interpolation of NaNs
        >>> nans, x= nan_helper(y)
        >>> y[nans]= np.interp(x(nans), x(~nans), y[~nans])
    """

    return np.isnan(y), lambda z: z.nonzero()[0]
    
def _double_exp_bleach_fit(sig, fs, tau_half_min = 0.1):
    """
    Fits a double exponential bleaching model to fluorescence data. Transients will be excluded from the model, i.e. the bleaching fit
    is done to the decaying baseline level.
    
    Parameters
    ----------
    sig : 1D numpy.nd array
        Fluorescence data. 
        
    fs : float
        Sampling rate in [Hz].
        
    tau_half_min : float
        Shortest fluorescence half-life bleaching time constant expected in [s].

    Returns
    -------
    tuple of (1D numpy array, dict)
        Returns fitted vector and fit parameters as dict with keys:
            'amp_fast': fast bleaching component amplitude.
            'tau_half_fast': fast component half-life in [s].
            'amp_slow': slow bleaching component amplitude.
            'tau_half_slow': slow component half-life in [s].
            'steady': state state component amplitude.
    """
    double_exp_bleach_fit_fn = lambda x,a_fast,k_fast,a_slow,k_slow,steady: a_fast*np.exp(x*k_fast)+\
                                   a_slow*np.exp(x*k_slow)+steady
    
    nan_helper = lambda x: (np.isnan(x), lambda z: z.nonzero()[0])

    # exclude NaN values in signal to be able to calculate a fit
    nans, x = nan_helper(sig)
    nan_idx = x(nans)
    sig_nonans = np.delete(sig, nan_idx)
    t = np.delete(np.arange(len(sig)).astype('float'), nan_idx)/fs
    
    k_initial = np.log(sig_nonans[-1]/sig_nonans[0])/t[-1]
    # enforce k_initial to be within bounds
    if not np.isfinite(k_initial) or k_initial < np.log(0.5)/tau_half_min:
        k_initial = np.log(0.5)/tau_half_min
    
    # note: for a_fast and a_slow changed min value from 0 to -max_sig on 9/27/20 to allow for
    # a strange exponential brightening of mRuby3 in certain recordings; don't know why this is happening.
    max_sig = np.max(sig_nonans)
    model_double_exp_bfit = lmfitModel(double_exp_bleach_fit_fn)
    model_double_exp_bfit.set_param_hint('a_fast', value = 0, min = -max_sig, max = max_sig)
    model_double_exp_bfit.set_param_hint('k_fast', value = k_initial, max = 0, min = np.log(0.5)/tau_half_min)
    model_double_exp_bfit.set_param_hint('a_slow', value = sig_nonans[0]-sig_nonans[-1], min = -max_sig, max = max_sig)
    model_double_exp_bfit.set_param_hint('k_slow', value = k_initial, max = 0, min = np.log(0.5)/tau_half_min)
    model_double_exp_bfit.set_param_hint('steady', value = sig_nonans[-1], min = 0, max = max_sig)
   
    result_double_exp_bfit = model_double_exp_bfit.fit(sig_nonans, x = t, method = 'least_squares')
   
    fit_fn = lambda t: result_double_exp_bfit.best_values['a_fast']*np.exp(t*result_double_exp_bfit.best_values['k_fast'])+\
                       result_double_exp_bfit.best_values['a_slow']*np.exp(t*result_double_exp_bfit.best_values['k_slow'])+\
                       result_double_exp_bfit.best_values['steady']
    # assemble fit parameters and convert rate to halving time constant; also ensure tau_fast < tau_slow
    fit_par = {}
    if result_double_exp_bfit.best_values['k_fast'] > result_double_exp_bfit.best_values['k_slow']: # if the slow rate is faster than the fast rate, swap
        fit_par['amp_fast'] = result_double_exp_bfit.best_values['a_slow']
        fit_par['tau_half_fast'] = np.log(0.5)/result_double_exp_bfit.best_values['k_slow']
        fit_par['amp_slow'] = result_double_exp_bfit.best_values['a_fast']
        fit_par['tau_half_slow'] = np.log(0.5)/result_double_exp_bfit.best_values['k_fast']
    else:
        fit_par['amp_fast'] = result_double_exp_bfit.best_values['a_fast']
        fit_par['tau_half_fast'] = np.log(0.5)/result_double_exp_bfit.best_values['k_fast']
        fit_par['amp_slow'] = result_double_exp_bfit.best_values['a_slow']
        fit_par['tau_half_slow'] = np.log(0.5)/result_double_exp_bfit.best_values['k_slow']
    fit_par['steady'] = result_double_exp_bfit.best_values['steady']

    return fit_fn(np.arange(len(sig)).astype('float')/fs), fit_par

def _natural_cubic_spline_bleach_fit(sig, fs, tstart = 0.1, nknots = 10):
    """
    Fits a natural cubic spline to fluorescence data subject to bleaching that is more complex than a bi-exponential model.

    Parameters
    ----------
    sig : 1D numpy.nd array
        Fluorescence data.     
    fs : float
        Sampling rate in [Hz].
    tstart : float
        Start time in [s] to consider for geometrically placing knot points. Cannot be 0.
    nknots : int
        Number of cubic spline knot points to control smoothess. For typical fluorescence bleaching use 5-15 points.
    """
    # exclude NaN values in signal to be able to calculate a fit
    nan_helper = lambda x: (np.isnan(x), lambda z: z.nonzero()[0])
    nans, x = nan_helper(sig)
    nan_idx = x(nans)
    sig_nonans = np.delete(sig, nan_idx)
    time = np.arange(len(sig))/fs
    time_nonans = np.delete(time, nan_idx)
        
    model = _get_natural_cubic_spline_model(time_nonans, sig_nonans, knots = list(np.geomspace(tstart, time_nonans[-1], nknots)))
    return model.predict(time)

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
    nan_itervals(tmp_sig, fs, exclude)
        
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

# deprecated, use debleach_fl instead
def debleach(sig, maxfev = 20000):
    """
    Removes fluorescence bleaching by normalizing it with a double exponential fit function with offset.
    
    Parameters
    ----------
    sig : array_like
        Fluorescence signal
        
    maxfev : int 
        Maximum number of function evaluations to fit a bleaching function.
    """    
    return np.array(sig)/bleach_fit(sig, maxfev)

# deprecated, use bleach_fit_fl
def bleach_fit(sig, maxfev = 50000, b_init = [1e-3, 1e-4, 1e-2, 1e-5, 1e-1, 1e-6, 0]):
    """
    Fits a double exponential bleaching curve to a fluorescence trace.
    Parameters
    ----------
    sig : array_like
        Fluorescence signal
        
    maxfev : int 
        Maximum number of function evaluations to fit a bleaching function.
    """
    # interpolate NaN values in signal to be able to calculate a fit
    """
    nans, x = nan_helper(sig)
    sig_nan_interp = sig.copy()
    sig_nan_interp[nans] = np.interp(x(nans), x(~nans), sig_nan_interp[~nans])
    fit = lambda x,a,b,c,d,e: a*np.exp(-b*x)+c*np.exp(-d*x)+e
    sig_fit_par, _ = curve_fit(fit, xdata = np.arange(0, len(sig)), ydata = sig_nan_interp, maxfev = maxfev, \
                              p0 = (sig_nan_interp[0], b_initial, 0, 0, sig_nan_interp[-1]))
    sig_fit = fit(np.arange(0, len(sig)), *tuple(sig_fit_par))
    """
    nans, x = nan_helper(sig)
    nan_idx = x(nans)
    sig_nonans = np.delete(sig, nan_idx)
    idx = np.delete(np.arange(0, len(sig)), nan_idx)
    fit = lambda x,a,b,c,d,e: a*np.exp(-b*x)+c*np.exp(-d*x)+e
    # try varying levels of fluorescence bleaching, 1e-3 - typical, 1e-4 - slow, 1e-2 fast
    err_ctr = 0
    for bi in b_init:
        try:
            sig_fit_par, _ = curve_fit(fit, xdata = idx, ydata = sig_nonans, maxfev = maxfev, \
                                      p0 = (sig_nonans[0]-sig_nonans[-1], bi, 0, 0, sig_nonans[-1]))
            sig_fit = fit(np.arange(0, len(sig)), *tuple(sig_fit_par))
            break
        except RuntimeError as e:
            print("skipping initial bleaching rate {}.\n".format(bi))
            err_ctr += 1
            # raise runtime error if there are no more initial bleaching rate values to try
            if err_ctr == len(b_init):
                raise e

    return sig_fit

def wrap_dist_on_belt(dist, belt_length):
    """
    Wraps distance traveled array on a circular belt to belt position.
    
    Parameters
    ----------
    dist : array_like
        Traveled distance.
        
    belt_length : float
        Length of belt.
        
    Returns
    -------
    1D numpy.ndarray
    """
    return np.array(dist)%belt_length
    
def cusum_old(x, mean_shift = (('std', 1), ('std', 1)), c = 0, x_mean = None, x_std = None):
    """
    Parameters
    ----------
    x : 1D np.array
        Waveform to use.
    c : float
        Initial cummulative sum parameter, usually set to 0 but if > 0, it can react
        faster at start.
    mean_shift : tuple of 2 tuples
        Specifies cusum shift sensitivity as anticipated shift size either in terms of signal
        standard deviations or as mean level. First tuple specifies upward sensitivity, while second
        downward. If anticipated mean_shift is specified as standard deviations 'std', second tuple element
        in each direction is the number of std away from mean, otherwise if using 'mean', it is the anticipated
        value.
    x_mean : None or float
        User provided mean for x waveform process. If not given, this is estimated from first 25 non-nan samples.
              
    x_std : None or float
        User provided standard deviation for x waveform process. If not given, this is estimated from first 25 non-nan samples.
        
    Output
    ------
    (c_up, c_down)
        Up and down cumulative sums normalized to the x vector standard deviation (either provided or estimated from first 25 samples).    
    """
    assert len(x)
    # estimate mean and std from first 25 non-nan samples if not user provided
    if x_mean is None:
        x_mean = np.mean(x[~np.isnan(x)][:25])
    if x_std is None:
        x_std = np.std(x[~np.isnan(x)][:25])
    c_up = np.empty((len(x),))
    c_down = np.empty((len(x),))
    c_up[0] = c
    c_down[0] = -c
    if mean_shift[0][0] == 'std':
        k_up = 0.5*mean_shift[0][1]*x_std
    elif mean_shift[0][0] == 'mean':
        k_up = 0.5*(mean_shift[0][1]-x_mean)
    else:
        raise ValueError
    if mean_shift[1][0] == 'std':
        k_down = 0.5*mean_shift[1][1]*x_std
    elif mean_shift[1][0] == 'mean':
        k_down = 0.5*(x_mean-mean_shift[1][1])
    else:
        raise ValueError
    
    # upward cusum
    for i in range(1, len(x)):
        c_up[i] = max(0, c_up[i-1]+x[i]-x_mean-k_up)
    # downward cusum
    for i in range(1, len(x)):
        c_down[i] = min(0, c_down[i-1]+x[i]-x_mean+k_down)
        
    return c_up/x_std, c_down/x_std

def cusum(x, target, shift_factor, c = 0, direction = 'fwd'):
    """
    Forward or backward CUSUM with target signal.
    
    Parameters
    ----------
    x : 1D np.array
        Waveform to use.
        
    fs : float
        Sampling frequency in [Hz].
        
    target : None, float or 1D numpy array
        User provided target signal for x waveform process. If None, target is the median of x.
            
    shift_factor : tuple of float
        Specifies cusum shift sensitivity as anticipated shift size as a multiplicative factor.

    c : float
        Initial cummulative sum parameter, usually set to 0 but if > 0, it can react faster at start.
        
    direction : str
        CUSUM calculation direction, choose 'fwd' for forward or 'bkwd' for backward in time; default 'fwd'.
                  
    Output
    ------
    1D numpy array
        Up or down cumulative sums normalized to the interquartile range of residuals x-target.
    """
    assert len(x)
    if target is None:
        target = bn.nanmedian(x)
    # calculate interquartile range
    iqr = stats.iqr(x-target, nan_policy = 'omit')
    if shift_factor > 1:
        c_up = np.empty((len(x),))
        if direction == 'fwd':
            c_up[0] = c
        elif direction == 'bkwd':
            c_up[len(x)-1] = c
        else:
            raise ValueError
            
        k_up = 0.5*target*(shift_factor-1)
        # upward cusum
        if isinstance(target, Iterable):
            if direction == 'fwd':  
                for i in range(1, len(x)):
                    c_up[i] = max(0, c_up[i-1]+x[i]-target[i]-k_up[i])
            elif direction == 'bkwd':
                for i in range(len(x)-2, 0, -1):
                    c_up[i] = max(0, c_up[i+1]+x[i]-target[i]-k_up[i])
            else:
                raise ValueError
                                
        else:
            if direction == 'fwd':
                for i in range(1, len(x)):
                    c_up[i] = max(0, c_up[i-1]+x[i]-target-k_up)
            elif direction == 'bkwd':
                for i in range(len(x)-2, 0, -1):
                    c_up[i] = max(0, c_up[i+1]+x[i]-target-k_up)
            else:
                raise ValueError
            
        return c_up/iqr

    elif shift_factor < 1:
        c_down = np.empty((len(x),))
        if direction == 'fwd':
            c_down[0] = -c
        elif direction == 'bkwd':
            c_down[len(x)-1] = -c
        else:
            raise ValueError
        
        k_down = 0.5*target*(1-shift_factor)
        # downward cusum
        if isinstance(target, Iterable):
            if direction == 'fwd':
                for i in range(1, len(x)):
                    c_down[i] = min(0, c_down[i-1]+x[i]-target[i]+k_down[i])
            elif direction == 'bkwd':
                for i in range(len(x)-2, 0, -1):
                    c_down[i] = min(0, c_down[i+1]+x[i]-target[i]+k_down[i])
            else:
                raise ValueError
        else:
            if direction == 'fwd':
                for i in range(1, len(x)):
                    c_down[i] = min(0, c_down[i-1]+x[i]-target+k_down)
            elif direction == 'bkwd':
                for i in range(len(x)-2, 0, -1):
                    c_down[i] = min(0, c_down[i+1]+x[i]-target+k_down)
            else:
                raise ValueError

        return c_down/iqr

    else:
        raise ValueError("cusum undefined for shift factor = 1, i.e. no shift")

def convert_bool_epochs_to_intervals(sig, fs = None):
    """
    Converts a bool iterable into time or array index itervals when a transition from False to True happens.

    Parameters
    ----------
    sig : array_like of bool
        Bool signal with False/True transitions.

    fs : float
        Sampling frequency in [Hz].
    
    Returns
    -------
    list of 2 element tuples
        epoch (<start>, <end>) where <start>, <end> is either array index or time in [s] if fs is given.
    """
    epochs = []
    idx = 0
    while idx < len(sig):
        # scan for epoch start index
        if not sig[idx]:
            idx += 1
            continue
        # count number of samples in the epoch
        n_samp = 0
        while idx+n_samp < len(sig) and sig[idx+n_samp]:
            n_samp += 1
        # convert to interval
        if fs is None:
            epochs.append([idx, idx+n_samp-1])
        else:
            epochs.append([idx/fs, (idx+n_samp-1)/fs])
        # jump over epoch
        idx += n_samp

    return epochs

def nan_itervals(sig, fs, intervals):
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
                
def intervals_to_bool(intervals, fs, nsamp):
    """
    Sets given array slices or time intervals to True in a bool array.
    
    Parameters
    ----------
    intervals : list of 2 element tuples
        (start, end) intervals in [s] to be set to True.

    fs : float or None
        If given, sampling frequency in [Hz].

    nsamp : int
        Number of samples in the output bool vector. Ensure number of samples is larger than given intervals.
        
    Returns
    -------
    1D numpy array of bool 
    """
    out = np.full((nsamp,), False)
    for interval in intervals:
        if interval:
            if fs is not None:
                if interval[1] is not None:
                    out[int(interval[0]*fs):int(interval[1]*fs)] = True
                else:
                    out[int(interval[0]*fs):] = True
            else:
                if interval[1] is not None:
                    out[interval[0]:interval[1]] = True
                else:
                    out[interval[0]:] = True
    return out
    
def remove_PLI(x, fs, M, B, P, W, f_ac = []):
    """
    Remove power line interference.

    This is an implementation of the proposed algorithm in,
    M. R. Keshtkaran and Z. Yang, "A fast, robust algorithm for power line 
    interference cancellation in neural recording," J. Neural Eng., vol. 11,
    no. 2, p. 026017, Apr. 2014.
    http://iopscience.iop.org/1741-2552/11/2/026017
    http://arxiv.org/abs/1402.6862

    Parameters
    ----------
    x : 1D numpy array
        Contaminated signal.

    fs : float
        Sample rate in [Hz].

    M : int
        Fundamental plus number of harmonics to remove, i.e. for 60 Hz and 120 Hz use 2.

    B : iterable of float, 3 elements
        Fundametal frequency notch filter settings. Contains three elements [B0,Binf,Bst]:
            - B0, initial notch bandwidth of the frequency estimator.
            - Binf, asymptotic notch bandwidth of the frequency estimator.
            - Bst, rate of convergence to 95% of the asymptotic bandwidth Binf.

    P : iterable of float, 3 elements
        Frequency estimator settling time. Contains three elements [P0,Pinf,Pst]: 
            - P0, initial settling time of the frequency estimator.
            - Pinf, asymptotic settling time of the frequency estimator.
            - Pst, rate of convergence to 95% of the asymptotic settling time.

    W : iterable of float
        Settling time in [s] of the amplitude and phase estimator. If 1 element, same settling time is applied to
        the fundamental and its harmonics. If >1, must have same number of elements as harmonics.

    f_ac : iterable of float, optional
        The nominal AC frequency if known (50 Hz or 60 Hz) specified as e.g. [60]. If providing 2 elements, these are the minimum
        and maximum fundamental line frequencies expected. If left empty, a range of 40 to 70 Hz is assumed.

    Returns
    -------
    s : 1D numpy array
        Clean signal.

    Example
    -------
    from __future__ import division
    from scipy import signal as scipysig
    import numpy as np
    from matplotlib import pyplot as plt

    fs = 500
    n = 120*fs # 2 min sequence 
    t = 2*np.pi*np.arange(n)/fs
    fline = 60 + np.random.randn() # random interference frequency
    s = scipysig.lfilter([1,1],[1,-0.99],100*np.random.randn(n)) #1/f^2 PSD
    p = 80*np.sin(fline*t+np.random.randn()) + 50*np.sin(2*fline*t+np.random.randn()) + \
        20*np.sin(3*fline*t+np.random.randn()) # interference 
    x = s + p
    sbar = remove_PLI(x, fs, 3, [100,0.01,4], [0.1,2,5], [3, 2.5, 2.5])

    f, P_x = scipysig.welch(x, fs, nperseg = 2049, average = 'median')
    f, P_sbar = scipysig.welch(sbar, fs, nperseg = 2049, average = 'median')
    plt.figure(figsize = (5,5))
    plt.semilogy(f, P_x, 'b-', alpha = 0.5)
    plt.semilogy(f, P_sbar, 'r-', alpha = 0.5)
    plt.xlabel('frequency [Hz]')
    plt.ylabel('PSD [$\mu V^2$/Hz]')

    Licence
    -------
    Downloaded from: https://github.com/mrezak/removePLI
    Author: Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
    Copyright (c) 2013, Mohammad Reza Keshtkaran <keshtkaran.github@gmail.com>
    All rights reserved.
    This program is provided "AS IS" for non-commercial, educational 
    and reseach purpose only. Any commercial use, of any kind, of 
    this program is prohibited. The Copyright notice should remain intact.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

    """
    assert M>0

    # removing the mean
    x -= np.mean(x)
    N = len(x)
    s = np.zeros((N,))

    # 3dB cutoff bandwidth
    alpha_f = (1-np.arctan(np.pi*B[0]/fs))/(1+np.arctan(np.pi*B[0]/fs)) # initial, \alpha_0
    alpha_inf = (1-np.tan(np.pi*B[1]/fs))/(1+np.tan(np.pi*B[1]/fs)) # asymptotic    
    alpha_st = np.exp(np.log(0.05)/(B[2]*fs+1)) # rate of change

    # frequency estimator's forgetting factors
    lambda_f = np.exp(np.log(0.05)/(P[0]*fs+1)) # initial
    lambda_inf = np.exp(np.log(0.05)/(P[1]*fs+1)) # asymptotic
    lambda_st = np.exp(np.log(0.05)/(P[2]*fs+1)) # rate of change
    # smoothing parameter (cut-off freq set at 90 Hz)
    gmma = (1-np.tan(0.5*np.pi*min(90,fs/2)/fs))/(1+np.tan(0.5*np.pi*min(90,fs/2)/fs))

    # phase/amplitude estimator forgetting factor
    lambda_a = np.exp(np.log(0.05)/(np.array(W)*fs+1))
    if len(lambda_a) == 1:
        lambda_a = lambda_a[0]*np.ones((M,))

    # initializing variables
    kappa_f = 0
    kappa_k = np.zeros((M+2,))
    D = 10
    C = 5
    f_n1 = 0
    f_n2 = 0

    # -- Alternative initialization:
    #    kappa_f = np.cos(55*2*np.pi/fs)
    # --

    # initializing the first oscillator
    u_kp = 1*np.ones((M,)) # u_k
    u_k = 1*np.ones((M,)) # u'_k

    # initializing the RLS parameters
    r1 = 10*np.ones((M,))
    r4 = 10*np.ones((M,))
    a = np.zeros((M,))
    b = np.zeros((M,))

    # IIR bandpass filtering:
    if f_ac: # if AC frequency is known
        if len(f_ac) == 2:
            Fc1 = f_ac[0]  # first cutoff frequency
            Fc2 = f_ac[1]  # second cutoff frequency               
        elif len(f_ac) == 1:
            # custom center frequency of pass band
            Fc1 = f_ac[0]-2  # first cutoff frequency
            Fc2 = f_ac[0]+2  # second cutoff frequency
        else:
            raise ValueError()
    else: # if AC frequency is not known
        # default 40--70 Hz pass band
        Fc1 = 40  # first cutoff frequency
        Fc2 = 70  # second cutoff frequency

    ordr = 4 # order
    x_f = scipysig.sosfilt(sos_butter_bp(Fc1, Fc2, fs, ordr), x) # band-pass filter
    x_f = np.insert(np.diff(x_f),0,0) # first difference to remove 1/f^2 PSD trend

    # --------- start of data processing
    for n in range(N):
        # lattice filter
        f_n = x_f[n] + kappa_f*(1+alpha_f)*f_n1 - alpha_f*f_n2
        
        # frequency estimation
        C = lambda_f*C+(1-lambda_f)*f_n1*(f_n+f_n2)
        D = lambda_f*D+(1-lambda_f)*2*f_n1**2
        kappa_t=C/D
        if kappa_t <-1:
            kappa_t = -1
        if kappa_t > 1:
            kappa_t= 1
        kappa_f = gmma*kappa_f + (1-gmma)*kappa_t
        
        f_n2 = f_n1
        f_n1 = f_n # updating lattice states 

        # bandwidth and forgetting factor updates
        alpha_f = alpha_st*alpha_f + (1-alpha_st)*alpha_inf
        lambda_f = lambda_st*lambda_f + (1-lambda_st)*lambda_inf
        
        # Discrete-Time Oscillators
        kappa_k[1] = 1
        kappa_k[0] = kappa_f
        
        e = x[n]
        for k in range(M): # for each harmonic do
            # calculating cos(kw) for k=1,2...
            kappa_k[k+2] = 2*kappa_f*kappa_k[k+1] - kappa_k[k] 
            
            # oscillator
            tmp = kappa_k[k+2]*(u_kp[k]+u_k[k])
            tmp2 = u_kp[k]
            u_kp[k] = tmp - u_k[k]
            u_k[k] = tmp + tmp2
            
            # gain Control
            G = 1.5 - (u_kp[k]**2 - (kappa_k[k+2]-1)/(kappa_k[k+2]+1)*u_k[k]**2)
            if G <= 0:
                G = 1
            u_kp[k] = G * u_kp[k]
            u_k[k] = G * u_k[k]

            # Phase/Amplitude Adaptation
            sincmp = a[k]*u_k[k] + b[k]*u_kp[k]
            e = e - sincmp
            # --- simplified RLS
            r1[k] = lambda_a[k]*r1[k] + u_k[k]**2
            r4[k] = lambda_a[k]*r4[k] + u_kp[k]**2
            a[k] = a[k] + u_k[k]*e/r1[k]
            b[k] = b[k] + u_kp[k]*e/r4[k]
            # ------
       
        s[n] = e

    return s

def crop_window(sig, fs, timepoint, window):
    """
    Crop a portion of a waveform centered on given time point.

    Parameters
    ----------
    sig : 1D numpy array
        Waveform to crop from.
    fs : float
        Sampling rate in [Hz].
    timepoint : float
        Timepoint in [s] around which to center the crop window.
    window : float
        Crop window duration in [s].
    

    Returns
    -------
    1D numpy array
        Cropped waveform. If crop window exceeds waveform extent, output is padded with nans.
    """
    # center crop window on chosen time point
    crop_window_center_idx = int(round(timepoint*fs))
    # desired crop window size
    crop_window_nsamples = util.round_to_nearest_odd_int(window*fs)
    # actual crop window is needed to account for cropping at edges
    actual_crop_window_left_idx = max(0, crop_window_center_idx - int((crop_window_nsamples-1)/2))
    actual_crop_window_right_idx = min(len(sig),crop_window_center_idx + int((crop_window_nsamples-1)/2)+1)
    n_left_nan_pad = abs(min(0, crop_window_center_idx - int((crop_window_nsamples-1)/2)))
    cropped_sig = sig[actual_crop_window_left_idx:actual_crop_window_right_idx]
    crop_buffer = np.full((crop_window_nsamples,), np.nan)
    crop_buffer[n_left_nan_pad:n_left_nan_pad+len(cropped_sig)] = cropped_sig

    return crop_buffer

def crop_around_idx(sig, idx, width):
    """
    Crops a signal around a given index and centers the cropped portion within a given window.
    
    Parameters
    ----------
    sig : 1D numpy.array
        Signal.
    idx : int
        Array index. If index and window combination falls out of array, returns empty numpy.array
    width : int
        Window width in samples. Number must be odd, >=3 such that the given index is centered within the window.
    """
    assert width>=3 and width%2
    w_halfwidth = int((width-1)/2)
    left_idx = idx - w_halfwidth
    right_idx = idx + w_halfwidth+1
    if left_idx<0 or right_idx>len(sig):
        return np.array([])
    else:
        return sig[left_idx:right_idx]
        
def crop_events(sig1, sig2, fs1, fs2, threshold, window, crop = (0, None), average = False):
    """
    Crop and average event waveforms triggered on event peaks that are larger than threshold.

    Parameters
    ----------
    sig1 : 1D numpy array
        Signal to use for triggering on positive going events.
    sig2 : 1D numpy array
        Signal to use for cropping.
    fs1, fs2 : float
        Sampling rates in [Hz] for signals 1 & 2.
    threshold : float
        Event threshold.
    window : float
        Crop window duration in [s].
    crop : iterable of 2 or iterable of iterable of 2
        Start and end signal time to crop in [s]. Can be same for all signals or specified for each signal.
    average : bool
        If true, return the average of cropped waveforms.

    Returns
    -------
    dict with keys:
        'times': list of float
            Event times in [s] with t = 0 measured from the beginning of sig1.
        'waveforms': 2D numpy.ndarray
            Cropped event waveforms, with row-idx indexing event number.
        'avg': 1D numpy.ndarray
            Average of all event waveforms
    """
    out = {}
    out['times'] = []
    cropped_evts_sig1 = []
    cropped_evts_sig2 = []
    # force 2D input, with first dimension indexing signal
    sig1 = np.atleast_2d(np.array(sig1))
    sig2 = np.atleast_2d(np.array(sig2))
    crop = np.atleast_2d(np.array(crop))
    # ensure same number of event and waveform crop signals
    nsignals = sig1.shape[0]
    assert sig1.shape[0] == sig2.shape[0]
    if crop.shape[0] == 1 and nsignals > 1:
        crop = np.tile(crop, (nsignals, crop.shape[1]))

    for sig_idx in range(nsignals):
        cropped_sig1, (sig1_tstart, sig1_tend) = tslice(sig = sig1[sig_idx], fs = fs1, t_slice = tuple(crop[sig_idx]), return_intervals = True)
        cropped_sig2, (sig2_tstart, sig2_tend) = tslice(sig = sig2[sig_idx], fs = fs2, t_slice = tuple(crop[sig_idx]), return_intervals = True)
        for pk_idx in scipysig.find_peaks(cropped_sig1, height = threshold)[0]:
            pk_time = pk_idx/fs1
            out['times'].append(pk_time+sig1_tstart)
            cropped_evts_sig1.append(crop_window(sig = cropped_sig1, fs = fs1, timepoint = pk_time, window = window))
            cropped_evts_sig2.append(crop_window(sig = cropped_sig2, fs = fs2, timepoint = pk_time, window = window))
    
    out['sig1'] = np.array(cropped_evts_sig1)
    out['sig2'] = np.array(cropped_evts_sig2)
    if average:
        out['avg1'] = np.nanmean(out['sig1'], axis = 0)
        out['avg2'] = np.nanmean(out['sig2'], axis = 0)

    return out

def crop_timepoints(tp, crop):
    """
    Restricts a series of time points within a given crop interval [tstart, tend) where tend may be = None.

    Parameters
    ----------
    tp: iterable of float
        Time points.
    crop : iterable of 2 elements
        Time interval.

    Returns
    -------
    list
    """
    return [tp[idx] for idx in range(len(tp)) if (crop[0]<=tp[idx]<crop[1] if crop[1] is not None else crop[0]<=tp[idx])]

def interp_resample(sig, fs_in, fs_out):
    """
    Resample signal with linear interpolation.

    Parameters
    ----------
    sig : 1D numpy.array
        Input signal.
    fs_in : float
        Input sampling rate in [Hz].
    fs_out : float
        Output sampling rate in [Hz].

    Returns
    -------
    (numpy.array, numpy.array)
        Tuple of resampled time points in [s] and signal.
    """
    if fs_out != fs_in:
        t_in = 1/fs_in*np.arange(len(sig))
        f = scipy_interp1d(t_in, sig, assume_sorted = True)
        t_out = np.linspace(0, t_in[-1], int(t_in[-1]*fs_out)+1)
        return t_out, f(t_out)
    else:
        return 1/fs_in*np.arange(len(sig)), sig

def centroid_time(sig, fs, baseline = 1):
    """
    Obtains the centroid timepoint of a waveform.

    Parameters
    ----------
    sig : numpy 2D nd.array
        Pulse waveforms to average, stacked along the 0-axis.
    fs : float
        Sampling rate in [Hz].
    baseline : float
        Pulse baseline.
    """
    # avoid issues when passing a 1D array
    if len(sig.shape) == 1:
        sig = np.reshape(sig, (1,len(sig)))
    return 1./fs*np.average(np.tile(np.arange(0, sig.shape[1]), (sig.shape[0], 1)),
                          weights = np.abs(sig-baseline), axis = 1)

def pulse_avg(sig, fs, method = 'peak', baseline = 1):
    """
    Averages pulses after temporal alignment.

    Parameters
    ----------
    sig : numpy 2D nd.array
        Pulse waveforms to average, stacked along the 0-axis.
    fs : float
        Sampling rate of input signal in [Hz].
    fs_out : float
        Sampling rate of output average in [Hz].
    method : str
        Alignment method, choose 'centroid', 'peak' or 'through'.
    """    
    # get pulse alignment timepoints
    if method == 'centroid':
        pt = centroid_time(sig, fs = fs, baseline = baseline)
    elif method == 'peak':
        pt = np.argmax(sig, axis = 1)/fs
    elif method == 'through':
        pt = np.argmin(sig, axis = 1)/fs
    else:
        raise ValueError()
    idx_offsets = ((pt-np.mean(pt))*fs).astype(int)
    # apply temporal offsets and average
    return np.nanmean(util.shift_rows(sig, idx_offsets), axis = 0)

def pulse_fwhm(sig, fs, baseline = 0):
    """
    Calculates the full width at half maximum of a pulse.
    Note: peak estimation is done by the pulse maximum or minimum, i.e.
    the pulse needs to have a distinguishable peak or through.
    
    Parameters
    ----------
    sig : 1D numpy.array
        Pulse waveform.
    fs : float
        Sampling rate in [Hz].
    baseline : float
        Pulse baseline.

    Returns
    -------
    float
        Pulse FWHM in [s].
    """
    # adjust signal to work with both positive and negative going pulses
    dt = 1/fs
    sig = np.absolute(sig-baseline)
    max_idx = np.argmax(sig)
    # half maximum
    hm = sig[max_idx]/2
    # start idx
    idx = max_idx
    while idx >= 0:
        if sig[idx] < hm:
            break
        idx -= 1
    start_idx = idx
    # left idx
    idx = max_idx
    while idx < len(sig):
        if sig[idx] < hm:
            break
        idx += 1
    end_idx = idx
    
    start_t = np.interp(x = hm, xp = [sig[start_idx], sig[start_idx+1]], fp = [start_idx*dt, (start_idx+1)*dt], left = np.nan, right = np.nan)
    end_t = np.interp(x = hm, xp = [sig[end_idx], sig[end_idx-1]], fp = [end_idx*dt, (end_idx-1)*dt], left = np.nan, right = np.nan)
    return end_t-start_t

def isolated_evt(evts, spacing):
    """
    Determines if an event index in an array is isolated from other events.

    Parameters
    ----------
    evts : iterable of int
        Array with monotonically increasing event indices.
    spacing : int
        Minimum spacing in units of array index between events to consider events isolated.

    Returns
    -------
    np.array of bool type
        True if event is isolated, False otherwise.
    """
    is_isolated = np.full((len(evts),), True)

    # first and last event
    if len(evts)>1:
        if evts[1]-evts[0]<spacing:
            is_isolated[0] = False
        else:
            is_isolated[0] = True
        if evts[-1]-evts[-2]<spacing:
            is_isolated[-1] = False
        else:
            is_isolated[-1] = True

    # in-between events
    for idx in range(1,len(evts)-1):
        if evts[idx]-evts[idx-1]<spacing or evts[idx+1]-evts[idx]<spacing:
            is_isolated[idx] = False
        else:
            is_isolated[idx] = True

    return is_isolated

def padded_shift(arr, n, fill_val = None):
    """
    Shift and pad an array.

    Parameters
    ----------
    arr : 1D numpy.array
        Input array.
    n : int
        Number of samples to shift. A positive value shifts the array to the right.
    fill_val : None or array element
        Value to use for padding missing elements. If None, the start and end of the array is padded
        with the first or last element value.

    Returns
    -------
    1D numpy.array
    """
    out = np.empty_like(arr)
    if n > 0:
        if fill_val is None:
            out[:n] = arr[0]
        else:
            out[:n] = fill_val
        out[n:] = arr[:-n]
    elif n < 0:
        if fill_val is None:
            out[n:] = arr[-1]
        else:
            out[n:] = fill_val
        out[:n] = arr[-n:]
    else:
        out[:] = arr
    return out

