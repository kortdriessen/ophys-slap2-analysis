from __future__ import division
from lmfit import Model as lmfitModel, conf_interval as lmfit_conf_int, fit_report as lmfit_report #	Non-Linear Least-Squares Minimization and Curve-Fitting
from lmfit import Minimizer, Parameters
from future.utils import iteritems # python 2&3 compatibility to retrieve (key,val) dict pairs
import numpy as np, math, pandas as pd, numpy.ma as ma, copy
import process as proc
from collections import OrderedDict, Iterable, namedtuple
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
from scipy import signal as scipysig, stats as scipy_stats
from numpy.fft import fft, ifft, fft2, ifft2, fftshift
from analysis_AN import util, process as proc
from scipy import optimize as scipyopt
from functools import partial
from multiprocessing import cpu_count
from multiprocess import Pool
import scipy.sparse as scipy_sparse
from seaborn.algorithms import bootstrap as sb_bootstrap
import warnings
import bottleneck as bn # speeds up certain operations on numpy arrays with nan

# multivariate regression
import statsmodels.api as sm

def circmean_with_ci(phi, alpha = 0.05, unit = 'deg', wrap_phase = False):
    """
    Calculates circular variable mean angle and calculates a confidence interval assuming no particular distribution.

    Parameters
    ----------
    phi : 1D numpy nd.array
        Circular variable, in units of [deg] or [rad] based on value of 'unit'.
    unit : str
        Angle unit. Choose 'deg' or 'rad'.
    alpha : float
        Significance level.
    wrap_phase : bool
        If True, wrap phase to (-180, 180] interval instead of (0,360].

    Returns
    -------
    (mean, low_ci, high_ci) : tuple
        Circular angle mean and 1-alpha confidence interval
            mean : float
                Mean.
            low_ci, high_ci : float
                Lower and upper confidence interval at alpha signicance level.
    If mean direction is ill defined, low_ci and high_ci will be numpy.nan
    """
    # convert to radians
    if unit == 'deg':
        phi *= np.pi/180.
    phi_mean = scipy_stats.circmean(phi) # mean in radians

    # ci calculation from Upton & Fingleton (1989) page 220
    # ------------------------------------------------------
    # resultant length
    r = (np.sum(np.cos(phi))**2+np.sum(np.sin(phi))**2)**0.5
    h = 1./len(phi)*(math.cos(2*phi_mean)*np.sum(np.cos(2*phi))+math.sin(2*phi_mean)*np.sum(np.sin(2*phi)))
    # circular variance
    w = len(phi)*(1-h)/(4.*r**2)
    u_alpha = scipy_stats.norm.ppf(1-alpha/2.) 
    
    try:
        ci_delta = math.asin(u_alpha*(2*w)**0.5) # if argument of arcsin > 1, CI is ill defined and returns np.nan
    except ValueError:
        ci_delta = np.nan

    # convert mean and ci's back to deg if needed
    if unit == 'deg':
        phi_mean *= 180/np.pi
        ci_delta *= 180/np.pi

    if wrap_phase:
        phi_mean = phase_wrap(phi_mean)[0]

    return phi_mean, phi_mean-ci_delta, phi_mean+ci_delta

def phase_wrap(phi):
    """
    Wraps phase to (-180, 180] deg interval.
    Parameters
    ----------
    phi : float
        phase in [deg].

    Returns
    -------
    (phase, wrapped) : tuple
        phase : float
        wrapped : bool
            True if wrapping was applied
    """
    if phi == -180:
        wrapped_phi = 180
        wrapped = True
    elif phi < -180:
        wrapped_phi = phi%180
        wrapped = True
    elif phi > 180:
        wrapped_phi = phi%-180
        wrapped = True
    else:
        wrapped_phi = phi
        wrapped = False
    return (wrapped_phi, wrapped)

def fit_lin_circ_data(x, phi, max_slope, phi_unit = 'deg', n_grid = 5000):
    """
    Fit a linear-circular model, method adapted from these references:
     * Kempter et. al., The Journal of Neuroscience, October 21, 2009, 29(42):13232-13241
     * Kempter et. al., Journal of Neuroscience Methods 207 (2012) 113-124 

    The function fits the following linear-circular model:
        Phi = s*X + phi0

    to the provided linear independent data vector x and circular dependent data vector phi.
    For X = 0, Phi = phi0.

    Parameters
    ----------
    x : 1D numpy.ndarray
        Linear data.

    phi : 1D numpy.ndarray
        Circular data, in [deg] or [rad] depending on phi_unit.

    max_slope : float
        Maximum positive or negative slope allowed to model the linear-circular relationship.

    phi_unit : str
        Circular variable phi measurement unit; choose between 'deg' and 'rad'.

    n_grid : int
        Number of grid points to use for brute force optimization and divide the range (-max_slope, max_slope)

    x must have the same number of elements as phi.

    Returns
    -------
    If successful, dict with keys:
        'slope' : float
            Linear-circular regression slope in [deg]|[rad]/[<unit of x>]
        'phi0' : float
            Linear-circular regression phase offsetin [deg]|[rad]
        'corr' : float
            Correlation coefficient in the interval [-1, 1]. Note that for n <= 10 samples,
            the correlation coefficient may be slightly outside [-1, 1]
        'pval' : float
            Significance level

    If fit fails, returns empty dict.
    """
    assert len(x) == len(phi)
    x = x.astype(float)
    # convert to radians
    if phi_unit == 'rad':
        phi = phi.astype(float)
    elif phi_unit == 'deg':
        phi = phi.astype(float)*np.pi/180.
        max_slope *= np.pi/180.
    else:
        raise ValueError("phi_unit can be either 'deg' or 'rad'.")
    
    def min_fn(params, x, phi):
        """
        Objective function to minimize that maximizes the mean resultant length.

        Parameters
        ----------
        s : float
            Linear fit slope measured in [rad]/[unit of x].

        x : 1D numpy.ndarray
            Linear data.

        phi : 1D numpy.ndarray
            Circular data, phase in radians.
        """
        s = params['s']
        # mean resultant vector
        r = ((np.sum(np.cos(phi-s*x))/len(x))**2+(np.sum(np.sin(phi-s*x)/len(x))**2))**0.5
        return -r

    # this thing is not optimizing reliably, i.e. the result is quite variable between different runs
    #opt_result = scipyopt.minimize(partial(min_fn, x = x , phi = phi), [0], args=(), method = 'L-BFGS-B', bounds = [(-abs(max_slope), abs(max_slope))], options = {'iprint':1})
    # brute_result = scipyopt.brute(func = min_fn, ranges = (-abs(max_slope), abs(max_slope)), args=(x,phi), Ns = 5000) 
    # Opt_result = namedtuple('opt_result', ['success', 'message', 'x'])
    # opt_result = Opt_result(True, 'success', brute_result[0])
    params = Parameters()
    params.add('s', value = 0, min = -abs(max_slope), max = abs(max_slope))

    fitter = Minimizer(min_fn, params, fcn_args = (x, phi))
    brute_result = fitter.minimize(method = 'brute', Ns = n_grid)

    #print(lmfit_report(brute_result))

    out = {}
    # slope and offset in radians
    slope = brute_result.params['s']
    phi0 = math.atan2(np.sum(np.sin(phi-slope*x)), np.sum(np.cos(phi-slope*x)))
    # calculate a linear-circular correlation coefficient
    # ---------------------------------------------------
    theta = np.mod(abs(slope)*x, 2*np.pi)
    phi_mean = math.atan2(np.sum(np.sin(phi)), np.sum(np.cos(phi)))
    theta_mean = math.atan2(np.sum(np.sin(theta)), np.sum(np.cos(theta)))
    out['corr'] = (np.sum(np.sin(phi-phi_mean)*np.sin(theta-theta_mean)))/(np.sum(np.sin(phi-phi_mean)**2)*np.sum(np.sin(theta-theta_mean)**2))**0.5

    # calculate p-value
    _lambda = lambda i,j: 1./len(x)*np.sum(np.sin(phi-phi_mean)**i * np.sin(theta-theta_mean)**j)
    z = out['corr']*(len(x)*_lambda(2,0)*_lambda(0,2)/_lambda(2,2))**0.5
    out['pval'] = 1-math.erf(abs(z)/2**0.5)

    # convert to desired unit
    if phi_unit == 'deg':
        out['slope'] = slope*180/np.pi
        out['phi0'] =  phi0*180/np.pi
    else:
        out['slope'] = slope
        out['phi0'] =  phi0
    return out
    
def sinfit(x, offset, amp, phi):
    """
    For fitting a sinusoidal model. 
    90 deg offset is needed because this is relative phase of fluorescence w.r.t LFP,
    a positive value means fluorescence is retarded w.r.t. LFP
    
    Parameters
    ----------
    x : 1D np.ndarray
        Phase in [deg].
    offset, amp, phi : float
        Offset in [%], pk-pk oscillation amplitude in [%] and phase in [deg].
    """
    return offset/100.+0.5*amp/100.*np.sin((x+phi-90)*np.pi/180.)

def sinfit2(x, offset, amp_sin, amp_cos):
    """
    Another approach to fit a sinusoidal model that should be more suitable for searching optimal
    parameters compared to the non-convex space of parameters (amp, phi) amp*sin(x+phi).
    """
    return offset/100.+0.5*amp_sin/100.*np.sin((x-90)*np.pi/180.)+0.5*amp_cos/100.*np.cos((x-90)*np.pi/180.)
    
def constfit(x, m):
    """
    Fits data with a constant value (if normally distributed, minimizing the least squares error would give the mean).
    
    Parameters
    ----------
    x : 1D np.ndarray
    
    m : float
      Constant value.
    """
    return m

def meas_LFP_fluorescence_phase(sig, phi_init = 0, mvmt = '', mint = 10):
    """
    Measures theta-band fluorescence phase w.r.t. local field potential (LFP).
    
    Parameters
    ----------
    sig: dict of dict
        Measurements to combine, with dict labels indicating measurement names and values being dict with labels:
            'fluo': 1D numpy.ndarray
                Fluorescence trace.
            'lfp': 1D numpy.ndarray
                Local field potential.
            'fs': float
                Sampling frequency in [Hz].
            'loc': optional, 1D numpy.ndarray
                Locomotion signal, i.e. sustained movement, during which phase is calculated.
            'still': optional, 1D numpy.ndarray
                Stillness signal, i.e. no movement, during which phase is calculated.
    
    phi_init : float
        Initial value for phase estimation.

    mvmt : str
        Mouse movement type for which phase is calculated. If '' (default), movement type is ignored and phase is calculated for the entire period.
        If 'loc', phase is calculated during locomotion periods only, and if 'still', phase is calculated during still periods only.

    mint : float
        Minimum total valid sample time used for phase estimation in [s]. If total time is less than this value, function returns None.

    Returns
    -------
    dict qith phase result keys or empty dict if there is not enough data to extract phase. When phase estimation is available, dict has:
        Analysis result with keys:
        'amp' : float
            pk-pk oscillation amplitude in [%].
        'amp_se' : float
            Oscillation amplitude standard error.
        'amp_ci' : dict
            If can be calculated, oscillation amplitude confidence intervals, otherwise this is not included.
        'phase' : float
            Oscillation phase relative to local field potential (LFP) in [deg]. Zero phase means the fluorescence oscillation through
            is aligned with the LFP through. Phase wraps around at +-180 deg.
        'phase_wrapped' : bool
            True if phase was wrapped to +- 180 deg, False otherwise. Confidence intervals are centered around the wrapped phase, but
            they are noy wrapped and may exceed +- 180 deg.
        'phase_se' : float
            Oscillation phase standard error.
        'phase_ci' : dict
            If can be calculated, phase confidence intervals, otherwise this is not included.
        'sin_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a sinusoidal model explaining the data w.r.t. a constant model.
        'const_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a constant model explaining the data w.r.t. a sinusoidal model.
        'total_valid_time' : float
            Total duration in [s] of valid samples used for phase estimation. This includes non-nan samples and movement-related gating.
        
    """
    # convert to ordered dict so as to not loose LFP and fluorescence pairing
    sig = OrderedDict(sig)
    # obtain theta oscillation through times
    theta_peaks = [proc.get_oscillation_extrema(meas['lfp'], meas['fs'], passband = (5,10), bp_order = 2, method = "peak") for (key, meas) in sig.iteritems()]
    fs = [meas['fs'] for (key, meas) in sig.iteritems()]
    fluo = [meas['fluo'] for (key, meas) in sig.iteritems()]

    # build locomotion or stillness gating signal
    if mvmt == 'loc':
        gating = [meas['loc'] for (key, meas) in sig.iteritems()] 
    elif mvmt == 'still':
        gating = [meas['still'] for (key, meas) in sig.iteritems()] 
    elif mvmt == '':
        gating = []
    else:
        raise ValueError("mvmt parameter should be 'loc', 'still' or ''.")    

    # assign LFP phase to fluorescence samples and exclude nan values
    try:    
        fl_phase, fl_int, total_valid_time = proc.assign_phase(fluo, fs, theta_peaks, gating)
    except:
        return {}

    if total_valid_time < mint:
        return {}
    # estimate fluorescence oscillation phase using a nonlinear least squares fit to a sinusoidal model
    model_sin = lmfitModel(sinfit)
    model_sin.set_param_hint('offset', value = 100)
    model_sin.set_param_hint('amp', value = 1, min = 0) # >0 min bound needed to avoid 180 deg phase flip if amplitude becomes negative 
    model_sin.set_param_hint('phi', value = phi_init, min = -360, max = 360)
    try:
        result_sin = model_sin.fit(fl_int, x = fl_phase)
    except ValueError:
        util.clrd_print('Cannot calculate phase for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}
    # calculate confidence intervals
    try:
        ci = result_sin.conf_interval(maxiter = 10000)
    except:
        util.clrd_print('Cannot determine confidence intervals.\n', 'warn')
        ci = None 
    # also do a nonlinear least squares fit to a constant model to compare Akaike information criteria and justify sinusoidal model selection over constant value
    model_const = lmfitModel(constfit)
    try:
        result_const = model_const.fit(fl_int, x = fl_phase, m = 1)
    except ValueError:
        util.clrd_print('Cannot calculate constant model for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}

    # output result
    out = {}
    # amplitude best fit
    out['amp'] = result_sin.best_values['amp']
    # amplitude standard error
    out['amp_se'] = result_sin.params['amp'].stderr
    # amplitude confidence intervals as 0, 1-sigma, 2-sigma, 3-sigma tuples of (probability, min, max)
    if ci is not None:
        out['amp_ci'] = {i: ci['amp'][int((len(ci['amp'])-1)/2)-i]+(ci['amp'][i+int((len(ci['amp'])-1)/2)][1],) for i in range(int((len(ci['amp'])-1)/2)+1)}
   
    # phase best fit w.r.t. LFP
    wrapped_phase, wrapped = phase_wrap(result_sin.best_values['phi'])
    out['phase'] = wrapped_phase
    out['phase_wrapped'] = wrapped
    # standard error of phase estimate
    out['phase_se'] = result_sin.params['phi'].stderr
    # phase confidence intervals as 0, 1-sigma, 2-sigma, 3-sigma tuples of (probability, min, max)
  
    if ci is not None:
        out['phase_ci'] = {i: (ci['phi'][int((len(ci['phi'])-1)/2)-i][0], wrapped_phase-(result_sin.best_values['phi']-ci['phi'][int((len(ci['phi'])-1)/2)-i][1]),
                               ci['phi'][i+int((len(ci['phi'])-1)/2)][1]-result_sin.best_values['phi']+wrapped_phase) for i in range(int((len(ci['phi'])-1)/2)+1)}
    
    # using Akaike's information criterion, determine the chance that the provided data can be better explained by a constant model compared to a sinusoidal model
    min_aic = min(result_sin.aic, result_const.aic)
    out['sin_model_rel_likelihood'] = np.exp((min_aic-result_sin.aic)/2.)
    out['const_model_rel_likelihood'] = np.exp((min_aic-result_const.aic)/2.)
    out['total_valid_time'] = total_valid_time
    
    return out 
  
def meas_LFP_theta_band_fluorescence_phase(sig, phi_init = 0, mvmt = '', mint = 10, max_amp = 20):
    """
    Measures theta-band fluorescence phase w.r.t. local field potential (LFP). Second implementation using sinfit2 that is a convex optimization problem.

    Parameters
    ----------
    sig: dict of dict
        Measurements to combine, with dict labels indicating measurement names and values being dict with labels:
            'fluo': 1D numpy.ndarray
                Fluorescence trace.
            'lfp': 1D numpy.ndarray
                Local field potential.
            'lfp_phase_offset' : float
                This is a phase offset between true LFP and measured LFP through electronics equipment that is subtracted from the
                determined fluorescence phase.
            'fs': float
                Sampling frequency in [Hz].
            'loc': optional, 1D numpy.ndarray
                Locomotion signal, i.e. sustained movement, during which phase is calculated.
            'still': optional, 1D numpy.ndarray
                Stillness signal, i.e. no movement, during which phase is calculated.
    
    phi_init : float
        Initial value for phase estimation.

    mvmt : str
        Mouse movement type for which phase is calculated. If '' (default), movement type is ignored and phase is calculated for the entire period.
        If 'loc', phase is calculated during locomotion periods only, and if 'still', phase is calculated during still periods only.

    mint : float
        Minimum total valid sample time used for phase estimation in [s]. If total time is less than this value, function returns None.

    max_amp : float
        Upper bound on % of amplitude change of fluorescence assumed for fitting.

    Returns
    -------
    dict qith phase result keys or empty dict if there is not enough data to extract phase. When phase estimation is available, dict has:
        Analysis result with keys:
        'amp' : float
            pk-pk oscillation amplitude in [%].
        'amp_se' : float
            Oscillation amplitude standard error.
        'phase' : float
            Oscillation phase relative to local field potential (LFP) in [deg]. Zero phase means the fluorescence oscillation through
            is aligned with the LFP through. Phase wraps around at +-180 deg.
        'phase_wrapped' : bool
            True if phase was wrapped to +- 180 deg, False otherwise. Confidence intervals are centered around the wrapped phase, but
            they are noy wrapped and may exceed +- 180 deg.
        'phase_se' : float
            Oscillation phase standard error.
        'sin_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a sinusoidal model explaining the data w.r.t. a constant model.
        'const_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a constant model explaining the data w.r.t. a sinusoidal model.
        'total_valid_time' : float
            Total duration in [s] of valid samples used for phase estimation. This includes non-nan samples and movement-related gating.
        
    """
    # convert to ordered dict so as to not loose LFP and fluorescence pairing
    sig = OrderedDict(sig)
    
    fs = [meas['fs'] for (key, meas) in sig.iteritems()]
    fluo = [meas['fluo'] for (key, meas) in sig.iteritems()]
    lfp_phase_offsets = [(meas['lfp_phase_offset'] if 'lfp_phase_offset' in meas else 0) for (key, meas) in sig.iteritems()]

    # build locomotion or stillness gating signal
    if mvmt == "loc":
        gating = [meas['loc'] for (key, meas) in sig.iteritems() if 'loc' in meas]
        # locomotion gating requested but locomotion data is not available
        if not gating:
           return {}
    elif mvmt == 'still':
        gating = [meas['still'] for (key, meas) in sig.iteritems() if 'still' in meas]
        # still gating requested but stillness data is not available
        if not gating:
           return {}
    elif mvmt == "":
        gating = []
    else:
        raise ValueError("mvmt parameter should be 'loc', 'still' or ''.")
       
    # obtain theta oscillation through times
    theta_peaks = [proc.get_oscillation_extrema(meas['lfp'], meas['fs'], passband = (5,10), bp_order = 10, method = "peak") for (key, meas) in sig.iteritems()]
    # assign LFP phase to fluorescence samples and exclude nan values
    fl_phase, fl_int, total_valid_time = proc.assign_phase(phase_reset_times = theta_peaks, fs = fs, sig = fluo, phase_offsets = lfp_phase_offsets, gating = gating)
    # make sure phase range is between -180 and 180 deg with 0 deg marking LFP troughs
    fl_phase = fl_phase - 180
    
    if total_valid_time < mint:
        return {}
    # estimate fluorescence oscillation phase using a nonlinear least squares fit to a sinusoidal model
    model_sin2 = lmfitModel(sinfit2)
    model_sin2.set_param_hint('offset', value = 100)
    model_sin2.set_param_hint('amp_cos', value = 0, min = -max_amp, max = max_amp)
    model_sin2.set_param_hint('amp_sin', value = 1, min = -max_amp, max = max_amp)
    try:
        result_sin2 = model_sin2.fit(fl_int, x = fl_phase)
    except ValueError:
        util.clrd_print('Cannot calculate phase for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}
    # calculate confidence intervals
    # don't know yet how to make use of confidence intervals when converting to a single amplitude and phase variable.
    #try:
    #    ci = result_sin2.conf_interval(maxiter = 10000)
    #except:
    #    util.clrd_print('Cannot determine confidence intervals.\n', 'warn')
    #    ci = None 
    # also do a nonlinear least squares fit to a constant model to compare Akaike information criteria and justify sinusoidal model selection over constant value
    model_const = lmfitModel(constfit)
    try:
        result_const = model_const.fit(fl_int, x = fl_phase, m = 1)
    except ValueError:
        util.clrd_print('Cannot calculate constant model for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}
    # output result
    out = {}
    # amplitude best fit
    amp_sin = result_sin2.best_values['amp_sin']
    amp_cos = result_sin2.best_values['amp_cos']
    phi = -math.atan2(amp_cos,amp_sin) # added a minus sign to define positive phases as a delay of fluorescence waveform w.r.t LFP.
    out['amp'] = amp_sin/math.cos(phi)
    # this is never needed as out['amp'] is always positive regardless of the sign of amp_sin and amp_cos
    # make amplitude positive and adjust phase since -a*sin(b) = a*sin(b+pi)
    #if out['amp'] < 0:
    #    out['amp'] = -out['amp']
    #    phi = phi+np.pi
    # amplitude standard error, using error propagation
    partial_amp_ampsin = (amp_cos**2/amp_sin**2+1)**-0.5
    partial_amp_ampcos = amp_cos/amp_sin*partial_amp_ampsin
    out['amp_se'] = ((partial_amp_ampcos*result_sin2.params['amp_cos'].stderr)**2+(partial_amp_ampsin*result_sin2.params['amp_sin'].stderr)**2)**0.5

    # phase best fit w.r.t. LFP
    # note: a positive phase means LFP is delayed w.r.t to the fluorescence signal
    out['phase'] = phi*180/np.pi
    out['phase_wrapped'] = False
    # standard error of phase estimate
    partial_phi_ampcos = -amp_sin/(amp_sin**2+amp_cos**2)
    partial_phi_ampsin = amp_cos/(amp_sin**2+amp_cos**2)
    out['phase_se'] = 180/np.pi*((partial_phi_ampcos*result_sin2.params['amp_cos'].stderr)**2+(partial_phi_ampsin*result_sin2.params['amp_sin'].stderr)**2)**0.5
    
    # using Akaike's information criterion, determine the chance that the provided data can be better explained by a constant model compared to a sinusoidal model
    min_aic = min(result_sin2.aic, result_const.aic)
    out['sin_model_rel_likelihood'] = np.exp((min_aic-result_sin2.aic)/2.)
    out['const_model_rel_likelihood'] = np.exp((min_aic-result_const.aic)/2.)
    out['total_valid_time'] = total_valid_time

    return out

# warning: results do not agree with v2, likely because of phase measurement method, i.e. here a hilbert transform is used for LFP phase
# measurement in the theta band
def meas_LFP_fluorescence_phase_v3(sig, phi_init = 0, mvmt = '', mint = 10, max_amp = 20):
    """
    Measures theta-band fluorescence phase w.r.t. local field potential (LFP). Second implementation using sinfit2 that is a convex optimization problem.

    Parameters
    ----------
    sig: dict of dict
        Measurements to combine, with dict labels indicating measurement names and values being dict with labels:
            'fluo': 1D numpy.ndarray
                Fluorescence trace.
            'lfp': 1D numpy.ndarray
                Local field potential.
            'fs': float
                Sampling frequency in [Hz].
            'loc': optional, 1D numpy.ndarray
                Locomotion signal, i.e. sustained movement, during which phase is calculated.
            'still': optional, 1D numpy.ndarray
                Stillness signal, i.e. no movement, during which phase is calculated.
    
    phi_init : float
        Initial value for phase estimation.

    mvmt : str
        Mouse movement type for which phase is calculated. If '' (default), movement type is ignored and phase is calculated for the entire period.
        If 'loc', phase is calculated during locomotion periods only, and if 'still', phase is calculated during still periods only.

    mint : float
        Minimum total valid sample time used for phase estimation in [s]. If total time is less than this value, function returns None.

    max_amp : float
        Upper bound on % of amplitude change of fluorescence assumed for fitting.

    Returns
    -------
    dict qith phase result keys or empty dict if there is not enough data to extract phase. When phase estimation is available, dict has:
        Analysis result with keys:
        'amp' : float
            pk-pk oscillation amplitude in [%].
        'amp_se' : float
            Oscillation amplitude standard error.
        'phase' : float
            Oscillation phase relative to local field potential (LFP) in [deg]. Zero phase means the fluorescence oscillation through
            is aligned with the LFP through. Phase wraps around at +-180 deg.
        'phase_wrapped' : bool
            True if phase was wrapped to +- 180 deg, False otherwise. Confidence intervals are centered around the wrapped phase, but
            they are noy wrapped and may exceed +- 180 deg.
        'phase_se' : float
            Oscillation phase standard error.
        'sin_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a sinusoidal model explaining the data w.r.t. a constant model.
        'const_model_rel_likelihood' : float
            Akaike's information criterion based relative likelihood of a constant model explaining the data w.r.t. a sinusoidal model.
        'total_valid_time' : float
            Total duration in [s] of valid samples used for phase estimation. This includes non-nan samples and movement-related gating.
        
    """
    # convert to ordered dict so as to not loose LFP and fluorescence pairing
    sig = OrderedDict(sig)
    
    #fs = [meas['fs'] for (key, meas) in sig.iteritems()]
    #fluo = [meas['fluo'] for (key, meas) in sig.iteritems()]

    # build locomotion or stillness gating signal
    if mvmt == "loc":
        for meas in sig.itervalues():
            if "loc" in meas:
                meas["gating"] = meas["loc"]
            else:
                return {}
        #gating = [meas['loc'] for (key, meas) in sig.iteritems() if 'loc' in meas]
        # locomotion gating requested but locomotion data is not available
        #if not gating:
        #    return {}
    elif mvmt == 'still':
        for meas in sig.itervalues():
            if "still" in meas:
                meas["gating"] = meas["still"]
            else:
                return {}
        #gating = [meas['still'] for (key, meas) in sig.iteritems() if 'still' in meas]
        # still gating requested but stillness data is not available
        #if not gating:
        #    return {}
    elif mvmt == "":
        #gating = []
        pass
    else:
        raise ValueError("mvmt parameter should be 'loc', 'still' or ''.")
       
    # obtain theta oscillation through times
    #theta_throughs = [proc.get_oscillation_throughs(meas['lfp'], meas['fs'], passband = (5,10)) for (key, meas) in sig.iteritems()]
    # assign LFP phase to fluorescence samples and exclude nan values
    #fl_phase, fl_int, total_valid_time = proc.assign_phase(fluo, fs, theta_throughs, gating)

    fl_phase, fl_int, total_valid_time = proc.assign_phase_v2(sigs = sig.itervalues(), fband = (5,10), filt_order = 2)
    
    if total_valid_time < mint:
        return {}
    # estimate fluorescence oscillation phase using a nonlinear least squares fit to a sinusoidal model
    model_sin2 = lmfitModel(sinfit2)
    model_sin2.set_param_hint('offset', value = 100)
    model_sin2.set_param_hint('amp_cos', value = 1, min = -max_amp, max = max_amp)
    model_sin2.set_param_hint('amp_sin', value = 0, min = -max_amp, max = max_amp)
    try:
        result_sin2 = model_sin2.fit(fl_int, x = fl_phase)
    except ValueError:
        util.clrd_print('Cannot calculate phase for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}
    # calculate confidence intervals
    # don't know yet how to make use of confidence intervals when converting to a single amplitude and phase variable.
    #try:
    #    ci = result_sin2.conf_interval(maxiter = 10000)
    #except:
    #    util.clrd_print('Cannot determine confidence intervals.\n', 'warn')
    #    ci = None 
    # also do a nonlinear least squares fit to a constant model to compare Akaike information criteria and justify sinusoidal model selection over constant value
    model_const = lmfitModel(constfit)
    try:
        result_const = model_const.fit(fl_int, x = fl_phase, m = 1)
    except ValueError:
        util.clrd_print('Cannot calculate constant model for signals (group): {}.'.format(', '.join(sig.keys())), 'error')
        return {}
    # output result
    out = {}
    # amplitude best fit
    amp_sin = result_sin2.best_values['amp_sin']
    amp_cos = result_sin2.best_values['amp_cos']
    phi = math.atan2(amp_cos,amp_sin)
    out['amp'] = amp_sin/math.cos(phi)
    # make amplitude positive and adjust phase since -a*sin(b) = a*sin(b+pi)
    if out['amp'] < 0:
        out['amp'] = -out['amp']
        phi = phi+np.pi
    # amplitude standard error, using error propagation
    partial_amp_ampsin = (amp_cos**2/amp_sin**2+1)**-0.5
    partial_amp_ampcos = amp_cos/amp_sin*partial_amp_ampsin
    out['amp_se'] = ((partial_amp_ampcos*result_sin2.params['amp_cos'].stderr)**2+(partial_amp_ampsin*result_sin2.params['amp_sin'].stderr)**2)**0.5

    # phase best fit w.r.t. LFP
    out['phase'] = phi*180/np.pi
    out['phase_wrapped'] = False
    # standard error of phase estimate
    partial_phi_ampcos = amp_sin/(amp_sin**2+amp_cos**2)
    partial_phi_ampsin = -amp_cos/(amp_sin**2+amp_cos**2)
    out['phase_se'] = 180/np.pi*((partial_phi_ampcos*result_sin2.params['amp_cos'].stderr)**2+(partial_phi_ampsin*result_sin2.params['amp_sin'].stderr)**2)**0.5
    
    # using Akaike's information criterion, determine the chance that the provided data can be better explained by a constant model compared to a sinusoidal model
    min_aic = min(result_sin2.aic, result_const.aic)
    out['sin_model_rel_likelihood'] = np.exp((min_aic-result_sin2.aic)/2.)
    out['const_model_rel_likelihood'] = np.exp((min_aic-result_const.aic)/2.)
    out['total_valid_time'] = total_valid_time

    return out

def meas_transient_events(sigs, flchan, tchan, excl_chan, peak_perc = 10, evt_direction = 'decrease', bs_alpha = 0.05, evt_dur_zscore = 1,
    evt_pk_zscore = 2, evt_dur_zscore_w = 5, avg_evt_waveform_dur_perc = 100, avg_evt_waveform_dur_perc_rng = 'lower', 
    evt_waveform_avg_window = 0.5, nbs = 1000, return_raw_meas = False):
    """
    Estimates transient amplitude from noisy fluorescence recordings.
    
    Parameters
    ----------
    sigs : list of tuple
        List of tuple pairs of dict signal channel and sampling frequnecy in [Hz]
        as (<signal channels>, <fs>).

    flchan : str
        Name of fluorescence channel to use for event measurement. E.g. 'Ch1' or 'Ch1/Ch2'.

    tchan : str
        Name of channel containing detected fluorescence transients E.g. 'Ch1/Ch2_tdetect'. Channel must be 1D numpy.ndarray of bool type.

    excl_chan : list of str
        Name of binary mask channels used to exclude samples from the average waveform transient analysis. Note that these are not applied again to
        the transient channels since transients must be free of artifacts and precomputed.
        
    peak_perc : float
        Estimates event amplitude using a percentile method. For events that decrease from baseline, the median of the lower percentile is calculated, otherwise, for
        events increasing from baseline, the median of the upper percentile is calculated. 
            
    alpha : float
        Confidence interval alpha level.
        
    evt_dur_zscore : float
        Event FWHM duration measurement z-score level for FWHM threshold to be crossed and trigger width measurement.

    evt_pk_zscore : float
        Event peak (median of lower or upper peak_perc) z-score threshold. Sign is ignored.
    
    evt_dur_zscore_w : float
        Event duration estimation moving-zscore window size in [s] if non-zero sampling frequency is provided with a signal
        or number of samples otherwise. This needs to be adjusted sensibly based on the bleaching time-scale.
        Set to 0 if no moving window should be used.

    avg_evt_waveform_dur_perc : float
        Event waveform averaging FWHM duration percentile that determines which events are to be averaged. Defaults to 10th percentile, i.e. fast events.

    evt_waveform_avg_window : float
        Duration of average waveform window in [s] to use for aligning events.

    nbs : int
        Number of iterations to use for bootstrap median estimation.

    return_raw_meas : bool
        If True, output will contain also raw measurement vectors of event amplitude, duration, etc.
            
    Notes
    -----
    - input signal baseline must be 1 and events normalized to baseline.
    - using 'min'/'max' tends to overestimate a decrease/increase in signal in a SNR dependent way.
    - using percentiles tends to underestimate amplitudes in a SNR dependent way for most part, and when SNR is very low, underestimates.
    
    Returns
    -------
    dict with keys
        'alpha' : float
            Alpha value chosen for confidence interval calculation.
        'peak_perc' : float or str
            Peak measurement lower or upper percentile.
        'amp_median' : float
            Mean event amplitude.
        'amp_median_lower_bound', 'amp_median_upper_bound' : float
            Mean event amplitude lower and upper confidence intervals with given alpha.
        'amp_values': 1D numpy.ndarray of event amplitudes that can be used to e.g. create a histogram if return_raw_meas is True.
        'dur_median' : float
            Mean event duration in [s] (if fs !=0) otherwise in samples.
        'dur_median_lower_bound', 'amp_median_upper_bound' : float
            Mean event duration lower and upper confidence intervals with given alpha.
        'dur_values' : 1D numpy.ndarray of event durations that can be used to e.g. create a histogram if return_raw_meas is True.
        'int_median' : float
            Mean event integral in [1*s] (if fs !=0) otherwise in samples.
        'int_median_lower_bound', 'amp_median_upper_bound' : float
            Mean event integral lower and upper confidence intervals with given alpha.
        'int_values' : 1D numpy.ndarray 
            Event integrals that can be used to e.g. create a histogram if return_raw_meas is True.
        'total_time' : float or None
            Total time considered for event analysis in [s] if fs provided, None otherwise.
        'avg_evt_waveform': 1D numpy.ndarray
            Average event waveform.
        'avg_evt_waveform_fs': float
            TEMPORARY: average sampling frequency of all sweeps; slight differences may arrise so it's best to make this more general.
    """
    # check input
    assert avg_evt_waveform_dur_perc_rng in ['lower', 'upper']
    
    out = {}
    event_amp = [] # collect estimated event amplitudes
    event_dur = [] # event FWHM duration in [s] if fs !=0 or number of samples otherwise
    peak_idx = [] # FWHM center timepoint to use for event alignment.
    tmp_event_dur = [] # temporary array to keep track of event durations used together with fwhm_center_timepoint 
    event_int = [] # event integral [1*s] if fs !=0 or [1] w.r.t baseline.
    total_time = 0 # total time considered for event analysis in [s].
    all_fs = [] # collect all sampling frequencies to average later
    
    # make a copy of the fluorescence signal and set exclude samples to nan
    flsigs = []
    for sig_idx, sig in enumerate(sigs):
        flsigs.append(np.copy(sig[0][flchan]))
        for exchan_name in excl_chan:
            flsigs[-1][sig[0][exchan_name]] = np.nan
    
    for sig_idx, sig in enumerate(sigs):

        peak_idx.append([])
        tmp_event_dur.append([])
        
        flsig = flsigs[sig_idx]
        evtsig = sig[0][tchan]

        # sampling frequency in [Hz]
        fs = sig[1]
        if fs is None or not fs:
            fs = 0
        else:
            all_fs.append(fs)
        
        if fs:
            nw = int(evt_dur_zscore_w*fs)
        else:
            nw = int(evt_dur_zscore_w)

        # add up total valid waveform time
        if fs and total_time is not None:
            total_time += np.count_nonzero(~np.isnan(flsig))/fs
        else:
            total_time = None

        # calculate a moving standard deviation to be used for each sample in the signal; use last valid std for last samples
        # -------------------------------------------------------------------------------------------------------------------
        # if signal is shorter than the window, do not use a moving window
        
        #if not evt_dur_zscore_w or evt_dur_zscore_w and len(flsig)<nw: 
        #    wstd = np.full((len(flsig),), np.nanstd(flsig))
        #else:
        #    wstd = np.zeros((len(flsig),))
        #    for i in range(len(flsig)-nw):
        #        wstd[i] = np.nanstd(flsig[i:i+nw])
        #    for i in range(len(flsig)-nw, len(flsig)):
        #        wstd[i] = np.nanstd(flsig[len(flsig)-nw:])
        
        flsig_no_evt = flsig[:len(evtsig)][~evtsig]
        if not evt_dur_zscore_w or evt_dur_zscore_w and len(flsig)<nw: 
            wstd = np.full((len(flsig),), np.nanstd(flsig_no_evt))
        else:
            wstd = np.zeros((len(flsig),))
            for i in range(len(flsig_no_evt)-nw):
                wstd[i] = np.nanstd(flsig_no_evt[i:i+nw])
            for i in range(len(flsig_no_evt)-nw, len(flsig)):
                wstd[i] = np.nanstd(flsig_no_evt[len(flsig_no_evt)-nw:])
        
        

        # estimate amplitude in each event interval
        # get indices for False->True transitions, index indicates the position of the False element
        t = (evtsig[:-1] < evtsig[1:]).nonzero()[0]
        if evtsig[0]:
            t = np.insert(t, 0, 0)
        # select intervals
        for t_idx in t:
            # get interval width
            w = 0
            while t_idx+w+1 < len(evtsig) and evtsig[t_idx+w+1]:
                w += 1
            # alias for part of fluorescence within the event window
            flsig_window = flsig[t_idx+1:t_idx+1+w]
            
            # calculate event peak value and index
            if evt_direction == 'decrease':
                sig_window_perc = np.nanpercentile(flsig_window, peak_perc)
                sig_perc_crop = flsig_window<=sig_window_perc
                fl_cropped_values = flsig_window[sig_perc_crop]
                
                # store event amplitude
                evt_amp_median = np.nanmedian(fl_cropped_values)
                
                # skip if significance cannot be achieved
                if evt_amp_median > 1-abs(evt_pk_zscore)*np.nanmean(wstd[t_idx+1:t_idx+1+w]):
                    continue
                
                
                # not sure why, but turns out this gives a bad waveform average
                #if min(flsig_window) > 1-2*np.nanmean(wstd[t_idx+1:t_idx+1+w]):
                #    continue
                    
                event_amp.append(evt_amp_median)
                
                # store event peak centroid index
                fl_cropped_indices = np.nonzero(sig_perc_crop)[0]
                sig_perc_wnd_sum = np.sum(sig_window_perc-fl_cropped_values)
                if sig_perc_wnd_sum:
                    # calculate center of mass index if all samples differ from the percentile
                    peak_idx[sig_idx].append(t_idx+1+int(np.sum(fl_cropped_indices*(sig_window_perc-fl_cropped_values))/sig_perc_wnd_sum))
                else:
                    # just get the center index
                    peak_idx[sig_idx].append(t_idx+1+int((fl_cropped_indices[-1]+fl_cropped_indices[0])/2.))
                            
            elif evt_direction == 'increase':
                sig_window_perc = np.nanpercentile(flsig_window, 100-peak_perc)
                sig_perc_crop = flsig_window>=sig_window_perc
                fl_cropped_values = flsig_window[sig_perc_crop]
                
                # store event amplitude (add option to calculate z-score for significance)
                evt_amp_median = np.nanmedian(fl_cropped_values)
                
                # skip if significance cannot be achieved
                if evt_amp_median < 1+abs(evt_pk_zscore)*np.nanmean(wstd[t_idx+1:t_idx+1+w]):
                    continue
                
                event_amp.append(evt_amp_median)
                
                # store event peak centroid index
                fl_cropped_indices = np.nonzero(sig_perc_crop)[0]
                sig_perc_wnd_sum = np.sum(fl_cropped_values-sig_window_perc)
                if sig_perc_wnd_sum:
                    # calculate center of mass index if all samples differ from the percentile
                    peak_idx[sig_idx].append(t_idx+1+int(np.sum(fl_cropped_indices*(fl_cropped_values-sig_window_perc))/sig_perc_wnd_sum))
                else:
                    # just get the center index
                    peak_idx[sig_idx].append(t_idx+1+int((fl_cropped_indices[-1]+fl_cropped_indices[0])/2.))
            
            # threshold all samples within the window around the FWHM amplitude value and get indices where this is true
            if evt_amp_median >= 1:
                dur_thresh_idx = np.where((flsig[t_idx+1:t_idx+1+w] > (1+evt_amp_median)/2.+abs(evt_dur_zscore)*np.nanmean(wstd[t_idx+1:t_idx+1+w]))==True)[0]
            else:
                dur_thresh_idx = np.where((flsig[t_idx+1:t_idx+1+w] < (1+evt_amp_median)/2.-abs(evt_dur_zscore)*np.nanmean(wstd[t_idx+1:t_idx+1+w]))==True)[0]
            
            # calculate event duration
            if len(dur_thresh_idx):
                duration_idx_diff = max(1, dur_thresh_idx[-1]-dur_thresh_idx[0]) # if indices are identical, then make the event at least one sample
                if fs:
                    event_dur.append(duration_idx_diff/fs)
                else:
                    event_dur.append(duration_idx_diff)
                tmp_event_dur[sig_idx].append(event_dur[-1])
            else:
                # for these events, duration could not be determined; it is important to keep peak_idx and event_dur aligned for later doing waveform averaging
                tmp_event_dur[sig_idx].append(0)
                
            # calculate event integral
            if fs:
                event_int.append(np.nansum((flsig[t_idx+1:t_idx+1+w]-1)/fs))
            else:
                event_int.append(np.nansum(flsig[t_idx+1:t_idx+1+w]-1))
            
    event_amp = np.array(event_amp)
    event_dur = np.array(event_dur)
    event_int = np.array(event_int)
    all_fs = np.array(all_fs)
    
    # WARNING!!!!!!!!!!!!!!!!!!!!!!!
    # MUST FIX!!!!!!!!!!!!!!!!!!!!!!
    # not sure why, the number of counted navg_traces for waveform averaging is different for the same 25 lower and 25 upper percentile!!!! WHY?
    
    # calculate event duration percentile to filter event
    lower_perc_dur = np.percentile(event_dur, avg_evt_waveform_dur_perc)
    upper_perc_dur = np.percentile(event_dur, 100-avg_evt_waveform_dur_perc)
    # average event waveforms by centering on the FWHM center timepoint
    # TEMPORARY!!!!!!!!!!!!!!!!!!
    # For now just assume all sampling frequencies are the same, and just pick out a slice from each sweep and average together.
    avg_fs = np.mean(all_fs)
    round_up_odd = lambda num: num - (num%2)+1
    n_avg_evt_waveform_samples = round_up_odd(int(avg_fs*evt_waveform_avg_window)) # define the middle sample of the average waveform as the one centered on the event peaks
    avg_evt_waveform = np.full((n_avg_evt_waveform_samples,), np.nan) # needed to keep track of whether signal samples were added, if initializing to 0, can distort signal.
    navg_samples = np.zeros((n_avg_evt_waveform_samples,))
    navg_traces = 0 # keep track of how many sweeps have been included in the average
    for sig_idx, flsig in enumerate(flsigs): 
        # add signals
        for evt_idx, evt_dur in enumerate(tmp_event_dur[sig_idx]):
            if evt_dur and (avg_evt_waveform_dur_perc_rng == 'lower' and evt_dur <= lower_perc_dur or avg_evt_waveform_dur_perc_rng == 'upper' and evt_dur >= upper_perc_dur): 
                # pick out waveform
                left_crop = peak_idx[sig_idx][evt_idx]-int((n_avg_evt_waveform_samples-1)/2)
                right_crop = peak_idx[sig_idx][evt_idx]+int((n_avg_evt_waveform_samples-1)/2)+1
                if left_crop >= 0 and right_crop <= len(flsig):
                    samples_to_add = ~np.isnan(flsig[left_crop:right_crop]) # keep track of valid samples in each time bin
                    navg_samples += samples_to_add
                    if np.any(samples_to_add):
                        navg_traces += 1
                    avg_evt_waveform = np.nansum(np.dstack((avg_evt_waveform,flsig[left_crop:right_crop])),2)[0] # this summation ignores nan, i.e. element value added to nan is same.
        
        
    avg_evt_waveform = np.divide(avg_evt_waveform, navg_samples)
    # bootstrap event peaks, durations and integrals
    bs_amp_result = bs.bootstrap(event_amp, stat_func = bs_stats.median, num_iterations = nbs, num_threads = 10, alpha = bs_alpha)
    bs_dur_result = bs.bootstrap(event_dur, stat_func = bs_stats.median, num_iterations = nbs, num_threads = 10, alpha = bs_alpha)
    bs_int_result = bs.bootstrap(event_int, stat_func = bs_stats.median, num_iterations = nbs, num_threads = 10, alpha = bs_alpha)
    # pack results
    out['alpha'] = bs_alpha
    out['peak_perc'] = peak_perc
    out['amp_median'] = bs_amp_result.value
    out['amp_median_lower_bound'] = bs_amp_result.lower_bound
    out['amp_median_upper_bound'] = bs_amp_result.upper_bound
    out['dur_median'] = bs_dur_result.value
    out['dur_median_lower_bound'] = bs_dur_result.lower_bound
    out['dur_median_upper_bound'] = bs_dur_result.upper_bound
    out['int_median'] = bs_int_result.value
    out['int_median_lower_bound'] = bs_int_result.lower_bound
    out['int_median_upper_bound'] = bs_int_result.upper_bound
    out['total_time'] = total_time
    out['avg_evt_waveform'] = avg_evt_waveform
    out['avg_evt_waveform_fs'] = avg_fs
    out['avg_evt_waveform_ntraces'] = navg_traces
    out['avg_evt_waveform_dur_perc'] = avg_evt_waveform_dur_perc
    out['avg_evt_waveform_dur_perc_rng'] = avg_evt_waveform_dur_perc_rng
    

    if return_raw_meas:
        out['amp_values'] = event_amp
        out['dur_values'] = event_dur
        out['int_values'] = event_int

    return out

def ripple_power(sig, fs, passband = [80, 250], fc = 25, clip = 0, return_ripple_band_LFP = False):
    """
    Generates time-varying ripple power from local-field potential signal (LFP).

    Parameters
    ----------
    sig : array_like
        Local field potential signal.
        
    fs : float
        Sampling frequency in [Hz].

    passband : tuple of float
        Ripple pass-band frequencies as (lower, upper) in [Hz].
        
    clip : float
        Ripple band signal clip, as S.D.

    return_ripple_band_LFP : bool
        If True, returns also the ripple band-passed LFP.
        
    Returns
    -------
    - if return_ripple_band_LFP is True:
        tuple of (1D numpy.ndarray, 1D numpy.ndarray)
            Time-varying ripple power and ripple-band LFP.
    - if return_ripple_band_LFP is False:
        1D numpy.ndarray
            Time-varying ripple power
    """
    ripple_band_LFP = proc.filt_butter_bp(sig = sig, fs = fs, lowcut = passband[0], highcut = passband[1], order = 2)
    #fc = (passband[0]+passband[1])/(2*np.pi)

    # clip extreme values
    std_ripple = np.std(ripple_band_LFP)
    if clip:
        ripple_band_LFP[ripple_band_LFP > clip*std_ripple] = clip*std_ripple
        ripple_band_LFP[ripple_band_LFP < -clip*std_ripple] = -clip*std_ripple          
    
    # instantaneous ripple power
    ripple_power = ripple_band_LFP**2
    # integrate/low-pass filter ripple power
    b, a = scipysig.butter(2, fc/(fs/2), 'low', analog = False) # filtfilt will double the order

    if return_ripple_band_LFP:
        return scipysig.filtfilt(b, a, ripple_power), ripple_band_LFP
    else:
        return scipysig.filtfilt(b, a, ripple_power)
    
def get_ripple_events(lfp, fs, gating = None, min_valid_gating_time = 2, iqr_score_level = 10, passband = [80, 250], ignore_dur = 15, join_dur = 15):
    """
    Returns binarized ripple epochs.
    
    Parameters
    ----------
    lfp : array_like
        Local field potential signal.
        
    fs : float
        Sampling frequency in [Hz].

    gating : None or 1D numpy bool array
        Gating signal during which ripples can be detected, e.g. during still periods marked by True values.

    min_valid_gating_time : float
        If gating provided, minimum total valid gating time in [s] to consider for calculating reliably Tukey's inter-quartile range score.

    iqr_score_level : float
        Ripple detection significance level for time-varying power within the ripple band. Ripple power is scored
        using Tukey's inter-quartile range score.
        
    passband : tuple of float
        Ripple pass-band frequencies as (lower, upper) in [Hz].
        
    ignore_dur : float
        Minimum ripple epoch duration in [ms] to consider as ripple.
        
    join_dur : float
        If duration between ripple epochs is smaller than this value, then join epochs.
        
    Returns
    -------
    None if cannot calculate ripple power dur to sampling rate and ripple band mismatch.

    tuple (df: pandas DataFrame, epoch: 1D bool numpy.ndarray, rp: 1D numpy.ndarray)
        df : pandas Dataframe
            Ripple interval info with ripple number on rows and columns:
                'interval' : tuple
                    Ripple epoch start and end time (t_start, t_end) tuple in [s].
                'duration' : float
                    Ripple epoch duration in [s].
                'peak-time' : float
                    Ripple peak time defined as the deepest ripple through within the ripple interval.
                'peak-power' : float
                    Ripple peak power in voltage squared units of lfp.

        epoch : 1D bool numpy.ndarray
            Ripple event binary segmentation signal, with True marking ripple period.

        rp : 1D numpy.ndarray
            Instantaneous ripple power.
        
    """
    empty_df = pd.DataFrame(columns = ['interval', 'duration', 'peak-time', 'peak-power'])
    empty_rp_epochs = np.full((len(lfp),), False)

    # if gating provided, check if minimum gating time condition is satisfied
    if gating is not None:
        gating = gating.astype(bool)
        if np.sum(gating)/fs < min_valid_gating_time:        
            return empty_df, empty_rp_epochs, np.array([])

    # null NaNs for the band-pass filter to work properly
    lfp[np.isnan(lfp)] = 0
    # try to obtain time-varying ripple power, if not able to, return None
    try:
        rp, rLFP = ripple_power(sig = lfp, fs = fs, passband = passband, clip = 5, return_ripple_band_LFP = True)
    except:
        return None

    # calculate iqr score during valid periods e.g. during stillness
    if gating is not None:
        rp_copy = np.copy(rp)
        rp_copy[~gating] = np.nan
        rp_iqr_score = iqr_score(rp_copy)
    else:
        rp_iqr_score = iqr_score(rp)
    
    rp_epochs = rp_iqr_score > iqr_score_level
    # apply gate
    if gating is not None:
        rp_epochs[~gating] = False
    
    # expand epochs until just 2 SD above baseline
    rp_idx = 0
    while rp_idx < len(rp):
        # scan for ripple epoch start index    
        if not rp_epochs[rp_idx]:
            rp_idx += 1
            continue
        # count number of samples in the ripple epoch within valid gating period
        n_rip_samp = 0    
        while rp_idx+n_rip_samp < len(rp) and rp_epochs[rp_idx+n_rip_samp]:
            n_rip_samp += 1                
        # discard epoch if too short
        if n_rip_samp < fs*ignore_dur*1e-3:
            rp_epochs[rp_idx:rp_idx+n_rip_samp] = False
            rp_idx += n_rip_samp
            continue
        
        # expand left until significance level drops below 2 SD
        for left_idx in range(rp_idx-1, 0, -1):
            if rp_iqr_score[left_idx] > 2:
                rp_epochs[left_idx] = True
            else:
                break

        # expand right
        for right_idx in range(rp_idx+n_rip_samp, len(rp)):
            if rp_iqr_score[right_idx] > 2:
                rp_epochs[right_idx] = True
            else:
                break
        # jump over ripple epoch
        rp_idx += n_rip_samp
        
    # join ripple epochs that are too close to each other
    rp_idx = 0
    while rp_idx < len(rp):
        # scan for ripple epoch start index
        if not rp_epochs[rp_idx]:
            rp_idx += 1
            continue
        # count number of samples in the ripple epoch
        n_rip_samp = 0
        while rp_idx+n_rip_samp < len(rp) and rp_epochs[rp_idx+n_rip_samp]:
            n_rip_samp += 1
            
        # measure gap between current and previous ripple epoch
        gap = 0
        while rp_idx-gap and not rp_epochs[rp_idx-1-gap]:
            gap += 1
        # join epochs if gap is too small and not first epoch with no previous epoch
        if gap < fs*join_dur*1e-3 and rp_idx-gap:
            rp_epochs[rp_idx-gap:rp_idx] = True
        
        # jump over ripple epoch
        rp_idx += n_rip_samp    
        
    epoch_intervals = proc.convert_bool_epochs_to_intervals(rp_epochs, fs)
    
    peak_times = []
    peak_powers = []
    for tstart, tend in epoch_intervals:
        start_idx = int(tstart*fs)
        end_idx = int(tend*fs)

        if start_idx == end_idx:
            print(tstart, tend)
        # pick LFP through closest to the peak power
        peak_times.append((start_idx+np.argmin(rLFP[start_idx:end_idx]*rp[start_idx:end_idx]))/fs)
        peak_powers.append(np.max(rp[start_idx:end_idx]))

    df = pd.DataFrame({'interval': epoch_intervals,
                       'duration': [ep_int[1]-ep_int[0] for ep_int in epoch_intervals],
                       'peak-time': peak_times,
                       'peak-power': peak_powers},
                       columns = ['interval', 'duration', 'peak-time', 'peak-power'])

    return df, rp_epochs, rp

def iqr_score(sig):
    """
    Calculates Tukey's interquartile score for clasification of outliers for non-normally distributed data.
    
    Parameters
    ----------
    sig : 1D numpy array
        Signal, can contain nan.
        
    Returns
    -------
    1D numpy array
        IQR-score.
    """
    return (sig - bn.nanmedian(sig))/scipy_stats.iqr(sig, nan_policy = 'omit')    

def nan_chan_excludes(channel, fs):
    """
    Applies exclusion intervals for signals in channels from aligned signals data structures and returns a new signal.

    Parameters
    ----------
    channel : dict
        Channel dict in the aligned signals dict. Minimally a channel has 'sig' and 'analysis' keys.

    fs : float
        Sampling frequency in [Hz].

    Returns
    -------
    1D numpy array or dict of 1D numpy array
        New signal with excluded intervals. Original signal is left unchanged.
    """
    exclude_intervals = []
    if 'analysis' in channel and 'exclude' in channel['analysis']:
        for ex_type, ex_val in iteritems(channel['analysis']['exclude']):
            if ex_val['times']:
                exclude_intervals.extend(ex_val['times'])
                
    if isinstance(channel['sig'], np.ndarray):
        ch_sig_copy = np.copy(channel['sig'].astype(float))
        proc.nan_itervals(ch_sig_copy, fs, exclude_intervals)
    elif isintance(channel['sig'], dict):
        ch_sig_copy = {}
        for sig_name in channel['sig']:
            ch_sig_copy[sig_name] = np.copy(channel['sig'][sig_name].astype(float))
            proc.nan_itervals(ch_sig_copy[sig_name], fs, exclude_intervals)
    else:
        raise TypeError

    return ch_sig_copy
    
def get_treadmill_mvmt(dist, fs, mvmt_vel_threshold = 2, locomotion_vel_threshold = 10, min_locomotion_dur = 1,
    still_settling_time = 1, tavg = 0.5):
    """
    Extracts treadmill movement parameters
    
    Parameters
    ----------
    dist: 1D numpy array
        Distance travelled on belt. Ensure this distance does not wrap around.
        
    fs : float
        Sampling frequency in [Hz].
        
    mvmt_vel_threshold : float
        Belt velocity above which movement is registered in [mm/s]. Below this velocity and after a given still_settling_time, the mouse is considered still.
        
    locomotion_vel_threshold : float
        Belt velocity in [mm/s] above which mouse is considered to be walking/running.
        
    min_locomotion_dur : float
        Minimum duration of sustained locomotion velocity above locomotion_vel_threshold to consider a locomotion epoch.
        
    still_settling_time : float
        Settling time from movement to stillness in [s].
    
    tavg : float
        Moving average window to smooth velocity trace in [s].
        
    Returns
    -------
    out : dict with keys:
      'sig': dict with keys:
        'vel' : 1D numpy float array
            Belt velocity in [mm/s].
        'mvmt': 1D numpy bool masked array
            Bool vector indicating if movement occured.
        'loc': 1D numpy bool array
            Signal marking where locomotion is occuring.
        'still': 1D numpy bool array
            Signal marking where mouse is still.
            
      'events': dict with keys:
        'mvmt_intervals': list of 2 element tuples
            Interval in [s] as (<start time>, <end time>) during which movement occurs.
        'mvmt_idxs': list of 2 element tuples
            Interval as array indices (<start idx>, <end idx>) during which movement occurs.     
        'loc_intervals': list of 2 element tuples
            Locomotion interval in [s] as (<start time>, <end time>). Note, this has a higher threshold than a movement and has to be sustained
            for at least min_locomotion_dur.
        'loc_interval_idxs': list of 2 element tuples
             Locomotion interval as array indices (<start idx>, <end idx>). Note, this has a higher threshold than a movement and has to be sustained
            for at least min_locomotion_dur.
        'still_intervals': list of 2 element tuples
            Interval in [s] as (<start time>, <end time>) during which mouse is considered still. These intervals are shorter than non-movement intervals
            because of a minimum still settling time.
        'still_interval_idxs': list of 2 element tuples
            Still intervals as list of tuples with array indices (<start idx>, <end idx>).
        'max_mvmt_intervals': 2 element tuple
            Longest interval in [s] as (<start time>, <end time>) during which movement occurs.
        'max_mvmt_idxs': 2 element tuple
            Longest interval as array indices (<start idx>, <end idx>) during which movement occurs.
        'max_loc_intervals': 2 element tuple
            Longest locomotion interval in [s] as (<start time>, <end time>). Note, this has a higher threshold than a movement and has to be sustained
            for at least min_locomotion_dur.
        'max_loc_interval_idxs': 2 element tuple
            Longest locomotion interval as array indices (<start idx>, <end idx>). Note, this has a higher threshold than a movement and has to be sustained
            for at least min_locomotion_dur.
        'max_still_intervals': 2 element tuple
            Longest interval in [s] as (<start time>, <end time>) during which mouse is considered still. These intervals are shorter than non-movement intervals
            because of a minimum still settling time.
        'max_still_interval_idxs': 2 element tuple
            Longest still intervals as list of tuples with array indices (<start idx>, <end idx>).

      'analysis': dict with keys:
        'loc_dur': float
            Total duration in [s] for locomotion state.
        'still_dur': float
            Total duration in [s] for stillness state.
        'total_dur': float
            Total duration of behavior observation in [s].
        
    """
    assert mvmt_vel_threshold > 0 and locomotion_vel_threshold > 0 and min_locomotion_dur >= 0 and still_settling_time >= 0 and tavg >= 0
    out = {}
    out['sig'] = {}
    out['events'] = {}
    out['analysis'] = {}
    if tavg:
        navg = int(np.ceil(tavg*fs))
    else:
        navg = 1
    out['sig']['vel'] = proc.mov_avg(np.diff(dist, append = np.nan)*fs, navg, no_delay = True)
    # get movement vectors
    out['sig']['mvmt'] = np.ma.array(np.logical_or(out['sig']['vel']>mvmt_vel_threshold, out['sig']['vel']<-mvmt_vel_threshold), mask = np.isnan(out['sig']['vel']))
    out['events']['mvmt_intervals'] = proc.convert_bool_epochs_to_intervals(out['sig']['mvmt'], fs)
    out['events']['mvmt_interval_idxs'] = proc.convert_bool_epochs_to_intervals(out['sig']['mvmt'])
    # longest movement interval start and end time
    if out['events']['mvmt_interval_idxs']:
        out['events']['max_mvmt_interval'] = out['events']['mvmt_interval_idxs'][np.argmax([i[1]-i[0] for i in out['events']['mvmt_interval_idxs']])]
    # longest movement interval start and end indices
    if out['events']['mvmt_interval_idxs']:
        out['events']['max_mvmt_interval_idxs'] = out['events']['mvmt_interval_idxs'][np.argmax([i[1]-i[0] for i in out['events']['mvmt_interval_idxs']])]
    # get still vectors
    out['events']['still_intervals'] = [(interval[0]+still_settling_time, interval[1]) for interval in proc.convert_bool_epochs_to_intervals(~out['sig']['mvmt'], fs) if interval[0]+still_settling_time < interval[1]]
    out['sig']['still'] = proc.intervals_to_bool(out['events']['still_intervals'], fs, len(dist))
    out['events']['still_interval_idxs'] = proc.convert_bool_epochs_to_intervals(out['sig']['still'])
    # longest still interval start and end time
    if out['events']['still_intervals']:
        out['events']['max_still_interval'] = out['events']['still_intervals'][np.argmax([i[1]-i[0] for i in out['events']['still_intervals']])]
    # longest still interval start and end indices
    if out['events']['still_interval_idxs']:
        out['events']['max_still_interval_idxs'] = out['events']['still_interval_idxs'][np.argmax([i[1]-i[0] for i in out['events']['still_interval_idxs']])]
    # get locomotion bool vector and convert to time intervals to check minimum duration
    out['events']['loc_intervals'] = [interval for interval in proc.convert_bool_epochs_to_intervals(np.ma.array(out['sig']['vel']>locomotion_vel_threshold,
                                      mask = np.isnan(out['sig']['vel'])), fs) if interval[1]-interval[0]>min_locomotion_dur]
    out['sig']['loc'] = proc.intervals_to_bool(out['events']['loc_intervals'], fs, len(dist))
    out['events']['loc_interval_idxs'] = proc.convert_bool_epochs_to_intervals(out['sig']['loc'])
    # longest locomotion interval start and end time
    if out['events']['loc_intervals']:
        out['events']['max_loc_interval'] = out['events']['loc_intervals'][np.argmax([i[1]-i[0] for i in out['events']['loc_intervals']])]
    # longest locomotion interval start and end indices
    if out['events']['loc_interval_idxs']:
        out['events']['max_loc_interval_idxs'] = out['events']['loc_interval_idxs'][np.argmax([i[1]-i[0] for i in out['events']['loc_interval_idxs']])]

    # calculate durations for locomotion and stillness behavior
    out['analysis']['loc_dur'] = np.sum(out['sig']['loc'])/fs
    out['analysis']['still_dur'] = np.sum(out['sig']['still'])/fs
    out['analysis']['total_dur'] = len(dist)/fs
    
    return out

def behavioral_psd(sig, psd_sig_path, fres = 0.5, norm = False):
    """
    Generates an average behavioral-dependent power spectral density of a selected signal, e.g. LFP using Welch's method.
    Two behavioral states are considered, which are read out from the treadmill: stillness and locomotion.

    Note: the sampling rate for the treadmill data is assumed to be the same as the sampling rate of the signal for which the PSD is calculated.

    Parameters
    ----------
    sig : dict
        Aligned signals.

    psd_sig_path : str
        Nested dict path within the provided aligned signals leading to a 1D waveform e.g. local field potential (LFP), fluorescence.
        E.g. fpr LFP, use 'ephys/LFP/sig'.

    fres : float
        Frequency resolution of the calculated PSD.

    norm : bool
        Normalize 1/f^2 PSD by applying a differentiator.

    Returns
    -------
    dict with keys:
        fbins : float
            PSD frequency bins in [Hz].
        loc : 1D numpy array
            PSD during locomotion epochs.
        still : 1D numpy array
            PSD during still epochs.
    """
    fbins = np.array([])
    PSDs_loc= []
    PSDs_still = []
    

    # PSDs for still intervals
    if util.chk_dict_path(sig, 'behavior/treadmill/events/still_interval_idxs'):
        nbins = round(sig['behavior']['treadmill']['fs']/(2.*fres))
        for i in util.get_dpath_val(sig, 'behavior/treadmill/events/still_interval_idxs'):
            if len(util.get_dpath_val(sig, psd_sig_path)[i[0]:i[1]]) >= nbins:
                w = util.get_dpath_val(sig, psd_sig_path)[i[0]:i[1]]
                if norm:
                    w = np.insert(np.diff(w),0,0)
                fbins, psd = scipysig.welch(w, sig['behavior']['treadmill']['fs'], nperseg = nbins, average = 'median')
                PSDs_still.append(psd)
    # PSDs for locomotion intervals
    if util.chk_dict_path(sig, 'behavior/treadmill/events/loc_interval_idxs'):
        nbins = round(sig['behavior']['treadmill']['fs']/(2.*fres))
        for i in util.get_dpath_val(sig, 'behavior/treadmill/events/loc_interval_idxs'):
            if len(util.get_dpath_val(sig, psd_sig_path)[i[0]:i[1]]) >= nbins:
                w = util.get_dpath_val(sig, psd_sig_path)[i[0]:i[1]]
                if norm:
                    w = np.insert(np.diff(w),0,0)
                fbins, psd = scipysig.welch(w, sig['behavior']['treadmill']['fs'], nperseg = nbins, average = 'median')
                PSDs_loc.append(psd)

    out = {}
    out['fbins'] = fbins
    if PSDs_loc:
        out['loc'] =  np.mean(PSDs_loc, axis = 0)
    else:
        out['loc'] = np.array([])

    if PSDs_still:
        out['still'] = np.mean(PSDs_still, axis = 0)
    else:
        out['still'] = np.array([])

    return out
    
def segment_cusum_peaks(fwd_cusum, bkwd_cusum, peak_height, peak_prominence = 1):
    """
    Find events based on forward and backard cusum signals. Since cusum is delayed from an event peak in a way that depends
    on the event amplitude w.r.t. estimated event size, thresholding both forward and backward cusums can be used to construct a window
    where the event peak will be.

    Returns
    -------
    tuple of:
        1D numpy array of bool
            Bool vector of detected events.
    """

    # forward and backward peaks
    assert len(fwd_cusum) == len(bkwd_cusum)
    nsamp = len(fwd_cusum)
    fwd_peaks, _ = scipysig.find_peaks(np.abs(fwd_cusum), height = abs(peak_height), prominence = abs(peak_prominence))
    bkwd_peaks, _ = scipysig.find_peaks(np.abs(bkwd_cusum), height = abs(peak_height), prominence = abs(peak_prominence))
    # pair forward and backward peaks, going forward, if forward peak is not followed by a backward peak, ignore
    epoch_mask = np.full((nsamp,), False)
    fill_mask = np.zeros((nsamp,))
    for pk_idx in bkwd_peaks:
        fill_mask[pk_idx] = 1
    for pk_idx in fwd_peaks:
        fill_mask[pk_idx] = 2
    idx1 = 0    
    while idx1 < nsamp:
        # find a bkwd peak
        if fill_mask[idx1] != 1:
            idx1 += 1
            continue
        # find fwd peak, and if bkwd peak found again, make this the new peak
        idx2 = idx1+1
        while idx2 < nsamp:
            if fill_mask[idx2] == 1:
                break
            elif fill_mask[idx2] == 2:
                epoch_mask[idx1:idx2+1] = True
                break
            else:
                idx2 += 1
        idx1 = idx2         

    return epoch_mask, fwd_peaks, bkwd_peaks

# 08-23-20 note: old implementation with incorrect false positive rate
def detect_cusum_events_old(fl_sig, fs, par):
    """
    Transient detection using cummulative sum method for single and ratiometric channels.

    Parameters
    ----------
    fl_sig : dict
        Fluorescence signals, keys:
            'Ch1' : 1D numpy array, mandatory
                Responsive raw fluorescence signal. If it has bleaching, provide 'Ch1_db_fit'.
            'Ch1_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch1, otherwise it is added as the mean signal level.
            'Ch2' : 1D numpy array, optional
                If provided, this is a control channel used to construct a ratiometric signal and reduce common artifacts between Ch1 and Ch2, e.g. motion modulation.
                If it has bleaching, provide 'Ch2_db_fit'.
            'Ch2_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch2, otherwise it is added as the mean signal level.

    par : dict
        Signal processing and detection parameters, keys:

            'prefilt' : list of dict, optional
                Filters applied prior to ratiometric cusum. See specification for process.filter_bank
            'cusum' : dict, mandatory, keys:
                'shift' : float or 2 element iterable, mandatory
                    Expected change of signal w.r.t. baseline, e.g. 0.8 if peak response decreases from 1. If <1, a downward cusum is used, otherwise if >1, upward.
                'fp_rate' : float, mandatory
                    Maximum rate of false positive events in [Hz].
                'pk_prominence' : float, default 1
                    CUSUM peak detection prominence. 
                'fp_sim_factor' : float, default 50
                    How many times longer should be the scrambled cusum signal compared to the acceptable false positive duration.
            "ch2_siglevel_cutoff" : float, mandatory if second channel is used for ratiometric signal
                Channel 2 signal level cutoff after normalization to avoid divide by zero for low SNR when calculating ratio of Ch1/Ch2.

    Returns
    -------
    out : dict with keys 'up' or 'down' if shift >1 or <1 respectively, and dict values of the form:
    {
        'fwd_cusum' : 1D numpy array
            Forward direction cusum based on input signals.
        'bkwd_cusum' : 1D numpy array
            Backward direction cusum based on input signals.
        'ctrl_fwd_cusum' : 1D numpy array
            Forward direction cusum based on control (scrambled, no events) signals generated from input signals.
        'event_epochs' : 1D numpy bool
            Event segmented signal, with True marking event epochs.
        'cusum_fwd_peaks' : list
            Forward CUSUM peak array indices.
        'cusum_bkwd_peaks' : list
            Backward CUSUM peak array indices.
        'event_intervals' : list of tuples
            CUSUM event epoch start and end times as (t_start, t_end) tuples in [s].
            Note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval.   
        'cusum_threshold' : float
            CUSUM threshold chosen to meet the requested false positive rate 'fp_rate'.
    }
    """
    out = {}
    # set default parameters
    util.set_default_keys({'prefilt': [], 'cusum': {}}, par)
    util.set_default_keys({'pk_prominence': 1, 'fp_sim_factor': 50}, par['cusum'])

    if isinstance(par['cusum']['shift'], Iterable):
        cusum_shifts = par['cusum']['shift']
        # consistency checks
        if len(cusum_shifts) > 2:
            raise ValueError("Cannot do more than one up and down cusum at once. Results are stored per cusum type.")
        elif len(cusum_shifts) == 2 and (cusum_shifts[0] > 1 and cusum_shifts[1] > 1 or cusum_shifts[0] < 1 and cusum_shifts[1] < 1):
            raise ValueError("Cannot do more than one cusum type at the same time. Results are stored per cusum type.") 
    else:
        cusum_shifts = [par['cusum']['shift']]

    for cusum_shift in cusum_shifts:        
        if cusum_shift > 1:
            cusum_sign = 1
            cusum_direction = 'up'
        elif cusum_shift < 1:
            cusum_sign = -1
            cusum_direction = 'down'
        else:
            raise ValueError

        out[cusum_direction] = {}

        # normalize signals to bleaching fits and fill nans with 1 (assume no transients when signal is invalid)
        if 'Ch1_db_fit' not in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch1_db_fit'] = bn.nanmean(fl_sig['Ch1'][np.abs(iqr_score(fl_sig['Ch1'])) < 2])
        if 'Ch2_db_fit' not in fl_sig and 'Ch2' in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch2_db_fit'] = bn.nanmean(fl_sig['Ch2'][np.abs(iqr_score(fl_sig['Ch2'])) < 2])

        Ch1_norm = fl_sig['Ch1']/fl_sig['Ch1_db_fit']
        Ch1_norm[np.isnan(Ch1_norm)] = 1
        if 'Ch2' in fl_sig:
            Ch2_norm = fl_sig['Ch2']/fl_sig['Ch2_db_fit']
            Ch2_norm[np.isnan(Ch2_norm)] = 1
        
        # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
        Ch1_norm_filt = proc.filter_bank(sig = Ch1_norm, fs = fs, filters = par['prefilt'])
        # HP filters bring mean to 0, others may leave it unchanged, so bring back mean to 1
        Ch1_norm_filt = Ch1_norm_filt - bn.nanmean(Ch1_norm_filt[np.abs(iqr_score(Ch1_norm_filt)) < 2])+1
        # enforce positivity
        Ch1_norm_filt[Ch1_norm_filt<0] = 0
        Ch1_norm_filt[np.isnan(Ch1_norm_filt)] = 1 # since some filters may pad with nans
        if 'Ch2' in fl_sig:
            Ch2_norm_filt = proc.filter_bank(sig = Ch2_norm, fs = fs, filters = par['prefilt'])
            Ch2_norm_filt = Ch2_norm_filt - bn.nanmean(Ch2_norm_filt[np.abs(iqr_score(Ch2_norm_filt)) < 2])+1
            Ch2_norm_filt[Ch2_norm_filt<par['ch2_siglevel_cutoff']] = par['ch2_siglevel_cutoff']
            Ch2_norm_filt[np.isnan(Ch2_norm_filt)] = 1 # since some filters may pad with nans
        
        # apply here signal loss criteria for Ch1_norm_filt and Ch2_norm_filt (better after filter bank to be more confident)


        # generate response signal
        resp = Ch1_norm_filt*fl_sig['Ch1_db_fit'] # watch out, this can have nans from the bleach fit signal
        # generate target signal by transfering the (motion artifact) modulation to the response channel without with matching exponential bleaching
        if 'Ch2' in fl_sig:
            target = fl_sig['Ch1_db_fit'] * Ch2_norm_filt # watch out, this can have nans from the bleach fit signal
        else:
            target = fl_sig['Ch1_db_fit']

        # calculate cusums in both direction (signal peaks will be between the two cusum peaks)
        out[cusum_direction]['fwd_cusum'] = proc.cusum(x = resp, target = target, shift_factor = cusum_shift, direction = 'fwd') 
        out[cusum_direction]['bkwd_cusum'] = proc.cusum(x = resp, target = target, shift_factor = cusum_shift, direction = 'bkwd')
   
        # calculate below here a control signal that has noise characteristics as close to the original signals as possible
        # use ratiometric signal to reduce common motion artifacts


        if 'Ch2' in fl_sig:
            ratio = Ch1_norm_filt/Ch2_norm_filt
        else:
            ratio = Ch1_norm_filt
        
        # apply a 1 Hz HP filter to reduce slow drift variability
        # nb. Ch1_norm and Ch2_norm below here are further modified, if needed later, make a copy of these before

        ratio = 1+proc.filt_butter_hp(sig = ratio, fs = fs, fc = 1.0, order = 2)
        Ch1_norm = 1+proc.filt_butter_hp(sig = Ch1_norm, fs = fs, fc = 1.0, order = 2)
        if 'Ch2' in fl_sig:
            Ch2_norm = 1+proc.filt_butter_hp(sig = Ch2_norm, fs = fs, fc = 1.0, order = 2)

        # generate scrambled control signals using the noise properties of the original signals without large outliers such as spikes
        # and put back their bleaching characteristic

        # this should not happen
        if np.isnan(ratio).all():
            raise Exception

        Ch1_norm[np.abs(iqr_score(ratio)) > 2] = 1
        if 'Ch2' in fl_sig:
            Ch2_norm[np.abs(iqr_score(ratio)) > 2] = 1

        ctrl_cusum_peaks = [] # collect here control cusum peaks
        # adjust number of times to repeat simulation to generate at least fp_sim_factor false positive events/s
        n_ctrl_repeat = int(math.ceil((1./par['cusum']['fp_rate']*par['cusum']['fp_sim_factor'])/(len(ratio)/fs)))
        for rnd_seed in range(n_ctrl_repeat):
            # scramble both channels in a similar way to be able to take their ratio again and reduce movement artifacts, then put back bleaching
            np.random.seed(rnd_seed)
            np.random.shuffle(Ch1_norm)
            if 'Ch2' in fl_sig:
                np.random.seed(rnd_seed)
                np.random.shuffle(Ch2_norm)
      
            # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
            shuffled_Ch1_norm_filt = proc.filter_bank(sig = Ch1_norm, fs = fs, filters = par['prefilt'])
            shuffled_Ch1_norm_filt = shuffled_Ch1_norm_filt - bn.nanmean(shuffled_Ch1_norm_filt)+1
            shuffled_Ch1_norm_filt[np.isnan(shuffled_Ch1_norm_filt)] = 1 # since some filters may pad with nans
            if 'Ch2' in fl_sig:
                shuffled_Ch2_norm_filt = proc.filter_bank(sig = Ch2_norm, fs = fs, filters = par['prefilt'])
                shuffled_Ch2_norm_filt = shuffled_Ch2_norm_filt - bn.nanmean(shuffled_Ch2_norm_filt)+1
                shuffled_Ch2_norm_filt[np.isnan(shuffled_Ch2_norm_filt)] = 1 # since some filters may pad with nans

            # !!! apply somewhere in this section signal loss criteria for shuffled_Ch1_norm_filt and shuffled_Ch2_norm_filt

            # generate shuffled target and response signals
            shuffled_resp = shuffled_Ch1_norm_filt*fl_sig['Ch1_db_fit'] # watch out, this can have nans from the bleach fit signal
            # generate target signal by transfering the (motion artifact) modulation to the response channel without with matching exponential bleaching
            if 'Ch2' in fl_sig:
                shuffled_target = fl_sig['Ch1_db_fit'] * shuffled_Ch2_norm_filt # watch out, this can have nans from the bleach fit signal
            else:
                shuffled_target = fl_sig['Ch1_db_fit']

            # calculate control cusums in both direction (signal peaks will be between the two cusum peaks)
            # note this is done in one direction for now, will check actual fp rate.
            out[cusum_direction]['ctrl_fwd_cusum'] = proc.cusum(x = shuffled_resp, target = shuffled_target, shift_factor = cusum_shift, direction = 'fwd') 
            #out['ctrl_bkwd_cusum_'+cusum_direction] = proc.cusum(x = shuffled_Ch1_filt, target = shuffled_Ch2_filt, shift_factor = cusum_shift, direction = 'bkwd')
            ctrl_fwd_cusum_abs = np.abs(out[cusum_direction]['ctrl_fwd_cusum'])
            # collect ctrl cusum peaks
            ctrl_cusum_peaks.append(ctrl_fwd_cusum_abs[scipysig.find_peaks(ctrl_fwd_cusum_abs, prominence = par['cusum']['pk_prominence'])[0]])


        # reset seed generator
        np.random.seed()
        ctrl_cusum_peaks = np.sort(np.concatenate(ctrl_cusum_peaks))
        # get cusum threshold to satisfy false-positive rate
        if len(ctrl_cusum_peaks):
            out[cusum_direction]['cusum_threshold'] = cusum_sign*ctrl_cusum_peaks[max(int(min(-1, -len(ratio)*n_ctrl_repeat/fs*par['cusum']['fp_rate'])),
                                                        -len(ctrl_cusum_peaks))] 
        else:
            out[cusum_direction]['cusum_threshold'] = 0

        out[cusum_direction]['event_epochs'], out[cusum_direction]['cusum_fwd_peaks'], out[cusum_direction]['cusum_bkwd_peaks'] = \
            segment_cusum_peaks(out[cusum_direction]['fwd_cusum'], out[cusum_direction]['bkwd_cusum'],
                peak_height = out[cusum_direction]['cusum_threshold'], peak_prominence = par['cusum']['pk_prominence'])
        # list of cusum epoch times as (start,end) tuples in [s]
        # note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval
        out[cusum_direction]['event_intervals'] = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['event_epochs'], fs)

        return out

# old implementation as of 4/12/21
def detect_censored_cusum_events_old(fl_sig, fs, par):
    """
    Transient detection using cummulative sum method for given signal and control channels and uses control channel events to censor signal channel events.

    Parameters
    ----------
    fl_sig : dict
        Fluorescence signals, keys:
            'Ch1' : 1D numpy array, mandatory
                Responsive raw fluorescence signal. If it has bleaching, provide 'Ch1_db_fit'.
            'Ch1_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch1, otherwise it is added as the mean signal level.
            'Ch2' : 1D numpy array, optional
                If provided, this is a control channel used to censor signal channel events due to e.g. motion artifacts.
                If it has bleaching, provide 'Ch2_db_fit'.
            'Ch2_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch2, otherwise it is added as the mean signal level.
    
    fs : float
        Sampling frequency in [Hz].
        
    par : dict
        Signal processing and detection parameters, keys:
            'stab_filter' : required, dict
                Natural cubic spline filter to stabilize signal prior to cusum detection in the low frequency range. Dict with keys:
                'knot_interval' : float
                    Natural cubic spline knot time interval spacing in [ms], default 750 ms.
                'npass' : int
                    Number of times to exclude signal outliers and refit, default 2.
                'outlier_perc': float
                    Percentile of lower and upper outliers to exclude from fitting, default 5, i.e. data points < 5th percentile or > 95th percentile are excluded.

            'prefilt' : list of dict, optional
                Filters applied prior to ratiometric cusum. See specification for process.filter_bank
                Note: by default a natural cubic spline filter is applied before all filters to stabilize the signal.

            'cusum' : dict, mandatory, keys:
                'shift' : float or 2 element iterable, mandatory
                    Expected change of signal w.r.t. baseline, e.g. 0.8 if peak response decreases from 1. If <1, a downward cusum is used, otherwise if >1, upward.
                'fp_rate' : float, mandatory
                    Maximum rate of false positive events in [Hz].
                'pk_prominence' : float, default 1
                    CUSUM peak detection prominence. 
                'fp_sim_factor' : float, default 200
                    How many times longer should be the scrambled cusum signal compared to the acceptable false positive duration.
                'rand_block_duration' : float, default 1
                    Duration of randomization blocks in [s] used to divide the signal to estimate shot noise change due to bleaching. This value should be larger than
                    the expected duration of an event, but shorter than the bleaching time constant.

    Returns
    -------
    out : dict with keys 'up' or 'down' if shift >1 or <1 respectively, and dict values of the form:
    {
        'Ch1_fwd_cusum', 'Ch2_fwd_cusum' : 1D numpy array
            Forward direction cusum based on input signals.
        'Ch1_bkwd_cusum', 'Ch2_bkwd_cusum' : 1D numpy array
            Backward direction cusum based on input signals.
        'Ch1_ctrl_fwd_cusum', 'Ch2_ctrl_fwd_cusum' : 1D numpy array
            Forward direction cusum based on control (scrambled, no events) signals generated from input signals.
        'event_epochs' : 1D numpy bool
            Event segmented signal, with True marking event epochs.
        'cusum_fwd_peaks' : list
            Forward CUSUM peak array indices.
        'cusum_bkwd_peaks' : list
            Backward CUSUM peak array indices.
        'event_intervals' : list of tuples
            CUSUM event epoch start and end times as (t_start, t_end) tuples in [s].
            Note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval.   
        'cusum_threshold' : float
            CUSUM threshold chosen to meet the requested false positive rate 'fp_rate'.
        'cusum_valid_time' : float
            Total valid time used for event detection. Use this together with the number of events to calculate an event rate.
        'cusum_event_rate' : float

    }
    If cusum detection fails e.g. due to not being able to do debleaching, cusum detection direction returns empty dict, i.e. out['up'] = {} and out['down'] = {}
    """
    out = {}
    # set default parameters
    util.set_default_keys({'prefilt': [], 'cusum': {}, 'stab_filter': {}}, par)
    util.set_default_keys({'pk_prominence': 1, 'fp_sim_factor': 200, 'rand_block_duration': 1}, par['cusum'])
    util.set_default_keys({'knot_interval': 750, 'npass': 2, 'outlier_perc': 5}, par['stab_filter'])

    if isinstance(par['cusum']['shift'], Iterable):
        cusum_shifts = par['cusum']['shift']
        # consistency checks
        if len(cusum_shifts) > 2:
            raise ValueError("Cannot do more than one up and down cusum at once. Results are stored per cusum type.")
        elif len(cusum_shifts) == 2 and (cusum_shifts[0] > 1 and cusum_shifts[1] > 1 or cusum_shifts[0] < 1 and cusum_shifts[1] < 1):
            raise ValueError("Cannot do more than one cusum type at the same time. Results are stored per cusum type.") 
    else:
        cusum_shifts = [par['cusum']['shift']]

    for cusum_shift in cusum_shifts:        
        if cusum_shift > 1:
            cusum_sign = 1
            cusum_direction = 'up'
        elif cusum_shift < 1:
            cusum_sign = -1
            cusum_direction = 'down'
        else:
            raise ValueError

        out[cusum_direction] = {}

        # normalize signals to bleaching fits and fill nans with 1 (assume no transients when signal is invalid)
        # ======================================================================================================
        if 'Ch1_db_fit' not in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch1_db_fit'] = bn.nanmean(fl_sig['Ch1'][np.abs(iqr_score(fl_sig['Ch1'])) < 2])
        if 'Ch2_db_fit' not in fl_sig and 'Ch2' in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch2_db_fit'] = bn.nanmean(fl_sig['Ch2'][np.abs(iqr_score(fl_sig['Ch2'])) < 2])

        Ch1_norm = fl_sig['Ch1']/fl_sig['Ch1_db_fit']
        # mark nan samples
        invalid_Ch1_samples = np.isnan(Ch1_norm)
        if 'Ch2' in fl_sig:
            Ch2_norm = fl_sig['Ch2']/fl_sig['Ch2_db_fit']
            # mark nan samples
            invalid_Ch2_samples = np.isnan(Ch2_norm)

        # stabilize signals
        # ============================
        debleach_success = True
        try:
            Ch1_norm = proc.filt_nat_cubic_spline(Ch1_norm, fs, settings = par['stab_filter'])
        except:
            util.clrd_print('Cannot debleach Ch1 (response) cusum signal. {} percent are invalid samples.'.format(np.sum(np.isnan(Ch1_norm))/len(Ch1_norm)), 'warn')
            debleach_success = False

        if 'Ch2' in fl_sig:
            try:
                Ch2_norm = proc.filt_nat_cubic_spline(Ch2_norm, fs, settings = par['stab_filter'])
            except:
                util.clrd_print('Cannot debleach Ch2 (control) cusum signal. {} percent are invalid samples.'.format(np.sum(np.isnan(Ch2_norm))/len(Ch2_norm)), 'warn')
                debleach_success = False

        if not debleach_success:
            continue

        # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
        # ===================================================================================================
        Ch1_norm_filt = proc.filter_bank(sig = Ch1_norm, fs = fs, filters = par['prefilt'])
        Ch1_norm_filt = Ch1_norm_filt - bn.nanmean(Ch1_norm_filt[np.abs(iqr_score(Ch1_norm_filt)) < 2])+1
        
        if 'Ch2' in fl_sig:
            Ch2_norm_filt = proc.filter_bank(sig = Ch2_norm, fs = fs, filters = par['prefilt'])
            Ch2_norm_filt = Ch2_norm_filt - bn.nanmean(Ch2_norm_filt[np.abs(iqr_score(Ch2_norm_filt)) < 2])+1
            
        # calculate cusums in both direction for each channel separately (signal peaks will be between the two cusum peaks)
        out[cusum_direction]['Ch1_fwd_cusum'] = proc.cusum(x = Ch1_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd') 
        out[cusum_direction]['Ch1_bkwd_cusum'] = proc.cusum(x = Ch1_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
        if 'Ch2' in fl_sig:
            out[cusum_direction]['Ch2_fwd_cusum'] = proc.cusum(x = Ch2_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd') 
            out[cusum_direction]['Ch2_bkwd_cusum'] = proc.cusum(x = Ch2_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')

        # determine cusum thresholds to obtain target false positive rates
        # by calculating a randomized control signal that has noise characteristics as
        # close as possible to the original signals
        # ================================================================
        # reuse original normalized signals and take out extreme values (this will slightly increase SNR but it does more good than harm)
        
        # WARNING: make sure nans are passed for iqr_score calculation to be correct
        Ch1_norm[np.abs(iqr_score(Ch1_norm)) > 3.5] = 1
        if 'Ch2' in fl_sig:
            Ch2_norm[np.abs(iqr_score(Ch2_norm)) > 3.5] = 1

        # integer number of times to generate randomized signals to have at least par['cusum']['fp_sim_factor'] false positive events on average,
        # as integer multiple of input signal duration, rounded up
        sig_dur = len(Ch1_norm)/fs
        n_ctrl_iter = int(math.ceil((1./par['cusum']['fp_rate']*par['cusum']['fp_sim_factor'])/sig_dur))
        nctrl_samples = int(sig_dur*fs)
        rand_block_nsamp = int(par['cusum']['rand_block_duration']*fs)

        rand_blocks_Ch1 = [Ch1_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch1_norm)/rand_block_nsamp))+1)]
        # remove last empty slice if blocks fit exactly
        if not len(rand_blocks_Ch1[-1]):
            del rand_blocks_Ch1[-1]
        if 'Ch2' in fl_sig:
            rand_blocks_Ch2 = [Ch2_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch2_norm)/rand_block_nsamp))+1)]
            # remove last empty slice if blocks fit exactly
            if not len(rand_blocks_Ch2[-1]):
                del rand_blocks_Ch2[-1]

        Ch1_ctrl_fwd_cusums = []
        Ch2_ctrl_fwd_cusums = []
        Ch1_ctrl_bkwd_cusums = []
        Ch2_ctrl_bkwd_cusums = []
        Ch1_max_ctrl_cusum = 0
        Ch2_max_ctrl_cusum = 0

        for rand_sig_idx in range(n_ctrl_iter):
            for block_idx in range(len(rand_blocks_Ch1)):
                # shuffle blocks in place, fix seed to apply same shuffle to both channels
                rnd_block_seed = rand_sig_idx*len(rand_blocks_Ch1)+block_idx
                np.random.seed(rnd_block_seed)
                np.random.shuffle(rand_blocks_Ch1[block_idx])
                if 'Ch2' in fl_sig:
                    np.random.seed(rnd_block_seed)
                    np.random.shuffle(rand_blocks_Ch2[block_idx])

            # join blocks and apply filters as they were applied to the original signals
            Ch1_norm_rand = np.concatenate(rand_blocks_Ch1)
            # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
            Ch1_norm_rand_filt = proc.filter_bank(sig = Ch1_norm_rand, fs = fs, filters = par['prefilt'])
            Ch1_norm_rand_filt = Ch1_norm_rand_filt - bn.nanmean(Ch1_norm_rand_filt)+1
            # collect forward and backward control cusums
            Ch1_ctrl_fwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
            Ch1_ctrl_bkwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
            Ch1_ctrl_fwd_cusums.append(Ch1_ctrl_fwd_cusum)
            Ch1_ctrl_bkwd_cusums.append(Ch1_ctrl_bkwd_cusum)
            Ch1_max_ctrl_cusum = max(max(np.abs(Ch1_ctrl_fwd_cusum)), max(np.abs(Ch1_ctrl_bkwd_cusum)), Ch1_max_ctrl_cusum)
            if 'Ch2' in fl_sig:
                Ch2_norm_rand = np.concatenate(rand_blocks_Ch2)
                Ch2_norm_rand_filt = proc.filter_bank(sig = Ch2_norm_rand, fs = fs, filters = par['prefilt'])
                Ch2_norm_rand_filt = Ch2_norm_rand_filt - bn.nanmean(Ch2_norm_rand_filt)+1
                # collect forward and backward control cusums
                Ch2_ctrl_fwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
                Ch2_ctrl_bkwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
                Ch2_ctrl_fwd_cusums.append(Ch2_ctrl_fwd_cusum)
                Ch2_ctrl_bkwd_cusums.append(Ch2_ctrl_bkwd_cusum)
                Ch2_max_ctrl_cusum = max(max(np.abs(Ch2_ctrl_fwd_cusum)), max(np.abs(Ch2_ctrl_bkwd_cusum)), Ch2_max_ctrl_cusum)
        
        # reset random seed
        np.random.seed()

        # determine Ch1 cusum threshold
        # ====================================
        # build threshold adjustment array
        if Ch1_max_ctrl_cusum > par['cusum']['pk_prominence']:
            Ch1_thresholds = [0, Ch1_max_ctrl_cusum/2, Ch1_max_ctrl_cusum]
            # do 7x threshold iterations
            for thresh_idx in range(7):
                # detect events for control signal
                n_fp_events = 0
                for rand_sig_idx in range(n_ctrl_iter):
                    ctrl_event_epochs, _, _ = \
                    segment_cusum_peaks(Ch1_ctrl_fwd_cusums[rand_sig_idx], Ch1_ctrl_bkwd_cusums[rand_sig_idx],
                        peak_height = Ch1_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
                    # calculate false positive rate
                    ctrl_event_intervals = proc.convert_bool_epochs_to_intervals(ctrl_event_epochs, fs)
                    n_fp_events += len(ctrl_event_intervals)   

                actual_fp_rate = n_fp_events/(n_ctrl_iter*(sig_dur-np.sum(invalid_Ch1_samples)/fs))
                if actual_fp_rate < par['cusum']['fp_rate']:
                    Ch1_thresholds = [Ch1_thresholds[0],(Ch1_thresholds[0]+Ch1_thresholds[1])/2.,Ch1_thresholds[1]]
                elif actual_fp_rate > par['cusum']['fp_rate']:
                    Ch1_thresholds = [Ch1_thresholds[1],(Ch1_thresholds[1]+Ch1_thresholds[2])/2.,Ch1_thresholds[2]]
                else:
                    break
                
            out[cusum_direction]['Ch1_cusum_threshold'] = Ch1_thresholds[1]

            out[cusum_direction]['event_epochs'], out[cusum_direction]['Ch1_cusum_fwd_peaks'], out[cusum_direction]['Ch1_cusum_bkwd_peaks'] = \
            segment_cusum_peaks(out[cusum_direction]['Ch1_fwd_cusum'], out[cusum_direction]['Ch1_bkwd_cusum'],
                peak_height = Ch1_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
            
        else:
            out[cusum_direction]['Ch1_cusum_threshold'] = None
            out[cusum_direction]['event_epochs'] = np.full((len(out[cusum_direction]['Ch1_fwd_cusum']),), False)
            out[cusum_direction]['Ch1_cusum_fwd_peaks'] = np.array([])
            out[cusum_direction]['Ch1_cusum_bkwd_peaks'] = np.array([])

        Ch1_event_intervals = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['event_epochs'], fs)


        # determine Ch2 cusum threshold
        # ====================================
        # build threshold adjustment array
        if 'Ch2' in fl_sig:
            if Ch2_max_ctrl_cusum > par['cusum']['pk_prominence']:
                Ch2_thresholds = [0, Ch2_max_ctrl_cusum/2, Ch2_max_ctrl_cusum]
                # do 5x threshold iterations
                for thresh_idx in range(7):
                    # detect events for control signal
                    n_fp_events = 0
                    for rand_sig_idx in range(n_ctrl_iter):
                        ctrl_event_epochs, _, _ = \
                        segment_cusum_peaks(Ch2_ctrl_fwd_cusums[rand_sig_idx], Ch2_ctrl_bkwd_cusums[rand_sig_idx],
                            peak_height = Ch2_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
                        # calculate false positive rate
                        ctrl_event_intervals = proc.convert_bool_epochs_to_intervals(ctrl_event_epochs, fs)
                        n_fp_events += len(ctrl_event_intervals)   

                    actual_fp_rate = n_fp_events/(n_ctrl_iter*(sig_dur-np.sum(invalid_Ch2_samples)/fs))
                    if actual_fp_rate < par['cusum']['fp_rate']:
                        Ch2_thresholds = [Ch2_thresholds[0],(Ch2_thresholds[0]+Ch2_thresholds[1])/2.,Ch2_thresholds[1]]
                    elif actual_fp_rate > par['cusum']['fp_rate']:
                        Ch2_thresholds = [Ch2_thresholds[1],(Ch2_thresholds[1]+Ch2_thresholds[2])/2.,Ch2_thresholds[2]]
                    else:
                        break
                out[cusum_direction]['Ch2_cusum_threshold'] = Ch2_thresholds[1]
                out[cusum_direction]['Ch2_event_epochs'], out[cusum_direction]['Ch2_cusum_fwd_peaks'], out[cusum_direction]['Ch2_cusum_bkwd_peaks'] = \
                segment_cusum_peaks(out[cusum_direction]['Ch2_fwd_cusum'], out[cusum_direction]['Ch2_bkwd_cusum'],
                    peak_height = Ch2_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
            else:
                out[cusum_direction]['Ch2_cusum_threshold'] = None
                out[cusum_direction]['Ch2_event_epochs'] = np.full((len(out[cusum_direction]['Ch2_fwd_cusum']),), False)
                out[cusum_direction]['Ch2_cusum_fwd_peaks'] = np.array([])
                out[cusum_direction]['Ch2_cusum_bkwd_peaks'] = np.array([])

            Ch2_event_intervals = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['Ch2_event_epochs'], fs)

        
        if 'Ch2' in fl_sig:      
            # censor Ch1 event intervals that overlap with Ch2
            censor_idx = []
            for idx1, ch1_interval in enumerate(Ch1_event_intervals):
                for ch2_interval in Ch2_event_intervals:
                    if ch2_interval[0] > ch1_interval[1]:
                        break
                    elif ch2_interval[1] > ch1_interval[0]:
                        censor_idx.append(idx1)   
            # censor Ch1 events
            util.del_list_idxs(Ch1_event_intervals, censor_idx)
            out[cusum_direction]['event_epochs'] = proc.intervals_to_bool(Ch1_event_intervals, fs, len(out[cusum_direction]['event_epochs']))

        # list of cusum epoch times as (start,end) tuples in [s]
        # note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval
        out[cusum_direction]['event_intervals'] = Ch1_event_intervals 
        
        # calculate valid cusum detection time
        if 'Ch2' in fl_sig:
            invalid_samples = np.logical_or(invalid_Ch1_samples, invalid_Ch2_samples)
        else:
            invalid_samples = invalid_Ch1_samples

        out[cusum_direction]['cusum_valid_time'] = (len(Ch1_norm)-np.sum(invalid_samples))/fs
        out[cusum_direction]['cusum_event_rate'] = len(out[cusum_direction]['event_intervals'])/out[cusum_direction]['cusum_valid_time']

    return out

# new implementation as of 4/12/21, adding noise model for a better event threshold estimation
def detect_censored_cusum_events(fl_sig, fs, par):
    """
    Transient detection using cummulative sum method for given signal and control channels and uses control channel events to censor signal channel events.

    Parameters
    ----------
    fl_sig : dict
        Fluorescence signals, keys:
            'Ch1' : 1D numpy array, mandatory
                Responsive raw fluorescence signal. If it has bleaching, provide 'Ch1_db_fit'.
            'Ch1_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch1, otherwise it is added as the mean signal level.
            'Ch2' : 1D numpy array, optional
                If provided, this is a control channel used to censor signal channel events due to e.g. motion artifacts.
                If it has bleaching, provide 'Ch2_db_fit'.
            'Ch2_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch2, otherwise it is added as the mean signal level.
    
    fs : float
        Sampling frequency in [Hz].
        
    par : dict
        Signal processing and detection parameters, keys:
            'stab_filter' : required, dict
                Natural cubic spline filter to stabilize signal prior to cusum detection in the low frequency range. Dict with keys:
                'knot_interval' : float
                    Natural cubic spline knot time interval spacing in [ms], default 750 ms.
                'npass' : int
                    Number of times to exclude signal outliers and refit, default 2.
                'outlier_perc': float
                    Percentile of lower and upper outliers to exclude from fitting, default 5, i.e. data points < 5th percentile or > 95th percentile are excluded.

            'prefilt' : list of dict, optional
                Filters applied prior to ratiometric cusum. See specification for process.filter_bank
                Note: by default a natural cubic spline filter is applied before all filters to stabilize the signal.

            'cusum' : dict, mandatory, keys:
                'shift' : float or 2 element iterable, mandatory
                    Expected change of signal w.r.t. baseline, e.g. 0.8 if peak response decreases from 1. If <1, a downward cusum is used, otherwise if >1, upward.
                'fp_rate' : float, mandatory
                    Maximum rate of false positive events in [Hz].
                'pk_prominence' : float, default 1
                    CUSUM peak detection prominence. 
                'fp_sim_factor' : float, default 200
                    How many times longer should be the scrambled cusum signal compared to the acceptable false positive duration.
                'rand_block_duration' : float, default 1
                    Duration of randomization blocks in [s] used to divide the signal to estimate shot noise change due to bleaching. This value should be larger than
                    the expected duration of an event, but shorter than the bleaching time constant.
                'noise_model' : str, default 'auto'
                    Signal noise model for adjusting event detection threshold. Choose between 'auto' and 'poisson'. Use 'auto' if noise level is estimated from the
                    the given signals or use 'poisson' if signals are calibrated photon counts and poisson statistics should be used.

    Returns
    -------
    out : dict with keys 'up' or 'down' if shift >1 or <1 respectively, and dict values of the form:
    {
        'Ch1_fwd_cusum', 'Ch2_fwd_cusum' : 1D numpy array
            Forward direction cusum based on input signals.
        'Ch1_bkwd_cusum', 'Ch2_bkwd_cusum' : 1D numpy array
            Backward direction cusum based on input signals.
        'Ch1_ctrl_fwd_cusum', 'Ch2_ctrl_fwd_cusum' : 1D numpy array
            Forward direction cusum based on control (scrambled, no events) signals generated from input signals.
        'event_epochs' : 1D numpy bool
            Event segmented signal, with True marking event epochs.
        'cusum_fwd_peaks' : list
            Forward CUSUM peak array indices.
        'cusum_bkwd_peaks' : list
            Backward CUSUM peak array indices.
        'event_intervals' : list of tuples
            CUSUM event epoch start and end times as (t_start, t_end) tuples in [s].
            Note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval.   
        'cusum_threshold' : float
            CUSUM threshold chosen to meet the requested false positive rate 'fp_rate'.
        'cusum_valid_time' : float
            Total valid time used for event detection. Use this together with the number of events to calculate an event rate.
        'cusum_event_rate' : float

    }
    If cusum detection fails e.g. due to not being able to do debleaching, cusum detection direction returns empty dict, i.e. out['up'] = {} and out['down'] = {}
    """
    out = {}
    # set default parameters
    util.set_default_keys({'prefilt': [], 'cusum': {}, 'stab_filter': {}}, par)
    util.set_default_keys({'pk_prominence': 1, 'fp_sim_factor': 200, 'rand_block_duration': 1, 'noise_model': 'auto'}, par['cusum'])
    util.set_default_keys({'knot_interval': 750, 'npass': 2, 'outlier_perc': 5}, par['stab_filter'])

    if isinstance(par['cusum']['shift'], Iterable):
        cusum_shifts = par['cusum']['shift']
        # consistency checks
        if len(cusum_shifts) > 2:
            raise ValueError("Cannot do more than one up and down cusum at once. Results are stored per cusum type.")
        elif len(cusum_shifts) == 2 and (cusum_shifts[0] > 1 and cusum_shifts[1] > 1 or cusum_shifts[0] < 1 and cusum_shifts[1] < 1):
            raise ValueError("Cannot do more than one cusum type at the same time. Results are stored per cusum type.") 
    else:
        cusum_shifts = [par['cusum']['shift']]

    for cusum_shift in cusum_shifts:        
        if cusum_shift > 1:
            cusum_sign = 1
            cusum_direction = 'up'
        elif cusum_shift < 1:
            cusum_sign = -1
            cusum_direction = 'down'
        else:
            raise ValueError

        out[cusum_direction] = {}

        # normalize signals to bleaching fits and fill nans with 1 (assume no transients when signal is invalid)
        # ======================================================================================================
        if 'Ch1_db_fit' not in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch1_db_fit'] = np.full((len(fl_sig['Ch1']),), bn.nanmean(fl_sig['Ch1'][np.abs(iqr_score(fl_sig['Ch1'])) < 2]))
        if 'Ch2_db_fit' not in fl_sig and 'Ch2' in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch2_db_fit'] = np.full((len(fl_sig['Ch2']),), bn.nanmean(fl_sig['Ch2'][np.abs(iqr_score(fl_sig['Ch2'])) < 2]))

        Ch1_norm = fl_sig['Ch1']/fl_sig['Ch1_db_fit']
        # mark nan samples
        invalid_Ch1_samples = np.isnan(Ch1_norm)
        if 'Ch2' in fl_sig:
            Ch2_norm = fl_sig['Ch2']/fl_sig['Ch2_db_fit']
            # mark nan samples
            invalid_Ch2_samples = np.isnan(Ch2_norm)

        # stabilize signals
        # ============================
        debleach_success = True
        try:
            Ch1_norm = proc.filt_nat_cubic_spline(Ch1_norm, fs, settings = par['stab_filter'])
        except:
            util.clrd_print('Cannot debleach Ch1 (response) cusum signal. {} percent are invalid samples.'.format(np.sum(np.isnan(Ch1_norm))/len(Ch1_norm)), 'warn')
            debleach_success = False

        if 'Ch2' in fl_sig:
            try:
                Ch2_norm = proc.filt_nat_cubic_spline(Ch2_norm, fs, settings = par['stab_filter'])
            except:
                util.clrd_print('Cannot debleach Ch2 (control) cusum signal. {} percent are invalid samples.'.format(np.sum(np.isnan(Ch2_norm))/len(Ch2_norm)), 'warn')
                debleach_success = False

        if not debleach_success:
            continue

        # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
        # ===================================================================================================
        Ch1_norm_filt = proc.filter_bank(sig = Ch1_norm, fs = fs, filters = par['prefilt'])
        Ch1_norm_filt = Ch1_norm_filt - bn.nanmean(Ch1_norm_filt[np.abs(iqr_score(Ch1_norm_filt)) < 2])+1
        
        if 'Ch2' in fl_sig:
            Ch2_norm_filt = proc.filter_bank(sig = Ch2_norm, fs = fs, filters = par['prefilt'])
            Ch2_norm_filt = Ch2_norm_filt - bn.nanmean(Ch2_norm_filt[np.abs(iqr_score(Ch2_norm_filt)) < 2])+1
            
        # calculate cusums in both direction for each channel separately (signal peaks will be between the two cusum peaks)
        out[cusum_direction]['Ch1_fwd_cusum'] = proc.cusum(x = Ch1_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd') 
        out[cusum_direction]['Ch1_bkwd_cusum'] = proc.cusum(x = Ch1_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
        if 'Ch2' in fl_sig:
            out[cusum_direction]['Ch2_fwd_cusum'] = proc.cusum(x = Ch2_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd') 
            out[cusum_direction]['Ch2_bkwd_cusum'] = proc.cusum(x = Ch2_norm_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')

        # determine cusum thresholds to obtain target false positive rates
        # by either calculating a randomized control signal that has noise characteristics as
        # close as possible to the original signals or use poisson shot noise level if signal is calibrated as photon counts
        # ================================================================

        Ch1_ctrl_fwd_cusums = []
        Ch2_ctrl_fwd_cusums = []
        Ch1_ctrl_bkwd_cusums = []
        Ch2_ctrl_bkwd_cusums = []
        Ch1_max_ctrl_cusum = 0
        Ch2_max_ctrl_cusum = 0

        # integer number of times to generate randomized signals to have at least par['cusum']['fp_sim_factor'] false positive events on average,
        # as integer multiple of input signal duration, rounded up
        sig_dur = len(Ch1_norm)/fs
        n_ctrl_iter = int(math.ceil((1./par['cusum']['fp_rate']*par['cusum']['fp_sim_factor'])/sig_dur))
        nctrl_samples = int(sig_dur*fs)

        if par['cusum']['noise_model'] == 'auto':

            # reuse original normalized signals and take out extreme values (this will slightly increase SNR but it does more good than harm)
            # WARNING: make sure nans are passed for iqr_score calculation to be correct
            Ch1_norm[np.abs(iqr_score(Ch1_norm)) > 3.5] = 1
            if 'Ch2' in fl_sig:
                Ch2_norm[np.abs(iqr_score(Ch2_norm)) > 3.5] = 1

            
            rand_block_nsamp = int(par['cusum']['rand_block_duration']*fs)

            rand_blocks_Ch1 = [Ch1_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch1_norm)/rand_block_nsamp))+1)]
            # remove last empty slice if blocks fit exactly
            if not len(rand_blocks_Ch1[-1]):
                del rand_blocks_Ch1[-1]
            if 'Ch2' in fl_sig:
                rand_blocks_Ch2 = [Ch2_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch2_norm)/rand_block_nsamp))+1)]
                # remove last empty slice if blocks fit exactly
                if not len(rand_blocks_Ch2[-1]):
                    del rand_blocks_Ch2[-1]

            for rand_sig_idx in range(n_ctrl_iter):
                for block_idx in range(len(rand_blocks_Ch1)):
                    # shuffle blocks in place, fix seed to apply same shuffle to both channels
                    rnd_block_seed = rand_sig_idx*len(rand_blocks_Ch1)+block_idx
                    np.random.seed(rnd_block_seed)
                    np.random.shuffle(rand_blocks_Ch1[block_idx])
                    if 'Ch2' in fl_sig:
                        np.random.seed(rnd_block_seed)
                        np.random.shuffle(rand_blocks_Ch2[block_idx])

                # join blocks and apply filters as they were applied to the original signals
                Ch1_norm_rand = np.concatenate(rand_blocks_Ch1)
                if 'Ch2' in fl_sig:
                    Ch2_norm_rand = np.concatenate(rand_blocks_Ch2)

                # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
                Ch1_norm_rand_filt = proc.filter_bank(sig = Ch1_norm_rand, fs = fs, filters = par['prefilt'])
                Ch1_norm_rand_filt = Ch1_norm_rand_filt - bn.nanmean(Ch1_norm_rand_filt)+1
                # collect forward and backward control cusums
                Ch1_ctrl_fwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
                Ch1_ctrl_bkwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
                Ch1_ctrl_fwd_cusums.append(Ch1_ctrl_fwd_cusum)
                Ch1_ctrl_bkwd_cusums.append(Ch1_ctrl_bkwd_cusum)
                Ch1_max_ctrl_cusum = max(max(np.abs(Ch1_ctrl_fwd_cusum)), max(np.abs(Ch1_ctrl_bkwd_cusum)), Ch1_max_ctrl_cusum)
                if 'Ch2' in fl_sig:
                    Ch2_norm_rand_filt = proc.filter_bank(sig = Ch2_norm_rand, fs = fs, filters = par['prefilt'])
                    Ch2_norm_rand_filt = Ch2_norm_rand_filt - bn.nanmean(Ch2_norm_rand_filt)+1
                    # collect forward and backward control cusums
                    Ch2_ctrl_fwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
                    Ch2_ctrl_bkwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
                    Ch2_ctrl_fwd_cusums.append(Ch2_ctrl_fwd_cusum)
                    Ch2_ctrl_bkwd_cusums.append(Ch2_ctrl_bkwd_cusum)
                    Ch2_max_ctrl_cusum = max(max(np.abs(Ch2_ctrl_fwd_cusum)), max(np.abs(Ch2_ctrl_bkwd_cusum)), Ch2_max_ctrl_cusum)
        
            # reset random seed
            np.random.seed()

        elif par['cusum']['noise_model'] == 'poisson':
            for rand_sig_idx in range(n_ctrl_iter):
                Ch1_norm_rand = np.random.poisson(fl_sig['Ch1_db_fit'])/fl_sig['Ch1_db_fit']
                if 'Ch2' in fl_sig:
                    Ch2_norm_rand = np.random.poisson(fl_sig['Ch2_db_fit'])/fl_sig['Ch2_db_fit']

                # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
                Ch1_norm_rand_filt = proc.filter_bank(sig = Ch1_norm_rand, fs = fs, filters = par['prefilt'])
                Ch1_norm_rand_filt = Ch1_norm_rand_filt - bn.nanmean(Ch1_norm_rand_filt)+1
                # collect forward and backward control cusums
                Ch1_ctrl_fwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
                Ch1_ctrl_bkwd_cusum = proc.cusum(x = Ch1_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
                Ch1_ctrl_fwd_cusums.append(Ch1_ctrl_fwd_cusum)
                Ch1_ctrl_bkwd_cusums.append(Ch1_ctrl_bkwd_cusum)
                Ch1_max_ctrl_cusum = max(max(np.abs(Ch1_ctrl_fwd_cusum)), max(np.abs(Ch1_ctrl_bkwd_cusum)), Ch1_max_ctrl_cusum)
                if 'Ch2' in fl_sig:
                    Ch2_norm_rand_filt = proc.filter_bank(sig = Ch2_norm_rand, fs = fs, filters = par['prefilt'])
                    Ch2_norm_rand_filt = Ch2_norm_rand_filt - bn.nanmean(Ch2_norm_rand_filt)+1
                    # collect forward and backward control cusums
                    Ch2_ctrl_fwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'fwd')
                    Ch2_ctrl_bkwd_cusum = proc.cusum(x = Ch2_norm_rand_filt, target = 1, shift_factor = cusum_shift, direction = 'bkwd')
                    Ch2_ctrl_fwd_cusums.append(Ch2_ctrl_fwd_cusum)
                    Ch2_ctrl_bkwd_cusums.append(Ch2_ctrl_bkwd_cusum)
                    Ch2_max_ctrl_cusum = max(max(np.abs(Ch2_ctrl_fwd_cusum)), max(np.abs(Ch2_ctrl_bkwd_cusum)), Ch2_max_ctrl_cusum)
        else:
            raise ValueError("cusum parameter 'noise_model' must be 'auto' or 'poisson'.")

        # determine Ch1 cusum threshold
        # ====================================
        # build threshold adjustment array
        if Ch1_max_ctrl_cusum > par['cusum']['pk_prominence']:
            Ch1_thresholds = [0, Ch1_max_ctrl_cusum/2, Ch1_max_ctrl_cusum]
            # do 7x threshold iterations
            for thresh_idx in range(7):
                # detect events for control signal
                n_fp_events = 0
                for rand_sig_idx in range(n_ctrl_iter):
                    ctrl_event_epochs, _, _ = \
                    segment_cusum_peaks(Ch1_ctrl_fwd_cusums[rand_sig_idx], Ch1_ctrl_bkwd_cusums[rand_sig_idx],
                        peak_height = Ch1_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
                    # calculate false positive rate
                    ctrl_event_intervals = proc.convert_bool_epochs_to_intervals(ctrl_event_epochs, fs)
                    n_fp_events += len(ctrl_event_intervals)   

                actual_fp_rate = n_fp_events/(n_ctrl_iter*(sig_dur-np.sum(invalid_Ch1_samples)/fs))
                if actual_fp_rate < par['cusum']['fp_rate']:
                    Ch1_thresholds = [Ch1_thresholds[0],(Ch1_thresholds[0]+Ch1_thresholds[1])/2.,Ch1_thresholds[1]]
                elif actual_fp_rate > par['cusum']['fp_rate']:
                    Ch1_thresholds = [Ch1_thresholds[1],(Ch1_thresholds[1]+Ch1_thresholds[2])/2.,Ch1_thresholds[2]]
                else:
                    break
                
            out[cusum_direction]['Ch1_cusum_threshold'] = Ch1_thresholds[1]

            out[cusum_direction]['event_epochs'], out[cusum_direction]['Ch1_cusum_fwd_peaks'], out[cusum_direction]['Ch1_cusum_bkwd_peaks'] = \
            segment_cusum_peaks(out[cusum_direction]['Ch1_fwd_cusum'], out[cusum_direction]['Ch1_bkwd_cusum'],
                peak_height = Ch1_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
            
        else:
            out[cusum_direction]['Ch1_cusum_threshold'] = None
            out[cusum_direction]['event_epochs'] = np.full((len(out[cusum_direction]['Ch1_fwd_cusum']),), False)
            out[cusum_direction]['Ch1_cusum_fwd_peaks'] = np.array([])
            out[cusum_direction]['Ch1_cusum_bkwd_peaks'] = np.array([])

        Ch1_event_intervals = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['event_epochs'], fs)


        # determine Ch2 cusum threshold
        # ====================================
        # build threshold adjustment array
        if 'Ch2' in fl_sig:
            if Ch2_max_ctrl_cusum > par['cusum']['pk_prominence']:
                Ch2_thresholds = [0, Ch2_max_ctrl_cusum/2, Ch2_max_ctrl_cusum]
                # do 5x threshold iterations
                for thresh_idx in range(7):
                    # detect events for control signal
                    n_fp_events = 0
                    for rand_sig_idx in range(n_ctrl_iter):
                        ctrl_event_epochs, _, _ = \
                        segment_cusum_peaks(Ch2_ctrl_fwd_cusums[rand_sig_idx], Ch2_ctrl_bkwd_cusums[rand_sig_idx],
                            peak_height = Ch2_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
                        # calculate false positive rate
                        ctrl_event_intervals = proc.convert_bool_epochs_to_intervals(ctrl_event_epochs, fs)
                        n_fp_events += len(ctrl_event_intervals)   

                    actual_fp_rate = n_fp_events/(n_ctrl_iter*(sig_dur-np.sum(invalid_Ch2_samples)/fs))
                    if actual_fp_rate < par['cusum']['fp_rate']:
                        Ch2_thresholds = [Ch2_thresholds[0],(Ch2_thresholds[0]+Ch2_thresholds[1])/2.,Ch2_thresholds[1]]
                    elif actual_fp_rate > par['cusum']['fp_rate']:
                        Ch2_thresholds = [Ch2_thresholds[1],(Ch2_thresholds[1]+Ch2_thresholds[2])/2.,Ch2_thresholds[2]]
                    else:
                        break
                out[cusum_direction]['Ch2_cusum_threshold'] = Ch2_thresholds[1]
                out[cusum_direction]['Ch2_event_epochs'], out[cusum_direction]['Ch2_cusum_fwd_peaks'], out[cusum_direction]['Ch2_cusum_bkwd_peaks'] = \
                segment_cusum_peaks(out[cusum_direction]['Ch2_fwd_cusum'], out[cusum_direction]['Ch2_bkwd_cusum'],
                    peak_height = Ch2_thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
            else:
                out[cusum_direction]['Ch2_cusum_threshold'] = None
                out[cusum_direction]['Ch2_event_epochs'] = np.full((len(out[cusum_direction]['Ch2_fwd_cusum']),), False)
                out[cusum_direction]['Ch2_cusum_fwd_peaks'] = np.array([])
                out[cusum_direction]['Ch2_cusum_bkwd_peaks'] = np.array([])

            Ch2_event_intervals = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['Ch2_event_epochs'], fs)

        
        if 'Ch2' in fl_sig:      
            # censor Ch1 event intervals that overlap with Ch2
            censor_idx = []
            for idx1, ch1_interval in enumerate(Ch1_event_intervals):
                for ch2_interval in Ch2_event_intervals:
                    if ch2_interval[0] > ch1_interval[1]:
                        break
                    elif ch2_interval[1] > ch1_interval[0]:
                        censor_idx.append(idx1)   
            # censor Ch1 events
            util.del_list_idxs(Ch1_event_intervals, censor_idx)
            out[cusum_direction]['event_epochs'] = proc.intervals_to_bool(Ch1_event_intervals, fs, len(out[cusum_direction]['event_epochs']))

        # list of cusum epoch times as (start,end) tuples in [s]
        # note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval
        out[cusum_direction]['event_intervals'] = Ch1_event_intervals 
        
        # calculate valid cusum detection time
        if 'Ch2' in fl_sig:
            invalid_samples = np.logical_or(invalid_Ch1_samples, invalid_Ch2_samples)
        else:
            invalid_samples = invalid_Ch1_samples

        out[cusum_direction]['cusum_valid_time'] = (len(Ch1_norm)-np.sum(invalid_samples))/fs
        out[cusum_direction]['cusum_event_rate'] = len(out[cusum_direction]['event_intervals'])/out[cusum_direction]['cusum_valid_time']

    return out

# 08-23-20 note: new implementation with correct false positive rate
def detect_ratio_cusum_events(fl_sig, fs, par):
    """
    Transient detection using cummulative sum method for single and ratiometric channels.

    Parameters
    ----------
    fl_sig : dict
        Fluorescence signals, keys:
            'Ch1' : 1D numpy array, mandatory
                Responsive raw fluorescence signal. If it has bleaching, provide 'Ch1_db_fit'.
            'Ch1_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch1, otherwise it is added as the mean signal level.
            'Ch2' : 1D numpy array, optional
                If provided, this is a control channel used to construct a ratiometric signal and reduce common artifacts between Ch1 and Ch2, e.g. motion modulation.
                If it has bleaching, provide 'Ch2_db_fit'.
            'Ch2_db_fit' : 1D numpy array, optional
                If provided, fluorescence bleaching fit that normalizes Ch2, otherwise it is added as the mean signal level.

    fs : float
        Sampling frequency in [Hz].

    par : dict
        Signal processing and detection parameters, keys:

            'prefilt' : list of dict, optional
                Filters applied prior to ratiometric cusum. See specification for process.filter_bank
            'cusum' : dict, mandatory, keys:
                'shift' : float or 2 element iterable, mandatory
                    Expected change of signal w.r.t. baseline, e.g. 0.8 if peak response decreases from 1. If <1, a downward cusum is used, otherwise if >1, upward.
                'fp_rate' : float, mandatory
                    Maximum rate of false positive events in [Hz].
                'pk_prominence' : float, default 1
                    CUSUM peak detection prominence. 
                'fp_sim_factor' : float, default 200
                    How many times longer should be the scrambled cusum signal compared to the acceptable false positive duration.
                'rand_block_duration': float, default 1
                    Duration of randomization blocks in [s] used to divide the signal to estimate shot noise change due to bleaching. This value should be larger than
                    the expected duration of an event, but shorter than the bleaching time constant.
            "ch2_siglevel_cutoff" : float, mandatory if second channel is used for ratiometric signal
                Channel 2 signal level cutoff after normalization to avoid divide by zero for low SNR when calculating ratio of Ch1/Ch2.

    Returns
    -------
    out : dict with keys 'up' or 'down' if shift >1 or <1 respectively, and dict values of the form:
    {
        'fwd_cusum' : 1D numpy array
            Forward direction cusum based on input signals.
        'bkwd_cusum' : 1D numpy array
            Backward direction cusum based on input signals.
        'ctrl_fwd_cusum' : 1D numpy array
            Forward direction cusum based on control (scrambled, no events) signals generated from input signals.
        'event_epochs' : 1D numpy bool
            Event segmented signal, with True marking event epochs.
        'cusum_fwd_peaks' : list
            Forward CUSUM peak array indices.
        'cusum_bkwd_peaks' : list
            Backward CUSUM peak array indices.
        'event_intervals' : list of tuples
            CUSUM event epoch start and end times as (t_start, t_end) tuples in [s].
            Note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval.   
        'cusum_threshold' : float
            CUSUM threshold chosen to meet the requested false positive rate 'fp_rate'.
        'cusum_valid_time' : float
            Total valid time used for event detection. Use this together with the number of events to calculate an event rate.
        'cusum_event_rate' : float

    }
    """
    out = {}
    # set default parameters
    util.set_default_keys({'prefilt': [], 'cusum': {}}, par)
    util.set_default_keys({'pk_prominence': 1, 'fp_sim_factor': 200, 'rand_block_duration': 1}, par['cusum'])

    if isinstance(par['cusum']['shift'], Iterable):
        cusum_shifts = par['cusum']['shift']
        # consistency checks
        if len(cusum_shifts) > 2:
            raise ValueError("Cannot do more than one up and down cusum at once. Results are stored per cusum type.")
        elif len(cusum_shifts) == 2 and (cusum_shifts[0] > 1 and cusum_shifts[1] > 1 or cusum_shifts[0] < 1 and cusum_shifts[1] < 1):
            raise ValueError("Cannot do more than one cusum type at the same time. Results are stored per cusum type.") 
    else:
        cusum_shifts = [par['cusum']['shift']]

    for cusum_shift in cusum_shifts:        
        if cusum_shift > 1:
            cusum_sign = 1
            cusum_direction = 'up'
        elif cusum_shift < 1:
            cusum_sign = -1
            cusum_direction = 'down'
        else:
            raise ValueError

        out[cusum_direction] = {}

        # normalize signals to bleaching fits and fill nans with 1 (assume no transients when signal is invalid)
        if 'Ch1_db_fit' not in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch1_db_fit'] = bn.nanmean(fl_sig['Ch1'][np.abs(iqr_score(fl_sig['Ch1'])) < 2])
        if 'Ch2_db_fit' not in fl_sig and 'Ch2' in fl_sig:
            # calculate mean signal level by leaving out sporadic outliers
            fl_sig['Ch2_db_fit'] = bn.nanmean(fl_sig['Ch2'][np.abs(iqr_score(fl_sig['Ch2'])) < 2])

        Ch1_norm = fl_sig['Ch1']/fl_sig['Ch1_db_fit']
        # mark nan samples
        nan_Ch1_norm = np.isnan(Ch1_norm)
        Ch1_norm[nan_Ch1_norm] = 1
        # keep track of invalid samples for cusum detection and event rate calculation
        invalid_samples = np.copy(nan_Ch1_norm)
        if 'Ch2' in fl_sig:
            Ch2_norm = fl_sig['Ch2']/fl_sig['Ch2_db_fit']
            # mark nan samples
            nan_Ch2_norm = np.isnan(Ch2_norm)
            Ch2_norm[nan_Ch2_norm] = 1
            invalid_samples = np.logical_or(invalid_samples, nan_Ch2_norm)
        
        # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
        Ch1_norm_filt = proc.filter_bank(sig = Ch1_norm, fs = fs, filters = par['prefilt'])
        # HP filters bring mean to 0, others may leave it unchanged, so bring back mean to 1
        Ch1_norm_filt = Ch1_norm_filt - bn.nanmean(Ch1_norm_filt[np.abs(iqr_score(Ch1_norm_filt)) < 2])+1
        # enforce positivity
        positivity_Ch1_norm_filt = Ch1_norm_filt<0
        Ch1_norm_filt[positivity_Ch1_norm_filt] = 0
        nan_Ch1_norm_filt = np.isnan(Ch1_norm_filt)
        Ch1_norm_filt[nan_Ch1_norm_filt] = 1 # since some filters may pad with nans
        invalid_samples = np.logical_or(np.logical_or(invalid_samples, nan_Ch1_norm_filt), positivity_Ch1_norm_filt)
        if 'Ch2' in fl_sig:
            Ch2_norm_filt = proc.filter_bank(sig = Ch2_norm, fs = fs, filters = par['prefilt'])
            Ch2_norm_filt = Ch2_norm_filt - bn.nanmean(Ch2_norm_filt[np.abs(iqr_score(Ch2_norm_filt)) < 2])+1
            sig_level_Ch2_norm_filt = Ch2_norm_filt<par['ch2_siglevel_cutoff']
            Ch2_norm_filt[sig_level_Ch2_norm_filt] = par['ch2_siglevel_cutoff']
            nan_Ch2_norm_filt = np.isnan(Ch2_norm_filt)
            Ch2_norm_filt[nan_Ch2_norm_filt] = 1 # since some filters may pad with nans
            invalid_samples = np.logical_or(np.logical_or(invalid_samples, nan_Ch2_norm_filt), sig_level_Ch2_norm_filt)

        # generate response signal
        resp = Ch1_norm_filt*fl_sig['Ch1_db_fit'] # watch out, this can have nans from the bleach fit signal
        # generate target signal by transfering the (motion artifact) modulation to the response channel without with matching exponential bleaching
        if 'Ch2' in fl_sig:
            target = fl_sig['Ch1_db_fit'] * Ch2_norm_filt # watch out, this can have nans from the bleach fit signal
        else:
            target = fl_sig['Ch1_db_fit']

        # calculate cusums in both direction (signal peaks will be between the two cusum peaks)
        out[cusum_direction]['fwd_cusum'] = proc.cusum(x = resp, target = target, shift_factor = cusum_shift, direction = 'fwd') 
        out[cusum_direction]['bkwd_cusum'] = proc.cusum(x = resp, target = target, shift_factor = cusum_shift, direction = 'bkwd')
   
        # calculate below here a control signal that has noise characteristics as close to the original signals as possible
        # use ratiometric signal to reduce common motion artifacts

        if 'Ch2' in fl_sig:
            ratio = Ch1_norm_filt/Ch2_norm_filt
        else:
            ratio = Ch1_norm_filt

        Ch1_norm[np.abs(iqr_score(ratio)) > 3.5] = 1
        if 'Ch2' in fl_sig:
            Ch2_norm[np.abs(iqr_score(ratio)) > 3.5] = 1

        # generate randomized control signals
        # ======================================================

        # integer number of times to generate randomized signals to have at least par['cusum']['fp_sim_factor'] false positive events on average,
        # as integer multiple of input signal duration, rounded up
        sig_dur = len(resp)/fs
        n_ctrl_iter = int(math.ceil((1./par['cusum']['fp_rate']*par['cusum']['fp_sim_factor'])/sig_dur))
        nctrl_samples = int(sig_dur*fs)
        rand_block_nsamp = int(par['cusum']['rand_block_duration']*fs)

        rand_blocks_Ch1 = [Ch1_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch1_norm)/rand_block_nsamp))+1)]
        # remove last empty slice if blocks fit exactly
        if not len(rand_blocks_Ch1[-1]):
            del rand_blocks_Ch1[-1]
        if 'Ch2' in fl_sig:
            rand_blocks_Ch2 = [Ch2_norm[block_idx*rand_block_nsamp:(block_idx+1)*rand_block_nsamp] for block_idx in range(int(math.ceil(len(Ch2_norm)/rand_block_nsamp))+1)]
            # remove last empty slice if blocks fit exactly
            if not len(rand_blocks_Ch2[-1]):
                del rand_blocks_Ch2[-1]

        ctrl_fwd_cusums = []
        ctrl_bkwd_cusums = []
        max_ctrl_cusum = 0

        for rand_sig_idx in range(n_ctrl_iter):
            for block_idx in range(len(rand_blocks_Ch1)):
                # shuffle blocks in place, fix seed to apply same shuffle to both channels
                rnd_block_seed = rand_sig_idx*len(rand_blocks_Ch1)+block_idx
                np.random.seed(rnd_block_seed)
                np.random.shuffle(rand_blocks_Ch1[block_idx])
                if 'Ch2' in fl_sig:
                    np.random.seed(rnd_block_seed)
                    np.random.shuffle(rand_blocks_Ch2[block_idx])

            # join blocks and apply filters as they were applied to the original signals
            Ch1_norm_rand = np.concatenate(rand_blocks_Ch1)
            # apply filter bank to each signal and bring back mean signal level to 1 in case HP filters were used
            Ch1_norm_rand_filt = proc.filter_bank(sig = Ch1_norm_rand, fs = fs, filters = par['prefilt'])
            Ch1_norm_rand_filt = Ch1_norm_rand_filt - bn.nanmean(Ch1_norm_rand_filt)+1
            Ch1_norm_rand_filt[np.isnan(Ch1_norm_rand_filt)] = 1 # since some filters may pad with nans
            if 'Ch2' in fl_sig:
                Ch2_norm_rand = np.concatenate(rand_blocks_Ch2)
                Ch2_norm_rand_filt = proc.filter_bank(sig = Ch2_norm_rand, fs = fs, filters = par['prefilt'])
                Ch2_norm_rand_filt = Ch2_norm_rand_filt - bn.nanmean(Ch2_norm_rand_filt)+1
                Ch2_norm_rand_filt[np.isnan(Ch2_norm_rand_filt)] = 1 # since some filters may pad with nans
        
            # generate randomized target and response signals with bleaching characteristics
            rand_resp = Ch1_norm_rand_filt*fl_sig['Ch1_db_fit'] # watch out, this can have nans from the bleach fit signal
            # generate target signal by transfering the (motion artifact) modulation to the response channel without with matching exponential bleaching
            if 'Ch2' in fl_sig:
                rand_target = Ch2_norm_rand_filt*fl_sig['Ch1_db_fit'] # watch out, this can have nans from the bleach fit signal
            else:
                rand_target = fl_sig['Ch1_db_fit']

            # collect forward and backward control cusums
            ctrl_fwd_cusum = proc.cusum(x = rand_resp, target = rand_target, shift_factor = cusum_shift, direction = 'fwd')
            ctrl_bkwd_cusum = proc.cusum(x = rand_resp, target = rand_target, shift_factor = cusum_shift, direction = 'bkwd')
            ctrl_fwd_cusums.append(ctrl_fwd_cusum)
            ctrl_bkwd_cusums.append(ctrl_bkwd_cusum)
            max_ctrl_cusum = max(max(np.abs(ctrl_fwd_cusum)), max(np.abs(ctrl_bkwd_cusum)), max_ctrl_cusum)

        # reset random seed
        np.random.seed()

        # build threshold adjustment array
        thresholds = [0,max_ctrl_cusum/2,max_ctrl_cusum]
        # do 5x threshold iterations
        for thresh_idx in range(5):
            # detect events for control signal
            n_fp_events = 0
            for rand_sig_idx in range(n_ctrl_iter):
                ctrl_event_epochs, _, _ = \
                segment_cusum_peaks(ctrl_fwd_cusums[rand_sig_idx], ctrl_bkwd_cusums[rand_sig_idx],
                    peak_height = thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
                # calculate false positive rate
                ctrl_event_intervals = proc.convert_bool_epochs_to_intervals(ctrl_event_epochs, fs)
                n_fp_events += len(ctrl_event_intervals)   

            actual_fp_rate = n_fp_events/(n_ctrl_iter*sig_dur)
            if actual_fp_rate < par['cusum']['fp_rate']:
                thresholds = [thresholds[0],(thresholds[0]+thresholds[1])/2.,thresholds[1]]
            elif actual_fp_rate > par['cusum']['fp_rate']:
                thresholds = [thresholds[1],(thresholds[1]+thresholds[2])/2.,thresholds[2]]
            else:
                break

        out[cusum_direction]['event_epochs'], out[cusum_direction]['cusum_fwd_peaks'], out[cusum_direction]['cusum_bkwd_peaks'] = \
            segment_cusum_peaks(out[cusum_direction]['fwd_cusum'], out[cusum_direction]['bkwd_cusum'],
                peak_height = thresholds[1], peak_prominence = par['cusum']['pk_prominence'])
        # list of cusum epoch times as (start,end) tuples in [s]
        # note: this is not the actual event duration, rather an event segmentation window with an event peak within this interval
        out[cusum_direction]['event_intervals'] = proc.convert_bool_epochs_to_intervals(out[cusum_direction]['event_epochs'], fs)
        out[cusum_direction]['cusum_threshold'] = thresholds[1]
        # calculate valid cusum detection time
        out[cusum_direction]['cusum_valid_time'] = (len(Ch1_norm)-np.sum(invalid_samples))/fs
        out[cusum_direction]['cusum_event_rate'] = len(out[cusum_direction]['event_intervals'])/out[cusum_direction]['cusum_valid_time']

    return out

def detect_poisson_single_exp_pulse(fl, baselineSNR, ml_thresh = None, dF_over_F = 0.05, fs = 20, tau = 0.15):
    """
    Detects single exponential fluorescence pulse events. Implementation follows:
    Wilt, B. A., Fitzgerald, J. E., & Schnitzer, M. J. (2013). Photon shot noise limits on optical detection of
    neuronal spikes and estimation of spike timing. Biophysical journal, 104(1), 51-62.

    Parameters
    ----------
    fl : 1D numpy array
        Detected fluorescence with each sample being the number of photons per integration time bin 1/fs.
    baselineSNR : float
        Baseline signal-to-noise ratio of the measured optical signal at fs sampling rate, i.e. with a 1/fs integration window.
    ml_thresh : float, None
        Maximum-likelihood event detection threshold.
    dF_over_F : float
        Expected fluorescence change, e.g. 0.2 for 20% increase over F0 baline.
    fs : float
        Sampling rate in [Hz].
    tau : float
        Single exponential pulse decay time-constant in [s].

    Returns
    -------   
    tuple (1D numpy nd.array float, 1D numpy nd.array bool or None)
        Tuple of maximum likelihood and detected events with given threshold.
    """
    F0 = baselineSNR**2*fs
    # n - 1-index time bin index
    inverse_fs = 1/fs
    inverse_tau_fs = 1/(tau*fs)
    const1 = dF_over_F*tau*(math.exp(inverse_tau_fs)-1)
    s_bar = np.vectorize(lambda n: F0*(inverse_fs + const1*math.exp(-n*inverse_tau_fs)))
    pulse = np.vectorize(lambda n: F0*const1*math.exp(-n*inverse_tau_fs))
    b_bar = F0/fs
    # number of samples to reach cutoff 1% amplitude pulse cutoff
    pulse_cutoff = 1e-2
    ns = int(math.ceil(tau*math.log(1./pulse_cutoff)*fs))
    # precalculate fixed terms
    s_bar_val = s_bar(np.arange(ns))
    pulse_val = pulse(np.arange(ns))
    t1 = np.log(s_bar_val/b_bar)
    t2 = np.sum(s_bar_val-b_bar)
    # likelihood ratio vector
    ml = np.empty((len(fl)-ns,))
    for idx in range(len(fl)-ns):
        ml[idx] = np.inner(fl[idx:idx+ns], t1) - t2
    # detect events by iteratively setting the fluorescence vector to 0 at the maximum of ml as long as ml_thresh is crossed
    if ml_thresh is not None:
        fl_copy = np.copy(fl)
        ml_copy = np.copy(ml) 
        events = np.full((len(fl),), False)
        while True:
            valid_idx = np.where(ml_copy > ml_thresh)[0] 
            if not len(valid_idx):
                break
            ml_peaks = np.unravel_index(valid_idx[ml_copy[valid_idx].argmax()], ml_copy.shape)
            for ml_pk in ml_peaks:
                events[ml_pk] = True
                fl_copy[ml_pk:ml_pk+len(pulse_val)] = fl_copy[ml_pk:ml_pk+len(pulse_val)]-pulse_val
                fl_copy[fl_copy<0] = 0 # ensure it is strict positive after peeling e.g. false positive

            for idx in range(len(fl_copy)-ns):
                ml_copy[idx] = np.inner(fl_copy[idx:idx+ns], t1) - t2
    else:
        events = None

    return ml, events

def generate_single_exp_pulse_train(baselineSNR, fs, pulse_rate, dF_over_F, tau, npulses):
    """
    Generates a series of regularly spaced, non-overlaping single-exponential pulses
    under Poisson shot noise conditions of the form:
        S(t) = F0+A*exp(-t/tau)*theta(t)
    where:
        F0 - time-independent background fluorescence rate
        A - pulse amplitude
        tau - decay time-constant
        theta(t) - step function, =0 if argument is negative and 1 otherwise.

    Another way to rewrite the above signal:
        S(t) = F0*(1+dF_over_F*exp(-t/tau)*theta(t))

    The connection between the baseline SNR and F0 is:
        baselineSNR^2 = F0/fs

    Parameters
    ----------
    baselineSNR : float
        Baseline signal-to-noise ratio of the measured optical signal at fs sampling rate, i.e. with a 1/fs integration window.
    fs : float
        Sampling rate in [Hz].
    pulse_rate: float
        Rate of pulses to be detected in [Hz].
    dF_over_F : float
        Normalized peak-amplitude of fluorescence events.
    tau : float
        Single-exponential decay time-constant in [s].
    npulses : int
        How many pulses to generate.

    Returns
    -------
    1D numpy.ndarray
        Evenly spaced pulse train under Poisson noise. First and last pulses are at 1/pulse_rate from beginning and end of array.
    """
    F0 = baselineSNR**2*fs
    inverse_tau_fs = 1/(tau*fs)
    const1 = dF_over_F*tau*(math.exp(inverse_tau_fs)-1)
    s_bar = np.vectorize(lambda n: F0*(1/fs + const1*math.exp(-n*inverse_tau_fs)))
    b_bar = F0/fs
    # number of samples to reach 1% pulse amplitude cutoff
    pulse_cutoff = 1e-2
    # number of samples to describe a pulse
    ns_pulse = int(math.ceil(tau*math.log(1./pulse_cutoff)*fs))
    # total number of samples within a pulse period
    ns_period = int(fs/pulse_rate)
    # ensure pulse is not longer than the pulse period
    assert ns_period >= ns_pulse
    # number of baseline samples to fill
    ns_baseline = ns_period-ns_pulse

    sig_mean = np.concatenate([np.full((ns_period,), b_bar), np.tile(np.concatenate([s_bar(np.arange(ns_pulse)),np.full((ns_baseline,), b_bar)]), npulses)])
    return np.random.poisson(sig_mean)
    
def calc_behavior_states_summary(df):
    """
    Aggregates behavior states information from multiple sweeps to calculate a percentage of time for each cell ID
    during which the animal is still, running, or other intermediate state.

    Parameters
    ----------
    df : pandas DataFrame
        Behavior state entries from multiple recordings for each cell that need to be summarized. The dataframe must have the
        following columns:
            'cells_cellID' : str
                Unique cell ID assigned from the DB.
            'still_dur' : float
                Still duration in [s] for a recording.
            'loc_dur' : float
                Locomotion duration in [s] for a recording.
            'total_dur' : float
                Total recording duration in [s].

    Returns
    -------
    pandas DataFrame
        Summary of behavior state per cell ID, DataFrame with cellID index and columns:
            'still_dur_perc' : float
                Percentage still duration to total duration of recordings.
            'loc_dur_perc' : float
                Percentage locomotion duration to toal duration of recordings.
            'other_dur_perc' : float
                Percentage other behavior states duration that are neither still nor locomotion.

        Note that 'still_dur_perc'+'loc_dur_perc'+'other_dur_perc' = 1
    """
    # calculate sum of times per cell, reindex with cellID
    df_summed = df.groupby(['cells_cellID'], as_index = False).sum()
    df_summed.set_index('cells_cellID', inplace = True)
    # normalize locomotion and still by total time
    out_df = pd.DataFrame(index = df_summed.index)
    out_df[['still_dur_perc', 'loc_dur_perc']] = df_summed[['still_dur', 'loc_dur']].div(df_summed['total_dur'], axis = 0)*100
    # for completeness add other behavior states such that normalized sum of still, locomotion and other adds up to 1
    out_df['other_dur_perc'] = 100 - out_df['still_dur_perc'] - out_df['loc_dur_perc']

    return out_df

def bootstrap_mean_waveform(waveforms, ci = 95, niter = 10000):
    """
    Parameters
    ----------
    waveforms : 2D numpy.ndarray
        Waveforms to use for mean and CI calculation. Trial waveforms are organized by rows.
    ci : float
        Confidence interval, default 95% CI.
    niter : int
        Number of bootstrap iterations.
    
    Returns
    -------
    tuple of 1D numpy.ndarray
        Mean, lower and upper CI waveforms as tuple (mean, lower_ci, upper_ci).
    """
    # ignore runtime warnings for slices containing only nans
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category = RuntimeWarning)
        mean = bn.nanmean(waveforms, axis = 0)
        bs = sb_bootstrap(waveforms, axis = 0, niter = niter, func = bn.nanmean)
        lower_ci, upper_ci = nan_ci(a = bs, which = ci, axis = 0)

    return mean, lower_ci, upper_ci

def nan_ci(a, which = 95, axis = None):
    """
    Return a percentile range from an array of values ignoring nans.
    """
    p = 50 - which / 2, 50 + which / 2
    return np.nanpercentile(a, p, axis)

def gaussian_pulse(fwhm_dur, window_dur, fs):
    """
    Generates a gaussian pulse woth peak amplitude normalized to 1.

    Parameters
    ----------
    fwhm_dur : float
        Gaussian pulse full width at half maximum (FWHM) in [s].
    window_dur : float
        Pulse window duration in [s]; this should be at least 2*fwhm_dur to avoid clipping too much.
    fs : float
        Sampling rate in [Hz].
    """
    std_nsamp = fwhm_dur/(2*(2*np.log(2))**0.5)*fs
    ns_window = util.round_up_to_odd_int(fs*window_dur)
    return scipysig.gaussian(ns_window, std_nsamp, sym = True)

def poisson_gaussian_pulse_train(baselineSNR, fs, pulse_rate, dF_over_F, fwhm_dur, npulses):
    """
    Generates a series of regularly spaced Gaussian pulses under poisson shot noise limit.

    Parameters
    ----------
    baselineSNR : float
        Baseline signal-to-noise ratio of the measured optical signal at fs sampling rate, i.e. with a 1/fs integration window.
    fs : float
        Sampling rate in [Hz].
    pulse_rate: float
        Rate of pulses to be detected in [Hz].
    dF_over_F : float
        Normalized peak-amplitude of fluorescence events.
    fwhm_dur : float
        Gaussian pulse full width at half maximum (FWHM) in [s].
    """
    return np.random.poisson(np.tile(baselineSNR**2*(1+dF_over_F*gaussian_pulse(fwhm_dur, 1./pulse_rate, fs)), npulses))

def get_event_interval_time_centroid(sig, event_intervals, fs):
    """
    Determines pulse-shaped events time centroid.

    Parameters
    ----------
    sig : 1D numpy array
        Signal containing event pulses. Signal must be normalized to the event-free baseline.
    event_intervals : iterable of 2 elements
        Event intervals in [s] as t_start and t_end pairs.
    fs : float
        Sampling rate in [Hz].

    Returns
    -------
    list
        Estimated time points in [s].
    """
    evt_interval_idxs = [(int(round(evt_interval[0]*fs)), int(round(evt_interval[1]*fs))) for \
                evt_interval in event_intervals]

    if any([idx[0] < 0 for idx in evt_interval_idxs]):
        raise Exception("Event startpoint is outside of signal bounds.")

    if any([idx[1] >= len(sig) for idx in evt_interval_idxs]):
        raise Exception("Event endpoint is outside of signal bounds.")

    out = [1./fs*np.average(np.arange(evt_idx[0],evt_idx[1]), weights = np.abs((1-sig[evt_idx[0]:evt_idx[1]]))) \
        for evt_idx in evt_interval_idxs]
    sig_duration = len(sig)/fs
    return out

def exp_rise_fall_pulse(t, tau_rise, tau_decay, amp = 1):
    """
    Pulse with single exponential rise and decay time.

    Parameters
    ----------
    t : 1D numpy.array
        Pulse sampling times in [s], evenly spaced.
    tau_rise, tau_decay : float
        Pulse rise and decay time constants in [ms].
    amp : float  
        Pulse amplitude. Pulse rises and falls to a 0 baseline.
    
    Returns
    -------
    numpy.ndarray
        Pulse rise starts at t=0 s.
    """
    out = np.piecewise(t, [t > 0, t <= 0], [lambda t: np.exp(-t/tau_decay*1e3)-np.exp(-t/tau_rise*1e3), 0])
    out /= np.max(out)
    return out*amp

def jittered_exp_rise_fall_pulse(t, tau_rise, tau_decay, amp = 1, n = 100):
    """
    Jittered sampling of pulses with single exponential rise and decay time within the uncertainty of sampling interval.
  
    Parameters
    ----------
    t : numpy.array
        Sampling times in [s].
    tau_rise, tau_decay : float
        Pulse rise and decay time constants in [ms].
    amp : float  
        Pulse amplitude. Pulse rises and falls to a 0 baseline.
    n : int
        Number of jittered pulses.

    Returns
    -------
    np.ndarray of shape (n,len(t))
        Jittered pulses, row-wise.
    """
    # calculate sampling rate from first two samples
    fs = 1/(t[1]-t[0])
    out = np.empty((n, len(t)))
    for i in range(n):
        out[i,:] = exp_rise_fall_pulse(t = t+np.random.uniform(-0.5/fs, 0.5/fs), tau_rise = tau_rise,
                         tau_decay = tau_decay, amp = amp)
    return out

def precision_recall(detected_times, true_times, max_dt):
    """
    Calculate precission and recall given ground truth event times and detected event times.

    Parameters
    ----------
    detected_times, true_times : 1D numpy array
        Detected and ground truth event times. Must have same time units.
    max_dt : float
        Maximum time difference between detected and ground truth times to consider as detected event. Same time units as detected_times and true_times.

    Returns
    -------
    tuple of floats
        (precision, recall)
    """

    cost_matrix = np.abs(detected_times.reshape(1,-1) - true_times.reshape(-1,1))
    true_idx, detected_idx = opt.linear_sum_assignment(cost_matrix)
    
    distances = cost_matrix[true_idx, detected_idx]
    true_idx = true_idx[distances < max_dt]
    detected_idx = detected_idx[distances < max_dt]
    
    precision = float(len(true_idx)) / len(detected_times) # positive predictive value
    recall = float(len(true_idx)) / len(true_times)
    
    return precision, recall

def sampled_exp_rise_fall_pulse(fs, tau_d, tau_r = None, amp = 1):
    """
    Generates a single pulse with decaying and rising (optionally) exponential time constants that
    is sampled by integration at sampling rate fs. Pulse baseline is 0.

    Parameters
    ----------
    fs : float
        Sampling rate in [Hz].
    tau_d,tau_r : float
        Pulse rise and decay time constants in [s].
    amp : float
        Pulse amplitude measured from 0 baseline.

    Returns
    -------
    1D numpy.array
        Pulse waveform.
    """
    inverse_tau_d_fs = 1/(tau_d*fs)
    const_d = fs*tau_d*(math.exp(inverse_tau_d_fs)-1)

    pulse_cutoff = 1e-2
    if tau_r is None:
        s_bar = np.vectorize(lambda n: const_d*math.exp(-n*inverse_tau_d_fs))
        # number of samples to reach pulse amplitude cutoff
        ns_pulse = int(math.ceil(tau_d*math.log(1./pulse_cutoff)*fs))
    else:
        # find peak for normalization
        def f(t, tr, td): return math.exp(-t/td)-math.exp(-t/tr)
        max_t = scipyopt.fmin(lambda t: -f(t, tau_r, tau_d), 0, disp = False)[0]
        max_amp = f(max_t, tau_r, tau_d)

        inverse_tau_r_fs = 1/(tau_r*fs)
        const_r = fs*tau_r*(math.exp(inverse_tau_r_fs)-1)/max_amp
        const_d /= max_amp
        s_bar = np.vectorize(lambda n: const_d*math.exp(-n*inverse_tau_d_fs) - const_r*math.exp(-n*inverse_tau_r_fs))
        # number of samples to reach pulse amplitude cutoff
        ns_pulse = int(math.ceil((tau_r+tau_d)*math.log(1./pulse_cutoff)*fs))

    return amp*s_bar(np.arange(1,1+ns_pulse))

def detect_mf_events(fl, fl_baseline, fl_fs, mf_template, pars, template_fs = None):
    """
    Event detection in noisy fluorescence signal under poisson noise using a matched-filter.

    Parameters
    ----------
    fl : 1D numpy.array
        Normalized fluorescence signal with unit mean (no bleaching).
    fl_baseline : float, 1D numpy.array
        Fluorescence baseline measured as photon counts per sample.
    fl_fs : float
        Sampling rate in [Hz].
    mf_template : 1D numpy.array
        Matched filter template. Template will be time-reversed and normalized to area under the curve.
        Start of event time-point within the template must be centered and number of samples must be odd.
    pars : dict
        Matched filter settings. Dict with keys and values:
            'event_direction' : mandatory, str
                Event direction. Choose 'rising' or 'falling'.
            'hp_cut' : optional, float
                If specified, apply an effective 4-th order zero-phase butter filter with given cut-off frequency in [Hz],
                after applying the matched filter and before thresholded event detection.
            'min_event_interval' : mandatory, float
                Minimum event interval used for detection in [ms]. Smaller events are removed first until the condition
                is fulfilled for all remaining events.
            'fp_rate' : mandatory, float
                False positive rate in [Hz].
            'min_threshold' : optional, float
                Minimum threshold to use for event detection. If event detection threshold for given fp_rate is below
                this value, it is adjusted to it and the fp_rate is recalculated.
            'fl_centroid_threshold' : optional, float, default 0
                Apply a fluorescence level threshold for centroid calculation of event time and centroid fluorescence amplitude.
                Note: threshold is always >= 0 regardless of event direction (fluorescence threshold relative to baseline).
            'isolated_event_interval' : mandatory, float
                Event is considered isolated if before or after it within this interval in [ms] there are no other events.

    template_fs : None or float
        Template sampling rate if different than fluorescence sampling rate.

    Returns
    -------
    out : dict
        Matched filter event detection, dict with keys and values:
            'events' : pandas.DataFrame
                Detected events, dataframe with columns:
                    'mf-peak-time' : float
                        Match-filter signal detected event peak time in [s].
                    'mf-peak-amplitude' : float
                        Match-filter signal detected event peak amplitude.
                    'photon-rate' : float
                        Number of measured photons for event baseline at time of event.
                    'is_isolated' : bool
                        If True, event is isolated from other events according to 'isolated_event_interval' parameter.
                    'fl-peak-time' : float
                        Fluorescence peak time in [s] measured within a window of 2*min_event_interval centered around the match-filter detected event.
                    'fl-peak-amplitude' : float
                        Fluorescence peak amplitude measured at 'fl-peak-time'.
                    'fl-centroid-time' : float
                        Fluorescence centroid within a window of 2*min_event_interval centered around the match-filter detected event for samples that
                        exceed 'fl_centroid_threshold'. If threshold is too high and centroid cannot be calculated, this will be set to NaN.
                    'fl-centroid-amplitude' : float
                        Fluorescence centroid amplitude measured at 'fl-centroid-time'. If threshold is too high and centroid cannot be calculated, this will be set to NaN.
            'threshold' : float
                Given or adjusted event detection threshold for the match-filtered signal.
            'template' : 1D numpy.array
                Given or resampled (sampling rate = fl_fs) normalized matched-filter template (non-time reversed).
            'valid_time': float
                Valid event detection time within the given signal in [s] when considering nan samples in 'fl'.
            'mf_sig' : 1D numpy.array
                Match-filtered signal. Inherits nans from 'fl'.
    """
    if not len(mf_template)%2:
        raise Exception("Template must have an odd number of samples.")
    
    # resample template if template sampling rate is different from the fluorescence signal
    if template_fs is not None and template_fs != fl_fs:
        rfactor = fl_fs/template_fs
        mf_template = scipysig.resample(mf_template, util.round_to_nearest_odd_int(len(mf_template)*rfactor))
        # recenter template
        mf_template = proc.padded_shift(arr = mf_template, n = int((len(mf_template)-1)/2)-np.argmax(mf_template))
    
    # normalize templates to area under template and reverse to prepare for convolution
    norm_mf_template = mf_template/np.sum(mf_template)
    rev_mf_template = np.flipud(norm_mf_template)

    # adjust fl_baseline to always have an array of same length as fl
    if np.isscalar(fl_baseline):
        fl_baseline = np.full((len(fl),), fl_baseline)
    
    # bring baseline to 0 to avoid convolution artifacts and ensure events are positive going to match the filter direction
    if pars['event_direction'] == 'rising':
        fl_sig = fl-1
    elif pars['event_direction'] == 'falling':
        fl_sig = 1-fl
    else:
        raise ValueError("Matched filter setting 'event_direction' must be either 'rising' or 'falling'.")

    # apply matched filter
    mf_sig = np.convolve(fl_sig, rev_mf_template, 'same') # nans are passed through
    mf_sig = mf_sig[:len(fl_sig)]

    # apply hp filter
    if 'hp_cut' in pars:
        mf_sig = proc.filt_butter_hp(sig = mf_sig, fs = fl_fs, fc = pars['hp_cut'], order = 2) # nans will pass through

    # calculate threshold for event detection using a poisson noise model that follows the mean photon counts of the bleach fit curve
    pks = np.array([]) # contains all peaks
    niter = 0
    while len(pks) < 1000: # ensure that there is a minimum number of peaks to be collected to build noise model
        ctrl_sig = np.random.poisson(fl_baseline)/fl_baseline-1
        
        # apply mf
        ctrl_sig = np.convolve(ctrl_sig, rev_mf_template, 'same')
        ctrl_sig = ctrl_sig[:len(fl_baseline)]

        # apply hp filter
        if 'hp_cut' in pars:
            ctrl_sig = proc.filt_butter_hp(sig = ctrl_sig, fs = fl_fs, fc = pars['hp_cut'], order = 2)

        pks = np.concatenate([pks, scipysig.find_peaks(ctrl_sig, height = 0, distance = math.ceil(pars['min_event_interval']*1e-3*fl_fs))[1]['peak_heights']])
        niter += 1

    # calculate probability of detecting false positive peaks given target false positive rate
    pk_det_prob = pars['fp_rate']/(len(pks)/(1./fl_fs*len(fl_baseline)*niter))

    # determine detection threshold using kernel density estimation
    KD = sm.nonparametric.KDEUnivariate(pks)
    KD.fit(bw = 'scott')
    thresh_idx = np.argwhere(KD.sf<pk_det_prob)[0][0]
    # take average
    thresh = 0.5*(KD.support[thresh_idx-1]+KD.support[thresh_idx])

    # if minimum threshold condition is applied, calculate expected false positive rate
    if 'min_threshold' in pars and pars['min_threshold'] is not None and thresh<pars['min_threshold']:    
        support_idxs = np.argwhere(KD.support>=pars['min_threshold'])
        if len(support_idxs):
            fp_rate = KD.sf[support_idxs[0][0]]
        else:
            # false positive rate cannot be determined accurately, thus it is approximated to 0
            fp_rate = 0
        thresh = pars['min_threshold']
    else:
        fp_rate = pars['fp_rate']

    # detect events
    events_df = pd.DataFrame()
    
    detected_evts = scipysig.find_peaks(mf_sig, height = thresh, distance = math.ceil(pars['min_event_interval']*1e-3*fl_fs))
    events_df['mf-peak-time'] = (1./fl_fs*detected_evts[0]).astype(np.float32)
    events_df['mf-peak-amplitude'] = (detected_evts[1]['peak_heights']).astype(np.float32)
    # baseline photon rate per 1/sig['imaging']['fs'] integration time
    events_df['photon-rate'] = fl_baseline[detected_evts[0]]
    # determine if event is isolated or not
    events_df['is_isolated'] = proc.isolated_evt(evts = detected_evts[0], spacing = int(math.ceil(pars['isolated_event_interval']*1e-3*fl_fs)))

    meas_window_nsamp = int(math.ceil(pars['min_event_interval']*1e-3*fl_fs))
    # extract fluorescence peak time and amplitude within the minimum event interval window by applying max to the samples within the window
    fl_peak_idxs = np.array([max(0, evt_idx-meas_window_nsamp)+ \
        np.nanargmax(fl_sig[max(0, evt_idx-meas_window_nsamp):min(evt_idx+meas_window_nsamp+1, len(fl_sig))]) for evt_idx in detected_evts[0]])
    events_df['fl-peak-time'] = (1./fl_fs*fl_peak_idxs).astype(np.float32)
    if len(fl_peak_idxs):
        events_df['fl-peak-amplitude'] = fl_sig[fl_peak_idxs].astype(np.float32)
    else:
        events_df['fl-peak-amplitude'] = np.array([], dtype = np.float32)

    # extract fluorescence centroid time within the minimum event interval; the centroid is calculated for fluorescence above 'min_threshold'
    if 'fl_centroid_threshold' in pars:
        fl_centroid_threshold = pars['fl_centroid_threshold']
    else:
        fl_centroid_threshold = 0
    fl_centroid_idxs = np.empty((len(detected_evts[0]),))
    for idx, evt_idx in enumerate(detected_evts[0]):
        crop_start_idx = max(0, evt_idx-meas_window_nsamp)
        crop_end_idx = min(evt_idx+meas_window_nsamp+1, len(fl_sig))
        cropped_fl = np.copy(fl_sig[crop_start_idx:crop_end_idx])
        # set invalid fluorescence samples to threshold level to avoid passing nans further
        cropped_fl[np.isnan(cropped_fl)] = fl_centroid_threshold
        weights = cropped_fl-fl_centroid_threshold
        # set weights to 0 for for samples below threshold
        weights[weights<0] = 0
        # calculate centroid index
        if sum(weights) > 0:
            fl_centroid_idxs[idx] = np.average(np.arange(crop_start_idx,crop_end_idx), weights = weights)
        else:
            # there are no samples that can be used for event centroid index estimation
            fl_centroid_idxs[idx] = np.nan

    # warning: should be used with care, may contain nans
    events_df['fl-centroid-time'] = (1./fl_fs*fl_centroid_idxs).astype(np.float32)
    events_df['fl-centroid-amplitude'] = np.array([(fl_sig[int(idx)] if not np.isnan(idx) else np.nan) for idx in fl_centroid_idxs], dtype = np.float32)

    # output
    out = {
        'events': events_df,
        'threshold': thresh, # given or adjusted detection threshold
        'template': norm_mf_template.astype(np.float32),
        'valid_time': (len(mf_sig)-np.sum(np.isnan(mf_sig)))/fl_fs,
        'mf_sig': mf_sig.astype(np.float32)
    }

    return out

def _cross_correlation_using_fft(x, y):
    f1 = fft(x)
    f2 = fft(np.flipud(y))
    cc = np.real(ifft(f1 * f2))
    return fftshift(cc)

def compute_shift(x, y):
    """
    Calculates shift in samples between x and y.
    
    x and y should have about the same lengths for this method to work.

    Parameters
    ----------
    x, y : 1D numpy.array

    Returns
    -------
    int
        Shift in samples.
        Shift < 0 means that y starts 'shift' time steps before x.
        Shift > 0 means that y starts 'shift' time steps after x.
    """
    ns = min(len(x), len(y))
    x = x[:ns]
    y = y[:ns]
    c = _cross_correlation_using_fft(x, y)
    assert len(c) == len(x)
    zero_index = int(len(x) / 2) - 1
    shift = zero_index - np.argmax(c)
    return shift



