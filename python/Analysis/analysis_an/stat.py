from __future__ import division
import numpy as np
from scipy import stats
import collections
import warnings
from seaborn.algorithms import bootstrap as sb_bootstrap
import bottleneck as bn # speeds up certain operations on numpy arrays with nan

def h_fpp(H):
    """Evaluate the significance of an H score.
    The H test is an extension of the Z_m^2 or Rayleigh tests for
    uniformity on the circle. These tests estimate the Fourier coefficients
    of the distribution and compare them with the values predicted for
    a uniform distribution, but they require the user to specify the number
    of harmonics to use. The H test automatically selects the number of
    harmonics to use based on the data. 
    Arguments
    ---------
    H : float
        The H value to evaluate
    Returns
    -------
    fpp : float
        The probability of an H score this large arising from sampling a
        uniform distribution.
    Reference
    ---------
    de Jager, O. C., Swanepoel, J. W. H, and Raubenheimer, B. C., "A
    powerful test for weak periodic signals of unknown light curve shape
    in sparse data", Astron. Astrophys. 221, 180-190, 1989.
    """
    # These values are obtained by fitting to simulations.
    a = 0.9999755
    b = 0.39802
    c = 1.210597
    d = 0.45901
    e = 0.0022900

    if H<=23:
        return a*np.exp(-b*H)
    elif H<50:
        return c*np.exp(-d*H+e*H**2)
    else:
        return 4e-8
        # This comes up too often to raise an exception
        raise ValueError("H=%g>50 not supported; false positive probability less than 4*10**(-8)" % H)

def h_test(events):
    """Apply the H test for uniformity on [0,1).
    The H test is an extension of the Z_m^2 or Rayleigh test for
    uniformity on the circle. These tests estimate the Fourier coefficients
    of the distribution and compare them with the values predicted for
    a uniform distribution, but they require the user to specify the number
    of harmonics to use. The H test automatically selects the number of
    harmonics to use based on the data. The returned statistic, H, has mean
    and standard deviation approximately 2.51, but its significance should
    be evaluated with the routine h_fpp. This is done automatically in this
    routine.
    Arguments
    ---------
    events : array-like
        events should consist of an array of values to be interpreted as
        values modulo 1. These events will be tested for statistically
        significant deviations from uniformity.
    Returns
    -------
    H : float
        The raw score. Larger numbers indicate more non-uniformity.
    M : int
        The number of harmonics that give the most significant deviation
        from uniformity.
    fpp : float
        The probability of an H score this large arising from sampling a
        uniform distribution.
    Reference
    ---------
    de Jager, O. C., Swanepoel, J. W. H, and Raubenheimer, B. C., "A
    powerful test for weak periodic signals of unknown light curve shape
    in sparse data", Astron. Astrophys. 221, 180-190, 1989.
    """
    max_harmonic = 20
    ev = np.reshape(events, (-1,))
    cs = np.sum(np.exp(2.j*np.pi*np.arange(1,max_harmonic+1)*ev[:,None]),axis=0)/len(ev)
    Zm2 = 2*len(ev)*np.cumsum(np.abs(cs)**2)
    Hcand = (Zm2 - 4*np.arange(1,max_harmonic+1) + 4)
    M = np.argmax(Hcand)+1
    H = Hcand[M-1]
    fpp = h_fpp(H)
    return (H, M, fpp)

def calc_poisson_process_rate(n_evts, trec, conf_level = 95):
    """
    Calculates the mean rate and confidence intervals of a poisson random process.

    Parameters
    ----------
    n_evts : int
        Number of observed events within the interval 'trec'.
    trec : float
        Observation duration in [s].
    conf_level : float
        Confidence level in percents (0, 100).

    Returns
    -------
    tuple of (float, tuple) 
        Mean rate and confidence interval.
    """
    return n_evts/trec, (stats.chi2.ppf((1-conf_level/100.)/2., 2*n_evts)/2./trec, stats.chi2.ppf(1-(1-conf_level/100.)/2., 2*(n_evts+1))/2./trec)

"""
    Module for computing The Hartigans' dip statistic
    The dip statistic measures unimodality of a sample from a random process.
    See: 
    Hartigan, J. A.; Hartigan, P. M. The Dip Test of Unimodality. The Annals 
    of Statistics 13 (1985), no. 1, 70--84. doi:10.1214/aos/1176346577. 
    http://projecteuclid.org/euclid.aos/1176346577.
"""

# ==================================================
# dip unimodality test
# implementation by Johannes Bauer: https://github.com/tatome/dip_test/

def _gcm_(cdf, idxs):
    work_cdf = cdf
    work_idxs = idxs
    gcm = [work_cdf[0]]
    touchpoints = [0]
    while len(work_cdf) > 1:
        distances = work_idxs[1:] - work_idxs[0]
        slopes = (work_cdf[1:] - work_cdf[0]) / distances
        minslope = slopes.min()
        minslope_idx = np.where(slopes == minslope)[0][0] + 1
        gcm.extend(work_cdf[0] + distances[:minslope_idx] * minslope)
        touchpoints.append(touchpoints[-1] + minslope_idx)
        work_cdf = work_cdf[minslope_idx:]
        work_idxs = work_idxs[minslope_idx:]
    return np.array(np.array(gcm)),np.array(touchpoints)

def _lcm_(cdf, idxs):
    g,t = _gcm_(1-cdf[::-1], idxs.max() - idxs[::-1])
    return 1-g[::-1], len(cdf) - 1 - t[::-1]

def _touch_diffs_(part1, part2, touchpoints):
    diff = np.abs((part2[touchpoints] - part1[touchpoints]))
    return diff.max(), diff

def dip(histogram = None, idxs = None):
    """
        Compute the Hartigans' dip statistic either for a histogram of
        samples (with equidistant bins) or for a set of samples.
    """
    if idxs is None:
        idxs = np.arange(len(histogram))
    elif histogram is None:
        h = collections.Counter(idxs)
        idxs = np.msort(h.keys())
        histogram = np.array([h[i] for i in idxs])
    else:
        if len(histogram) != len(idxs):
            raise ValueError("Need exactly as many indices as histogram bins.")
        if len(idxs) != len(set(idxs)):
            raise ValueError("idxs must be unique if histogram is given.")
        if not np.array_equal(np.msort(idxs), idxs):
            idxs_s = np.argsort(idxs)
            idx = np.asarray(idxs)[idxs_s]
            histogram = np.asarray(histogram)[idxs_s]

    cdf = np.cumsum(histogram, dtype=float)
    cdf /= cdf[-1]

    work_idxs = idxs
    work_histogram = np.asarray(histogram, dtype=float) / np.sum(histogram)
    work_cdf = cdf

    D = 0
    left = [0]
    right = [1]

    while True:
        left_part, left_touchpoints   = _gcm_(work_cdf - work_histogram, work_idxs)
        right_part, right_touchpoints = _lcm_(work_cdf, work_idxs)

        d_left, left_diffs   = _touch_diffs_(left_part, right_part, left_touchpoints)
        d_right, right_diffs = _touch_diffs_(left_part, right_part, right_touchpoints)

        if d_right > d_left:
            xr = right_touchpoints[d_right == right_diffs][-1]
            xl = left_touchpoints[left_touchpoints <= xr][-1]
            d  = d_right
        else:
            xl = left_touchpoints[d_left == left_diffs][0]
            xr = right_touchpoints[right_touchpoints >= xl][0]
            d  = d_left

        left_diff  = np.abs(left_part[:xl+1] - work_cdf[:xl+1]).max()
        right_diff = np.abs(right_part[xr:]  - work_cdf[xr:] + work_histogram[xr:]).max()

        if d <= D or xr == 0 or xl == len(work_cdf):
            the_dip = max(np.abs(cdf[:len(left)] - left).max(), np.abs(cdf[-len(right)-1:-1] - right).max())
            return the_dip/2, (cdf, idxs, left, left_part, right, right_part)
        else:
            D = max(D, left_diff, right_diff)

        work_cdf = work_cdf[xl:xr+1]
        work_idxs = work_idxs[xl:xr+1]
        work_histogram = work_histogram[xl:xr+1]

        left[len(left):] = left_part[1:xl+1]
        right[:0] = right_part[xr:-1]

def plot_dip(histogram=None, idxs=None):
    from matplotlib import pyplot as plt

    d,(cdf,idxs,left,left_part,right,right_part) = dip(histogram,idxs)


    plt.plot(idxs[:len(left)], left, color='red')
    plt.plot(idxs[len(left)-1:len(left)+len(left_part) - 1], left_part, color='gray')
    plt.plot(idxs[-len(right):], right, color='blue')
    plt.plot(idxs[len(cdf) - len(right) + 1 - len(right_part):len(cdf) - len(right) + 1], right_part, color='gray')

    the_dip = max(np.abs(cdf[:len(left)] - left).max(), np.abs(cdf[-len(right)-1:-1] - right).max())
    l_dip_idxs = np.abs(cdf[:len(left)] - left) == the_dip
    r_dip_idxs = np.abs(cdf[-len(right)-1:-1] - right) == the_dip
    print(the_dip/2, d)

    plt.vlines(x=idxs[:len(left)][l_dip_idxs], ymin=cdf[:len(left)][l_dip_idxs], ymax = cdf[:len(left)][l_dip_idxs] - the_dip)
    plt.vlines(x=idxs[-len(right):][r_dip_idxs], ymin=cdf[-len(right)-1:-1][r_dip_idxs], ymax = cdf[-len(right)-1:][r_dip_idxs] + the_dip)

    plt.plot(np.repeat(idxs,2)[1:], np.repeat(cdf,2)[:-1], color='black')
    plt.scatter(idxs, cdf)

    plt.show()

def crit_points(random_function, quantiles, sample_size, n_samples):
    """
        Compute the quantiles for the dip statistic for n_samples 
        samples of size sample_size from the random process given by 
        random_function.
        Parameters:
        random_function : a paramter-free function which returns random values.
        quantiles : a sequence of values between 0 and 1
        sample_size : the size of the samples to draw from random_function
        n_samples : the number of samples for which to compute dips
        Returns: a list such that the i'th value is the greatest dip observed
        such that the fraction of dips less than or equal to that value is less
        than the i'th value from quantiles.
    """
    data = [[random_function() for _ in range(sample_size)] for __ in range(n_samples)]
    dips = np.array([dip(idxs=samples)[0] for samples in data])
    
    return np.percentile(dips, [p * 100 for p in quantiles])
    
def _nan_ci(a, which = 95, axis = None):
    """
    Return a percentile range from an array of values ignoring nans.
    """
    p = 50 - which / 2, 50 + which / 2
    return np.nanpercentile(a, p, axis)

def bs_mean(a, ci = 95, niter = 10000):
    """
    Calculate mean value and boostrap confidence intervals.

    Parameters
    ----------
    a : 1D numpy.array
        Array of values.
    ci : float
        Confidence interval, default 95% CI.
    niter : int
        Number of bootstrap iterations.
    
    Returns
    -------
    tuple 
        mean, lower and upper confidence intervals as (mean, lower_ci, upper_ci).
    """
    # ignore runtime warnings for slices containing only nans
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category = RuntimeWarning)
        mean = bn.nanmean(a)
        bs = sb_bootstrap(a, niter = niter, func = bn.nanmean)
        lower_ci, upper_ci = _nan_ci(a = bs, which = ci)

    return mean, lower_ci, upper_ci

