"""Extract neural activity from a fluorescence trace using OASIS,
an active set method for sparse nonnegative deconvolution
Created on Mon Apr 4 18:21:13 2016
@author: Johannes Friedrich
"""

import numpy as np
cimport numpy as np
cimport cython
from libc.math cimport sqrt, log, exp, fmax, fabs, isnan
from scipy.optimize import fminbound, minimize
from cpython cimport bool
from libcpp.vector cimport vector

ctypedef np.float_t DOUBLE

cdef struct Pool:
    DOUBLE v
    DOUBLE w
    Py_ssize_t t
    Py_ssize_t l


@cython.cdivision(True)
def oasisAR1(np.ndarray[DOUBLE, ndim=1] y, DOUBLE g, DOUBLE lam=0, DOUBLE s_min=0):
    """ Infer the most likely discretized spike train underlying an AR(1) fluorescence trace

    Solves the sparse non-negative deconvolution problem
    min 1/2|c-y|^2 + lam |s|_1 subject to s_t = c_t-g c_{t-1} >=s_min or =0

    Parameters
    ----------
    y : array of float
        One dimensional array containing the fluorescence intensities with
        one entry per time-bin.
    g : float
        Parameter of the AR(1) process that models the fluorescence impulse response.
    lam : float, optional, default 0
        Sparsity penalty parameter lambda.
    s_min : float, optional, default 0
        Minimal non-zero activity within each bin (minimal 'spike size').

    Returns
    -------
    c : array of float
        The inferred denoised fluorescence signal at each time-bin.
    s : array of float
        Discretized deconvolved neural activity (spikes)

    References
    ----------
    * Friedrich J and Paninski L, NIPS 2016
    * Friedrich J, Zhou P, and Paninski L, PLOS Computational Biology 2017
    """

    cdef:
        Py_ssize_t i, j, k, t, T
        DOUBLE tmp
        np.ndarray[DOUBLE, ndim = 1] c, s
        vector[Pool] P
        Pool newpool

    T = len(y)
    # [value, weight, start time, length] of pool
    newpool.v, newpool.w, newpool.t, newpool.l = y[0] - lam * (1 - g), 1, 0, 1
    P.push_back(newpool)
    i = 0  # index of last pool
    t = 1  # number of time points added = index of next data point
    while t < T:
        # add next data point as pool
        newpool.v = y[t] - lam * (1 if t == T - 1 else (1 - g))
        newpool.w, newpool.t, newpool.l = 1, t, 1
        P.push_back(newpool)
        t += 1
        i += 1
        while (i > 0 and (isnan(P[i].v) or  # backtrack until violations fixed
               (P[i - 1].v / P[i - 1].w * g**P[i - 1].l + s_min > P[i].v / P[i].w))):
            i -= 1
            # merge two pools
            if not isnan(P[i + 1].v):
                P[i].v += P[i + 1].v * g**P[i].l
                P[i].w += P[i + 1].w * g**(2 * P[i].l)
            P[i].l += P[i + 1].l
            P.pop_back()
    # construct c and s
    c = np.empty(T)
    s = np.zeros(T)
    for j in range(i + 1):
        tmp = fmax(P[j].v, 0) / P[j].w
        if j > 0:
            k = P[j].t - 1
            s[P[j].t] = tmp - c[k] * g
            while isnan(y[k]):
                c[k] = np.nan
                s[k] = np.nan
                k -= 1
        for k in range(P[j].l):
            c[k + P[j].t] = tmp
            tmp *= g
    return c, s
