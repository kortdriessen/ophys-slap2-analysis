from __future__ import division
import numpy as np

def ideal_transfer_function(f, fL, fH, fDSP = 0):
    """
    Calculates ideal transfer function for the intan RHD2000 amplifier system.

    Parameters
    ----------
    f : float, 1D iterable
        Frequencies in [Hz] at which to calculate the response.
    fL: float
        Low-pass filter cutoff frequency in [Hz].
    fH : float
        High-pass filter cutoff frequency in [Hz].

    Returns
    -------
    tuple (gain, phase)

        Transfer function gain and phase in [deg].
    """
    f = np.atleast_1d(f)
    assert len(f.shape) == 1
    gain = 1./(np.sqrt(1+np.power(f/fH, 6))*np.sqrt(1+np.power(fL/f, 2)))
    phase = np.arctan(fL/f)-np.arctan((2*(f/fH)-np.power((f/fH), 3))/(1-2*np.power((f/fH), 2)))
    if fDSP:
        gain = gain/np.sqrt(1 + np.power(fDSP/f,2))
        phase = phase+np.arctan(fDSP/f);


    gain = 20 * np.log10(gain) # convert gain to decibels
    phase = np.unwrap(phase) # unwrap phase
    phase = (180./np.pi) * phase # convert phase from radians to degrees

    return gain, phase
    