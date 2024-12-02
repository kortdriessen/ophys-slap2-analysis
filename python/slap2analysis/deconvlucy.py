import numpy as np
from scipy import fft
from typing import Union, Tuple, List, Optional
from pathlib import Path

def deconvlucy(image: np.ndarray, 
               psf: np.ndarray, 
               num_iter: int = 10,
               dampar: Union[float, np.ndarray] = 0,
               weight: Optional[np.ndarray] = None,
               readout: Union[float, np.ndarray] = 0,
               subsmpl: int = 1) -> Union[np.ndarray, List[np.ndarray]]:
    """
    Deblur image using Lucy-Richardson method.
    
    Parameters:
    -----------
    image : np.ndarray
        Input image (N-dimensional array)
    psf : np.ndarray
        Point spread function
    num_iter : int, optional
        Number of iterations (default is 10)
    dampar : float or np.ndarray, optional
        Damping parameter for noise threshold (default is 0)
    weight : np.ndarray, optional
        Pixel quality weights (default is array of ones)
    readout : float or np.ndarray, optional
        Additive noise (default is 0)
    subsmpl : int, optional
        Subsampling factor (default is 1)
    
    Returns:
    --------
    np.ndarray or List[np.ndarray]
        Deconvolved image or list of intermediate results
    """
    
    # Input validation and initialization
    if not isinstance(image, np.ndarray):
        raise TypeError("Image must be a numpy array")
    
    # Convert image to float64 if needed
    if image.dtype != np.float64:
        image = image.astype(np.float64)
    
    # Initialize weights if not provided
    if weight is None:
        weight = np.ones_like(image)
    
    # Get image dimensions
    size_i = np.array(image.shape)
    size_psf = np.array(psf.shape)
    
    # Find non-singleton dimensions of PSF
    num_ns_dim = np.where(size_psf != 1)[0]
    
    # Prepare PSF and create OTF
    size_otf = size_i.copy()
    size_otf[num_ns_dim] = subsmpl * size_i[num_ns_dim]
    H = psf2otf(psf, size_otf)
    
    # Initialize J cell structure equivalent
    J = [image.copy(),  # Original image
         image.copy(),  # Current iteration
         np.zeros_like(image),  # Previous iteration
         np.zeros((np.prod(size_i) * subsmpl**len(num_ns_dim), 2))]  # Internal state
    
    # Create indexes for image according to sampling rate
    idx = tuple(slice(None) for _ in range(len(size_i)))  # Changed this line
    
    # If subsampling is needed, modify the relevant indices
    if subsmpl > 1:
        for k in num_ns_dim:
            temp_idx = np.repeat(np.arange(size_i[k]), subsmpl)
            idx = idx[:k] + (temp_idx,) + idx[k+1:]
    
    # Prepare parameters for iterations
    w_i = np.maximum(weight * (readout + J[0]), 0)
    J[1] = J[1][idx]  # Now using the properly constructed idx tuple
    scale = np.real(fft.ifftn(np.conj(H) * fft.fftn(weight[idx]))) + np.sqrt(np.finfo(float).eps)
    dampar22 = (dampar ** 2) / 2
    
    # L-R Iterations
    lambda_val = 2 * np.any(J[3] != 0)
    for k in range(int(lambda_val), int(lambda_val + num_iter)):
        # Make image predictions for next iteration
        if k > 2:
            lambda_val = (J[3][:,0].T @ J[3][:,1]) / (J[3][:,1].T @ J[3][:,1] + np.finfo(float).eps)
            lambda_val = np.clip(lambda_val, 0, 1)
        
        Y = np.maximum(J[1] + lambda_val * (J[1] - J[2]), 0)
        
        # Make core for the LR estimation
        CC = corelucy(Y, H, dampar22, w_i, readout, subsmpl, idx)
        
        # Determine next iteration image & apply positivity constraint
        J[2] = J[1].copy()
        J[1] = np.maximum(Y * np.real(fft.ifftn(np.conj(H) * CC)) / scale, 0)
        J[3] = np.column_stack([J[1].flatten() - Y.flatten(), J[3][:,0]])
    
    return J[1]

def psf2otf(psf: np.ndarray, outSize: Optional[Union[Tuple, np.ndarray]] = None) -> np.ndarray:
    """
    Convert point-spread function to optical transfer function.
    
    Parameters:
    -----------
    psf : np.ndarray
        The point-spread function array
    outSize : tuple or np.ndarray, optional
        Desired output size of the OTF array. Must not be smaller than PSF size
        in any dimension. If None, uses same size as PSF.
    
    Returns:
    --------
    np.ndarray
        The optical transfer function array
    
    Notes:
    ------
    To ensure the OTF is not altered due to PSF off-centering, this function:
    1. Post-pads the PSF array with zeros to match outSize
    2. Circularly shifts the PSF array so the central pixel reaches (0,0)
    3. Computes the FFT to get the OTF
    """
    
    # Input validation
    if not isinstance(psf, np.ndarray):
        raise TypeError("PSF must be a numpy array")
    
    psf = psf.astype(np.float64)
    if not np.all(np.isfinite(psf)):
        raise ValueError("PSF must contain only finite values")
    
    psfSize = np.array(psf.shape)
    
    # Handle outSize
    if outSize is None:
        outSize = psfSize
    else:
        outSize = np.array(outSize)
        if not np.all(np.isfinite(outSize)):
            raise ValueError("outSize must contain only finite values")
        if np.any(outSize < 0):
            raise ValueError("outSize elements must be non-negative")
            
        # Pad lengths to match
        if len(psfSize) > len(outSize):
            outSize = np.pad(outSize, (0, len(psfSize) - len(outSize)), 
                           mode='constant', constant_values=1)
        elif len(outSize) > len(psfSize):
            psfSize = np.pad(psfSize, (0, len(outSize) - len(psfSize)), 
                           mode='constant', constant_values=1)
            psf = psf.reshape(psfSize)
            
        if np.any(outSize < psfSize):
            raise ValueError("outSize cannot be smaller than PSF size in any dimension")
    
    # If PSF is all zeros, return zeros
    if np.all(psf == 0):
        return np.zeros(outSize)
    
    # Pad PSF to outSize
    padSize = outSize - psfSize
    padSize = [(0, p) for p in padSize]  # Convert to format needed by np.pad
    psf = np.pad(psf, padSize, mode='constant')
    
    # Circularly shift PSF
    shift = -(psfSize // 2)
    psf = np.roll(psf, shift, range(len(psfSize)))
    
    # Compute OTF
    otf = fft.fftn(psf)
    
    # Estimate number of operations in FFT
    nElem = np.prod(psfSize)
    nOps = 0
    for k in range(psf.ndim):
        nffts = nElem / psfSize[k]
        nOps += psfSize[k] * np.log2(psfSize[k]) * nffts
    
    # Discard imaginary part if it's within roundoff error
    if np.max(np.abs(np.imag(otf))) / np.max(np.abs(otf)) <= nOps * np.finfo(float).eps:
        otf = np.real(otf)
        
    return otf

def corelucy(Y: np.ndarray, 
             H: np.ndarray, 
             DAMPAR22: float, 
             wI: np.ndarray, 
             READOUT: Union[float, np.ndarray], 
             SUBSMPL: int, 
             idx: List[slice]) -> np.ndarray:
    """
    Accelerated Damped Lucy-Richardson Operator.
    
    Calculates function that when used with the scaled projected array
    produces the next iteration array that maximizes the likelihood that
    the entire suite satisfies the Poisson statistics.
    
    Parameters:
    -----------
    Y : np.ndarray
        Current estimate of the image
    H : np.ndarray
        Optical transfer function (OTF)
    DAMPAR22 : float
        Damping parameter squared divided by 2
    wI : np.ndarray
        Weighted input image
    READOUT : float or np.ndarray
        Readout noise
    SUBSMPL : int
        Subsampling factor
    idx : List[slice]
        Index list for subsampling
    
    Returns:
    --------
    np.ndarray
        Fourier transform of the ratio image
    """
    
    # Calculate reblurred image
    ReBlurred = np.real(fft.ifftn(H * fft.fftn(Y)))
    
    # 1. Resampling if needed
    if SUBSMPL != 1:
        # Reshape for binning
        shape = ReBlurred.shape
        new_shape = []
        for i, s in enumerate(shape):
            if i in idx:
                new_shape.extend([SUBSMPL, s // SUBSMPL])
            else:
                new_shape.append(s)
        ReBlurred = ReBlurred.reshape(new_shape)
        
        # Calculate mean along subsampled dimensions
        axes_to_mean = range(0, len(new_shape), 2)
        for axis in axes_to_mean:
            ReBlurred = np.mean(ReBlurred, axis=axis)
    
    # 2. Estimate for the next step
    ReBlurred = ReBlurred + READOUT
    ReBlurred[ReBlurred == 0] = np.finfo(float).eps
    AnEstim = wI / ReBlurred + np.finfo(float).eps
    
    # 3. Damping if needed
    if DAMPAR22 == 0:  # No Damping
        ImRatio = AnEstim[tuple(idx)]
    else:  # Damping of the image relative to DAMPAR22 = (N*sigma)^2
        gm = 10
        g = (wI * np.log(AnEstim) + ReBlurred - wI) / DAMPAR22
        g = np.minimum(g, 1)
        G = (g ** (gm-1)) * (gm - (gm-1) * g)
        ImRatio = 1 + G[tuple(idx)] * (AnEstim[tuple(idx)] - 1)
    
    return fft.fftn(ImRatio)