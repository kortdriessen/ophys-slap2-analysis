import torch
import numpy as np

def reconstruct(A, phi, B, beta, subsampleMatrixInds, sparseHInds, sparseHVals, uniqueMotion, motInds, framesToUse=[]):
    numCycles = phi.shape[0]
    numSuperPixels = subsampleMatrixInds.shape[1]

    if not isinstance(framesToUse, np.ndarray):
        framesToUse = np.arange(numCycles)

    (numZs, dmdPixelsPerColumn, dmdPixelsPerRow) = B.shape
    nPixels = dmdPixelsPerColumn * dmdPixelsPerRow

    dataEst = torch.zeros(numSuperPixels,framesToUse.shape[0])

    refPixs = torch.from_numpy(subsampleMatrixInds[0])

    d = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
    c = torch.div((refPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
    r = refPixs % dmdPixelsPerColumn

    mIdxs = motInds[framesToUse]

    for m in np.unique(mIdxs):
        i = np.where(mIdxs == m)[0]
        t = framesToUse[i]

        newR = r+uniqueMotion[m,0]
        newC = c+uniqueMotion[m,1]
        newD = uniqueMotion[m,2] - 1

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[m,0].astype(int) * dmdPixelsPerRow + uniqueMotion[m,1].astype(int)

        H = torch.sparse_coo_tensor(sparseHIndsShifted,sparseHVals,(numSuperPixels,nPixels),dtype=torch.float32)

        X = torch.sparse.mm(H, A)

        dataEst[:,i] = X.float() @ (beta[t].reshape((1,len(t))) * phi[t,:].T.reshape((-1,len(t)))).float() + torch.as_tensor(B[newD.astype(np.uint8)][newR.int(),newC.int()].reshape((-1,1))) @ beta[t].reshape((1,len(t))).to(torch.float32)

    return torch.maximum(dataEst, torch.zeros_like(dataEst))

def reconstruct_patch(A, phi, B, beta, subsampleMatrixInds, sparseHInds, sparseHVals, uniqueMotion, motInds, patchCenter, patchLength, framesToUse=[]):
    numCycles = phi.shape[0]
    numSuperPixels = subsampleMatrixInds.shape[0]

    if not isinstance(framesToUse, np.ndarray):
        framesToUse = np.arange(numCycles)

    (numZs, dmdPixelsPerColumn, dmdPixelsPerRow) = B.shape

    dataEst = torch.zeros(numSuperPixels,framesToUse.shape[0])

    buffer = 5

    refPixs = torch.from_numpy(subsampleMatrixInds[:,0]-1)
    refD = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
    refC = torch.div((refPixs - refD * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
    refR = refPixs % dmdPixelsPerColumn

    openPixs = torch.from_numpy(sparseHInds[1,:])
    openR = torch.div(openPixs, dmdPixelsPerRow, rounding_mode='floor') - (patchCenter[0]-patchLength-buffer)
    openC = openPixs % dmdPixelsPerRow - (patchCenter[1]-patchLength-buffer)

    mIdxs = motInds[framesToUse]

    A_crop = A.reshape((dmdPixelsPerColumn, dmdPixelsPerRow,-1))[patchCenter[0]-patchLength-buffer:patchCenter[0]+patchLength+buffer+1,patchCenter[1]-patchLength-buffer:patchCenter[1]+patchLength+buffer+1,:].reshape(((patchLength*2+buffer*2+1)**2,-1))

    for m in np.unique(mIdxs):
        i = np.where(mIdxs == m)[0]
        t = framesToUse[i]

        newR = refR+uniqueMotion[m,0]
        newC = refC+uniqueMotion[m,1]
        newD = uniqueMotion[m,2] - 1

        newOpenR = openR + uniqueMotion[m,0].astype(int)
        newOpenC = openC + uniqueMotion[m,1].astype(int)

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = newOpenR * (patchLength*2+buffer*2+1) + newOpenC

        validInds = torch.nonzero((newOpenR >= 0) & (newOpenR < patchLength*2+buffer*2+1) & (newOpenC >= 0) & (newOpenC < patchLength*2+buffer*2+1), as_tuple=True)

        H = torch.sparse_coo_tensor(sparseHIndsShifted[:,validInds[0]],sparseHVals[validInds[0]],(numSuperPixels,(patchLength*2+buffer*2+1)**2),dtype=torch.float32)

        X = torch.sparse.mm(H, A_crop)

        dataEst[:,i] = X.float() @ (beta[t].reshape((1,len(t))) * phi[t,:].T.reshape((-1,len(t)))).float() + torch.as_tensor(B[newD.astype(np.uint8)][newR.int(),newC.int()].reshape((-1,1))) @ beta[t].reshape((1,len(t))).to(torch.float32)

    return torch.maximum(dataEst, torch.zeros_like(dataEst))

# create function make2DGaussian(mean,var_x,var_y,theta,dmdPixelsPerColumn,dmdPixelsPerRow)
# returns a dmdPixelsPerColumn by dmdPixelsPerRow matrix with a 2D Gaussian with the given parameters normalized to a peak of one and truncated below values of 0.5
def make2DGaussian(mean,var_x,var_y,theta,dmdPixelsPerColumn,dmdPixelsPerRow):
    theta = min(max(theta,-90),90) * torch.pi / 180
    x = torch.arange(dmdPixelsPerRow)
    y = torch.arange(dmdPixelsPerColumn)
    X,Y = torch.meshgrid(x,y, indexing='ij')
    X = X - mean[0]
    Y = Y - mean[1]
    a = torch.cos(theta)**2/(2*var_x) + torch.sin(theta)**2/(2*var_y)
    b = -torch.sin(2*theta)/(4*var_x) + torch.sin(2*theta)/(4*var_y)
    c = torch.sin(theta)**2/(2*var_x) + torch.cos(theta)**2/(2*var_y)
    tmp = -(a*X**2 + 2*b*X*Y + c*Y**2)
    gauss = torch.maximum(torch.exp(tmp) - 1e-3, torch.zeros(1))

    return gauss.T / torch.sum(gauss.T)