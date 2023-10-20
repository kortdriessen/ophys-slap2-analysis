import numpy as np
import torch
import matplotlib.pyplot as plt

from scipy import io as spio
from scipy import signal

from skimage import io as skimio

import h5py

import cv2

import reconstruct
# import os

####################################################
## Read in data and set relevant global variables ##
####################################################
dir = 'Z:\\ophys\\SLAP2\\exp data\\Mice\\slap2_664318_2023-08-01_11-07-05\\FOV1\\simulation\\bright_0\\activity_0\\'

print('Reading in data... ',end='')

# load true phis
f = spio.loadmat(dir + '..\\..\\simulatedActivity_phi_0.mat')
truePhi = f['phi']

# read in simulated recording data
f = spio.loadmat(dir + 'simulated_SLAP2_data_activity.mat')
Y = f['spData'] / 100
sparseMaskInds = np.int64(f['sparseMaskInds'])
trueMotion = f['motFinal']

# read in inferMotion outputs
f = h5py.File(dir + 'INFER_MOTION_OUT.mat','r')
motion = np.array(f.get('inferMotionOut/motion')).T
brightness = np.array(f.get('inferMotionOut/brightness')).T
expectedData = np.array(f.get('inferMotionOut/expectedMatrix')).T
f.close()

# read in background reference stack
B = skimio.imread(dir + '..\\..\\..\\refStack_20230801_112813_DMD1_CONFIG3-REFERENCE_CH2.tif') / 100

# read in dilation 5 to 15 filter
f = spio.loadmat(dir + '..\\..\\..\\Dil5To15Filter.mat')
filter5To15 = f['filterfinal']

# set relevant global variables
numCycles = Y.shape[1]
numSuperPixels = Y.shape[0]
spSizes = np.unique(sparseMaskInds[:,1], return_counts = True)[1]

dmdPixelsPerColumn = B.shape[1]
dmdPixelsPerRow = B.shape[2]
numZs = B.shape[0]

nPixels = (dmdPixelsPerColumn)*(dmdPixelsPerRow)

uniqueMotion, motInds = np.unique(motion,axis=0,return_inverse=True)

# extract the median open pixel for each superpixel defined by the sparseMaskInds
subsampleMatrixInds = np.zeros((numSuperPixels,2))

for spIdx in range(numSuperPixels):
    currSpInds = np.where(sparseMaskInds[:,1] == spIdx+1)[0]
    currSpOpenPixs = sparseMaskInds[currSpInds,0]
    spRefPix = currSpOpenPixs[np.floor(len(currSpOpenPixs)/2).astype(int)]
    subsampleMatrixInds[spIdx,0] = spRefPix
    subsampleMatrixInds[spIdx,1] = spIdx+1

refPixs = torch.from_numpy(subsampleMatrixInds[:,0])
refD = torch.div(refPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
refC = torch.div((refPixs - refD * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
refR = refPixs % dmdPixelsPerColumn

# extract all open pixels and locations
openPixs = torch.from_numpy(sparseMaskInds[:,0]-1)
d = torch.div(openPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
c = torch.div((openPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
r = openPixs % dmdPixelsPerColumn

# create sparse convolution matrix, H
filterSize = filter5To15.shape[0]*filter5To15.shape[1]
sparseHInds = np.zeros((2,numSuperPixels*filterSize))
sparseHVals = np.zeros((numSuperPixels*filterSize,))

for spIdx in range(subsampleMatrixInds.shape[0]):
    tmpMap = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow))
    tmpMap[refR[spIdx].int()-filter5To15.shape[0]//2:refR[spIdx].int()+filter5To15.shape[0]//2+1,refC[spIdx].int()-filter5To15.shape[1]//2:refC[spIdx].int()+filter5To15.shape[1]//2+1] = torch.from_numpy(filter5To15)

    sparseHInds[0,spIdx*filterSize:(spIdx+1)*filterSize] = subsampleMatrixInds[spIdx,1] - 1
    sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize] = np.where(tmpMap.flatten() > 0)[0]

    sparseHVals[spIdx*filterSize:(spIdx+1)*filterSize] = tmpMap.flatten()[sparseHInds[1,spIdx*filterSize:(spIdx+1)*filterSize].astype(np.uint32)]

print('done')

###########################
## Algorithm starts here ##
###########################

maxK = 5000
Afinal = torch.empty((nPixels,0))
phiFinal = torch.empty((numCycles,0))
# phiVarFinal = torch.empty((numCycles,0))
trueIdxs = np.array([])

baseline = reconstruct.reconstruct(Afinal,phiFinal,B,torch.from_numpy(brightness),subsampleMatrixInds,sparseHInds,sparseHVals,uniqueMotion,motInds).detach().numpy()
Yhat = baseline

for k in range(maxK):
    print(f'Finding peak k={k}... ',end='')
    R = Y - Yhat
    R_norm = R / np.sqrt(baseline)

    gaussianGuess = np.expand_dims(np.exp(-(np.linspace(-3,3,7))**2/(2*2)),1)
    expoFilter = np.expand_dims(np.exp(np.linspace(-6,0,7)/2),0)

    R_filt = signal.convolve2d(R_norm,gaussianGuess,mode='same')
    R_filt = signal.convolve2d(R_filt,expoFilter,mode='same')

    th = max(0,np.percentile(R_filt[:],75))
    Rfilt_thresh = R_filt - th
    Rfilt_thresh[Rfilt_thresh <= 0] = 0

    energy = Rfilt_thresh ** 2
    energyBackproj = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow,numCycles))

    for m in np.unique(motInds):
        t = np.where(motInds == m)[0]

        newR = r+np.round(uniqueMotion[m,0])
        newC = c+np.round(uniqueMotion[m,1])
        newD = uniqueMotion[m,2] - 1

        sparseHMaskInds = torch.cat((torch.unsqueeze(torch.from_numpy(sparseMaskInds[:,1]-1),1),torch.unsqueeze(newR*(dmdPixelsPerRow)+newC,1)),1)

        validInds = torch.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow), as_tuple=True)

        H = torch.sparse_coo_tensor(sparseHMaskInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))

        tmp = torch.sparse.mm(torch.transpose(H,0,1),torch.from_numpy(energy[:,t] / np.expand_dims(spSizes,1)).to(torch.float32))

        energyBackproj[:,:,t] = torch.squeeze(tmp).reshape((dmdPixelsPerColumn,dmdPixelsPerRow,-1)).numpy()

    v = np.sum(energyBackproj,axis=2)

    # plt.figure(figsize=(10,10))
    # plt.imshow(v)
    # plt.colorbar()
    # plt.show()

    # if there are no more sources, stop adding sources
    if np.max(v) < 200:
        break

    # find peak of energy from the residual backprojected to the pixel space
    (colCenter, rowCenter) = np.unravel_index(np.argmax(v),v.shape)
    print(f'found at ({colCenter},{rowCenter})')

    # trueIdx = np.argmax(np.squeeze(trueA.reshape((11,dmdPixelsPerColumn,dmdPixelsPerRow,200))[6,rowCenter,colCenter,:]))
    # trueIdxs = np.append(trueIdxs,trueIdx)
    #
    # A_groundTruth = trueA.reshape((11,dmdPixelsPerRow,dmdPixelsPerColumn,200))[6,:,:,trueIdx].T.reshape((nPixels,1))
    # phi_groundTruth = truePhi[:,trueIdx]

    initProfileWidth = 2
    # initProfile = np.exp(-np.sum((np.mgrid[0:dmdPixelsPerColumn, 0:dmdPixelsPerRow] - np.expand_dims(np.array([colCenter, rowCenter]), (1, 2))) ** 2, 0) / (2 * initProfileWidth))
    #
    # initProfile[initProfile < 1e-5] = 0

    # update spatial profile using gradient descent
    # gradient descent parameters
    lr_A0 = 1e0 # A position
    lr_A1 = 1e1 # A length and width
    lr_A2 = 1e-1 # A tilt

    lr_phi = 8e1 # phi

    nEpochs = 5000
    ep = 1e-4
    beta1 = 0.9
    beta2 = 0.999

    trackedLosses_full = np.zeros((nEpochs,1))
    # A_validation = np.zeros((nEpochs,1))
    # phi_validation = np.zeros((nEpochs,1))

    bleach = torch.from_numpy(brightness.copy()).requires_grad_(True)

    params = torch.tensor([rowCenter,colCenter,3,3,0], dtype=torch.float32, requires_grad=True)

    # A = torch.from_numpy(initProfile.reshape(nPixels,1) / np.sum(initProfile)).float()
    # A.requires_grad = True

    allPhi = torch.zeros((numCycles, k+1), requires_grad=True)
    initPhi = torch.zeros((numCycles, k+1))
    spOverlapMetric = torch.zeros((numCycles,1)) # measures overlap of A and superpixels (SPs)

    maskSize = 15
    # mask = torch.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))
    # mask[colCenter-maskSize:colCenter+maskSize+1,rowCenter-maskSize:rowCenter+maskSize+1] = 1
    # mask = mask.reshape((nPixels,1))
    #
    # centroidCol = colCenter
    # centroidRow = rowCenter

    G_A = torch.zeros_like(params)
    Gsum_A = torch.zeros_like(params)

    G_phi = torch.zeros_like(allPhi)
    Gsum_phi = torch.zeros_like(allPhi)

    # lf = torch.nn.PoissonNLLLoss(log_input=False)

    j = 10
    miniBatchSize = 600

    def lf(theta, data):
        epsilon = 1e-8
        return torch.mean(theta - data * torch.log(theta+epsilon))

    print(f"Fitting component k={k}... ",end='')

    for epoch in range(nEpochs):

        A = reconstruct.make2DGaussian((params[0], params[1]), params[2],params[3],params[4],dmdPixelsPerColumn,dmdPixelsPerRow).reshape((nPixels,1))

        fullACurr = torch.cat((Afinal,A),dim=1)

        with torch.no_grad():
            if epoch == 0:
                epMat = 1e-8 * torch.eye(fullACurr.shape[1])
                print(f"Epoch {epoch}")
                for m in np.unique(motInds):
                    t = np.where(motInds == m)[0]

                    sparseHIndsShifted = sparseHInds.copy()
                    sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[m,0].astype(int) * dmdPixelsPerRow + uniqueMotion[m,1].astype(int)
                    H = torch.sparse_coo_tensor(sparseHIndsShifted,sparseHVals,(numSuperPixels,nPixels),dtype=torch.float32)

                    tmp = torch.sparse.mm(H,fullACurr)
                    tmpDotInv = torch.linalg.inv(tmp.T @ tmp)

                    allPhi[t,:] = torch.from_numpy((Y[:,t]-baseline[:,t])/brightness[t].reshape((-1,len(t)))).T.float() @ tmp @ tmpDotInv.T

                initPhi = allPhi.detach().clone()

            if epoch % j == 0:
                predFull = reconstruct.reconstruct_patch(fullACurr,allPhi,B,bleach,subsampleMatrixInds,sparseHInds,sparseHVals,uniqueMotion,motInds,(colCenter,rowCenter),maskSize)
                trackedLosses_full[epoch // j] = lf(predFull,torch.from_numpy(Y).float()).numpy()
                print(f"Loss Full = {trackedLosses_full[epoch // j]}")

                if epoch // j > 0 and (abs(trackedLosses_full[epoch // j] - trackedLosses_full[epoch // j - 1]) < 1e-7):
                    print(f"CONVERGED! ({epoch} epochs)")
                    break

        # goodFrames = np.argwhere(np.squeeze(spOverlapMetric.numpy()) > np.percentile(spOverlapMetric.numpy(),20)).T[0]

        framesToUse = np.random.choice(numCycles, miniBatchSize, replace=False)
        pred = reconstruct.reconstruct_patch(fullACurr,allPhi,B,bleach,subsampleMatrixInds,sparseHInds,sparseHVals,uniqueMotion,motInds,(colCenter,rowCenter),maskSize,framesToUse)
        # pred[pred==0] = torch.from_numpy(Y[:,framesToUse])[pred==0].float()

        loss = lf(pred,torch.from_numpy(Y[:,framesToUse]).float())

        loss.backward()

        # Agrad = torch.autograd.grad(loss,A)[0]

        with torch.no_grad():
            if epoch % 10 < 5:
                # print(A.grad[A.grad > 0])
                Gsum_A = Gsum_A * beta1 + (1-beta1) * params.grad
                G_A = G_A * beta2 + (1-beta2) * (params.grad) ** 2
                A_update = torch.tensor([lr_A0,lr_A0,lr_A1,lr_A1,lr_A2], dtype=torch.float32) * Gsum_A / (torch.sqrt(G_A)+ep)

                params = params - A_update
                params.requires_grad_(True)
                allPhi.grad.zero_()
                # A_validation[epoch] = np.corrcoef(A.numpy().T,A_groundTruth.T)[0,1]
                # phi_validation[epoch] = np.corrcoef(allPhi[:,k].numpy()[goodFrames],phi_groundTruth[goodFrames])[0,1]

                # A = torch.maximum(A - A_update * mask, torch.zeros_like(A))
                # output = cv2.connectedComponentsWithStats((A.numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)) > 1e-5).astype(np.int8),connectivity=4)
                # A[output[1].reshape(A.shape) != np.argmax(output[2][1:,cv2.CC_STAT_AREA])+1] = 0
                #
                # A = (A / torch.sum(A)).requires_grad_(True)
                # # A -= A_update * mask
                # # A.grad.zero_()
                # # print(A_update[A_update != 0])
                #
                # centroidRow = np.round(np.sum((np.sum(A.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)),axis=0) * np.arange(dmdPixelsPerRow)) \
                #                               / np.sum(A.detach().numpy()))).astype(int)
                # centroidCol = np.round(np.sum((np.sum(A.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)),axis=1) * np.arange(dmdPixelsPerRow)) \
                #                               / np.sum(A.detach().numpy()))).astype(int)

                # mask = torch.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))
                # mask[centroidCol-maskSize:centroidCol+maskSize+1,centroidRow-maskSize:centroidRow+maskSize+1] = 1
                # mask = mask.reshape((nPixels,1))
            if epoch % 10 >= 5:
                Gsum_phi = Gsum_phi * beta1 + (1-beta1) * allPhi.grad
                G_phi = G_phi * beta2 + (1-beta2) * (allPhi.grad) ** 2
                phi_update = lr_phi * Gsum_phi / (torch.sqrt(G_phi)+ep)

                allPhi = allPhi - phi_update
                allPhi.requires_grad_(True)
                params.grad.zero_()

    A = reconstruct.make2DGaussian((params[0], params[1]), params[2],params[3],params[4],dmdPixelsPerColumn,dmdPixelsPerRow).reshape((nPixels,1))

    Afinal = torch.cat((Afinal,A.detach()),dim=1)
    phiFinal = allPhi.detach()

    Yhat = reconstruct.reconstruct(Afinal,phiFinal,B,bleach,subsampleMatrixInds,sparseHInds,sparseHVals,uniqueMotion,motInds).detach().numpy()
    # phiVarFinal = torch.zeros_like(phiFinal)

    # epMat = 1e-8 * torch.eye(Afinal.detach().shape[1])
    # for t in range(numCycles):
    #
    #     newR = r+np.round(motion[t,0])
    #     newC = c+np.round(motion[t,1])
    #
    #     sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(dmdPixelsPerRow)+newC,1)),1)
    #
    #     validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))
    #
    #     H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))
    #
    #     # H = torch.load(f"runtime_H_matrices/H{motInds[t]}.pt")
    #
    #     tmp = torch.sparse.mm(H,Afinal.detach())
    #     tmpDot = tmp.T @ tmp
    #     tmpDotInv = torch.linalg.inv(tmpDot + epMat)
    #
    #     phiVarFinal[t,:] = torch.diag(tmpDotInv @ tmp.T @ torch.diag(torch.from_numpy(Yhat[:,t])).float() @ tmp @ tmpDotInv.T) / (brightness[t]**2)

    phiCorrs = np.zeros((truePhi.shape[-1],1))

    for idx in range(truePhi.shape[-1]):
        phiCorrs[idx] = np.corrcoef(truePhi[:,idx],phiFinal[:,k].detach().numpy())[1,0]

    trueIdx = np.argmax(phiCorrs)

    phiLS = torch.zeros((numCycles, Afinal.shape[1]))
    for m in np.unique(motInds):
        t = np.where(motInds == m)[0]

        sparseHIndsShifted = sparseHInds.copy()
        sparseHIndsShifted[1,:] = sparseHIndsShifted[1,:] + uniqueMotion[m,0].astype(int) * dmdPixelsPerRow + uniqueMotion[m,1].astype(int)
        H = torch.sparse_coo_tensor(sparseHIndsShifted,sparseHVals,(numSuperPixels,nPixels),dtype=torch.float32)

        tmp = torch.sparse.mm(H,A.detach())
        tmpDotInv = torch.linalg.inv(tmp.T @ tmp)

        phiLS[t,:] = torch.from_numpy(R[:,t]/brightness[t].reshape((-1,len(t)))).T.float() @ tmp @ tmpDotInv.T

    plt.figure(figsize=(20,10))
    plt.subplot(241)
    plt.imshow(Afinal[:,k].detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow))[colCenter-20:colCenter+21,rowCenter-20:rowCenter+21]);
    plt.title('estimate')
    plt.colorbar(orientation='horizontal')

    ax_loss = plt.subplot(2,4,(2,4))
    ax_loss.plot(trackedLosses_full[0:epoch // j],color='blue')
    ax_loss.set_ylabel('Poisson NLL Loss (on all frames)',color='blue')

    ax_phi = plt.subplot(2,4,(5,8))
    ax_phi.set_xlabel('frame')
    ax_phi.set_ylabel('phi')
    ax_phi.plot(truePhi[:,trueIdx] / np.max(truePhi[:200,trueIdx]) * np.max(phiFinal[:200,k].detach().numpy().T),color='black',alpha=0.6)
    ax_phi.plot(phiLS[:,k].detach().numpy().T,color='green',alpha=0.6)
    ax_phi.plot(initPhi[:,k].detach().numpy().T,color='blue',alpha=0.6)
    ax_phi.plot(phiFinal[:,k].detach().numpy().T,color='red',alpha=0.6)
    plt.legend(('ground truth','LS estimate','initial guess','Poisson estimate'))

    # save current plot as a png
    plt.savefig(dir + f'source{k}.png',dpi=300)
    plt.close()

    # plt.show()

# dir = "runtime_H_matrices/patches/"
# for f in os.listdir(dir):
#     os.remove(os.path.join(dir, f))

print('Found all sources! Saving results... ',end='')
# np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/trueIdxs_profiling.npy',trueIdxs)
np.savez(dir + '/extractedSources.npz',phi=phiFinal.detach().numpy(),A=Afinal.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow,-1)))
# np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/phiVar_profiling.npy',phiVarFinal.detach().numpy())
# np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/A_profiling.npy',Afinal.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow,-1)))
print('FINISHED')