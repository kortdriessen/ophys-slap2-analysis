import numpy as np
import torch
import matplotlib.pyplot as plt

from scipy import io as spio
from scipy import signal

import h5py

import cv2
# import os


####################################################
## Read in data and set relevant global variables ##
####################################################
dmdPixelsPerColumn = 768
dmdPixelsPerRow = 768

print('Reading in data... ',end='')
f = h5py.File('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/INFER_MOTION_OUT.mat','r')

motionFull = np.array(f.get('inferMotionOut/motion')).T
brightness = np.array(f.get('inferMotionOut/brightness')).T
f.close()

f = h5py.File('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/trueA.mat','r')
trueA = np.array(f.get('A')).T
f.close()

f = spio.loadmat('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/simulated_activity_SLAP2_data_triStrip.mat')
truePhi = f['phi'].T
Y = f['spData'] / 100
sparseMaskInds = np.int64(f['sparseMaskInds'])
B = f['B'].T.reshape((11,dmdPixelsPerColumn,dmdPixelsPerRow)).astype(float)

print('done')

numCycles = Y.shape[1]
numSuperPixels = Y.shape[0]
spSizes = np.expand_dims(np.unique(sparseMaskInds[:,1], return_counts = True)[1],1)

nPixels = dmdPixelsPerColumn*dmdPixelsPerRow

motion = motionFull[:numCycles,:]
uniqueMotion, motInds = np.unique(motion,axis=0,return_inverse=True)

openPixs = sparseMaskInds[:,0]-1
d = openPixs // (dmdPixelsPerColumn*dmdPixelsPerRow)
c = (openPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)) // dmdPixelsPerColumn
r = openPixs % dmdPixelsPerColumn

###################################################
## Function definitions for movie reconstruction ##
###################################################

def reconstruct(A, phi, beta, framesToUse=[]):
    if not isinstance(framesToUse, np.ndarray):
        framesToUse = np.arange(numCycles)

    dataEst = torch.zeros(numSuperPixels,framesToUse.shape[0])

    # openPixs = torch.from_numpy(sparseMaskInds[:,0]-1)
    #
    # d = torch.div(openPixs, (dmdPixelsPerColumn*dmdPixelsPerRow),rounding_mode='floor')
    # c = torch.div((openPixs - d * (dmdPixelsPerColumn*dmdPixelsPerRow)), dmdPixelsPerColumn, rounding_mode='floor')
    # r = openPixs % dmdPixelsPerColumn

    mIdxs = motInds[framesToUse]

    for m in np.unique(mIdxs):
        i = np.where(mIdxs == m)[0]
        t = framesToUse[i]

        newR = r+uniqueMotion[m,0]
        newC = c+uniqueMotion[m,1]

        sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(dmdPixelsPerRow)+newC,1)),1)

        validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))

        H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))

        # H = torch.load(f"runtime_H_matrices/H{m}.pt")

        X = A @ (beta[t].reshape((1,len(t))) * phi[t,:].T.reshape((-1,len(t)))).float() + torch.from_numpy(B[6].T.reshape((dmdPixelsPerRow*dmdPixelsPerColumn,1))) @ beta[t].reshape((1,len(t)))

        dataEst[:,i] = torch.sparse.mm(H, X.to(torch.float32))

    return torch.maximum(dataEst, torch.zeros_like(dataEst))

def reconstruct_patch(A, phi, beta, patchCenter, patchLength, framesToUse=[]):
    if not isinstance(framesToUse, np.ndarray):
        framesToUse = np.arange(numCycles)

    dataEst = torch.zeros(numSuperPixels,framesToUse.shape[0])

    # openPixs = torch.from_numpy(sparseMaskInds[:,0]-1)

    buffer = 10

    # d = openPixs // (dmdPixelsPerColumn*dmdPixelsPerRow)
    c_patch = c - (patchCenter[1]-patchLength-buffer)
    r_patch = r - (patchCenter[0]-patchLength-buffer)

    mIdxs = motInds[framesToUse]

    A_crop = A.reshape((dmdPixelsPerColumn, dmdPixelsPerRow,-1))[patchCenter[0]-patchLength-buffer:patchCenter[0]+patchLength+buffer+1,patchCenter[1]-patchLength-buffer:patchCenter[1]+patchLength+buffer+1,:].reshape(((patchLength*2+buffer*2+1)**2,-1))

    for m in np.unique(mIdxs):
        i = np.where(mIdxs == m)[0]
        t = framesToUse[i]

        # try:
        #     H = torch.load(f"runtime_H_matrices/H{m}_{int(patchCenter[0])}x{int(patchCenter[1])}_{int(patchLength)}.pt")
        # except IOError:
        newR = r_patch+uniqueMotion[m,0];
        newC = c_patch+uniqueMotion[m,1];

        sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(patchLength*2+buffer*2+1)+newC,1)),1)

        validInds = np.nonzero((newR >= 0) & (newR < patchLength*2+buffer*2+1) & (newC >= 0) & (newC < patchLength*2+buffer*2+1))

        H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,(patchLength*2+buffer*2+1)**2))
        # H = torch.zeros([numSuperPixels,(patchLength*2+buffer*2+1)**2]).index_put_(tuple(k for k in torch.from_numpy(sparseHInds[validInds[0]].T.astype(int))), torch.ones(validInds[0].shape[0]))

            # torch.save(H, f"runtime_H_matrices/H{m}_{int(patchCenter[0])}x{int(patchCenter[1])}_{int(patchLength)}.pt")

        X = A_crop @ (beta[t].reshape((1,len(t))) * phi[t,:].T.reshape((-1,len(t)))).float() + torch.from_numpy(B[6].T[patchCenter[0]-patchLength-buffer:patchCenter[0]+patchLength+buffer+1,patchCenter[1]-patchLength-buffer:patchCenter[1]+patchLength+buffer+1].reshape(((patchLength*2+buffer*2+1)**2,1))) @ beta[t].reshape((1,len(t)))

        dataEst[:,i] = torch.sparse.mm(H, X.to(torch.float32))

    return torch.maximum(dataEst, torch.zeros_like(dataEst))

#############################
## Precalculate H matrices ##
#############################

# for m in np.unique(motInds):
#     t = np.where(motInds == m)[0]
#
#     newR = r+uniqueMotion[m,0]
#     newC = c+uniqueMotion[m,1]
#
#     sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(dmdPixelsPerRow)+newC,1)),1)
#
#     validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))
#
#     H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))
#
#     torch.save(H,f"runtime_H_matrices/H{m}.pt")

###########################
## Algorithm starts here ##
###########################

maxK = 5000
Afinal = torch.empty((nPixels,0))
phiFinal = torch.empty((numCycles,0))
phiVarFinal = torch.empty((numCycles,0))
trueIdxs = np.array([])

baseline = reconstruct(Afinal,phiFinal,torch.from_numpy(brightness)).detach().numpy()
Yhat = baseline

for k in range(maxK):
    print(f'Finding peak k={k}... ',end='')
    R = Y - Yhat
    R_norm = R / np.sqrt(baseline)

    gaussianGuess = np.expand_dims(np.exp(-(np.linspace(-3,3,7))**2/(2*2)),1)

    R_filt = signal.convolve2d(R_norm,gaussianGuess,mode='same')

    th = max(0,np.percentile(R_filt[:],75))
    Rfilt_thresh = R_filt - th
    Rfilt_thresh[Rfilt_thresh <= 0] = 0

    energy = Rfilt_thresh ** 2
    energyBackproj = np.zeros((dmdPixelsPerColumn,dmdPixelsPerRow,numCycles))

    for m in np.unique(motInds):
        t = np.where(motInds == m)[0]

        newR = r+uniqueMotion[m,0]
        newC = c+uniqueMotion[m,1]

        sparseHInds = np.concatenate((np.expand_dims(newR*(dmdPixelsPerRow)+newC,1),np.expand_dims(sparseMaskInds[:,1]-1,1)),1)

        validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))

        HT = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(nPixels,numSuperPixels))

        # HT = torch.transpose(torch.load(f"runtime_H_matrices/H{m}.pt"),0,1)

        energyBackproj[:,:,t] = torch.sparse.mm(HT,torch.from_numpy(energy[:,t] / spSizes).to(torch.float32)).reshape((dmdPixelsPerColumn,dmdPixelsPerRow,-1)).numpy()

    v = np.sum(energyBackproj,axis=2)

    # plt.figure(figsize=(10,10))
    # plt.imshow(v)
    # plt.colorbar()
    # plt.show()

    # if there are no more sources, stop adding sources
    if np.max(v) < 500:
        break

    # find peak of energy from the residual backprojected to the pixel space
    (colCenter, rowCenter) = np.unravel_index(np.argmax(v),v.shape)
    print(f'found at ({colCenter},{rowCenter})')

    trueIdx = np.argmax(np.squeeze(trueA.reshape((11,dmdPixelsPerColumn,dmdPixelsPerRow,200))[6,rowCenter,colCenter,:]))
    trueIdxs = np.append(trueIdxs,trueIdx)

    A_groundTruth = trueA.reshape((11,dmdPixelsPerRow,dmdPixelsPerColumn,200))[6,:,:,trueIdx].T.reshape((nPixels,1))
    phi_groundTruth = truePhi[:,trueIdx]

    initProfileWidth = 2
    initProfile = np.exp(-np.sum((np.mgrid[0:dmdPixelsPerColumn, 0:dmdPixelsPerRow] - np.expand_dims(np.array([colCenter, rowCenter]), (1, 2))) ** 2, 0) / (2 * initProfileWidth))

    initProfile[initProfile < 1e-5] = 0

    # update spatial profile using gradient descent
    # gradient descent parameters
    lr_A = 1e-4
    nEpochs = 5000
    ep = 1e-4
    beta1 = 0.9
    beta2 = 0.9

    trackedLosses_full = np.zeros((nEpochs,1))
    # A_validation = np.zeros((nEpochs,1))
    # phi_validation = np.zeros((nEpochs,1))

    bleach = torch.from_numpy(brightness.copy()).requires_grad_(True)

    A = torch.from_numpy(initProfile.reshape(nPixels,1) / np.sum(initProfile)).float()
    A.requires_grad = True

    allPhi = torch.zeros((numCycles, k+1))
    spOverlapMetric = torch.zeros((numCycles,1)) # measures overlap of A and superpixels (SPs)

    maskSize = 5
    mask = torch.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))
    mask[colCenter-maskSize:colCenter+maskSize+1,rowCenter-maskSize:rowCenter+maskSize+1] = 1
    mask = mask.reshape((nPixels,1))

    centroidCol = colCenter
    centroidRow = rowCenter

    G_A = torch.zeros_like(A)
    Gsum_A = torch.zeros_like(A)

    # lf = torch.nn.PoissonNLLLoss(log_input=False)

    j = 10
    miniBatchSize = 100

    def lf(theta, data):
        epsilon = 1e-8
        return torch.mean(theta - data * torch.log(theta+epsilon))

    print(f"Fitting component k={k}... ",end='')

    for epoch in range(nEpochs):

        fullACurr = torch.cat((Afinal,A),dim=1)

        with torch.no_grad():
            if epoch % j == 0:
                epMat = 1e-8 * torch.eye(fullACurr.shape[1])
                # print(f"Epoch {epoch}")
                for m in np.unique(motInds):
                    t = np.where(motInds == m)[0]

                    newR = r+uniqueMotion[m,0]
                    newC = c+uniqueMotion[m,1]

                    sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(dmdPixelsPerRow)+newC,1)),1)

                    validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))

                    H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))

                    # H = torch.load(f"runtime_H_matrices/H{m}.pt")

                    tmp = torch.sparse.mm(H,fullACurr)
                    tmpDot = tmp.T @ tmp
                    spOverlapMetric[t] = torch.diag(tmpDot)[k]

                    allPhi[t,:] = (torch.from_numpy((Y[:,t]-baseline[:,t])/brightness[t].reshape((-1,len(t)))).T.float() @ tmp @ torch.linalg.inv(tmpDot + epMat).T)

                predFull = reconstruct_patch(fullACurr,allPhi,bleach,(colCenter,rowCenter),maskSize)
                trackedLosses_full[epoch // j] = lf(predFull[predFull > 0],torch.from_numpy(Y)[predFull > 0].float()).numpy()
                # print(f"Loss Full = {trackedLosses_full[epoch // j]}")

                if epoch // j > 0 and (abs(trackedLosses_full[epoch // j] - trackedLosses_full[epoch // j - 1]) < 1e-5):
                    print(f"CONVERGED! ({epoch} epochs)")
                    break

        # goodFrames = np.argwhere(np.squeeze(spOverlapMetric.numpy()) > np.percentile(spOverlapMetric.numpy(),20)).T[0]

        framesToUse = np.random.choice(numCycles, miniBatchSize, replace=False)
        pred = reconstruct_patch(fullACurr,allPhi,bleach,(centroidCol,centroidRow),maskSize,framesToUse)
        # pred[pred==0] = torch.from_numpy(Y[:,framesToUse])[pred==0].float()

        loss = lf(pred[pred > 0],torch.from_numpy(Y[:,framesToUse])[pred > 0].float())

        loss.backward()

        with torch.no_grad():
            # print(A.grad[A.grad > 0])
            Gsum_A = Gsum_A * beta1 + (1-beta1) * A.grad
            G_A = G_A * beta2 + (1-beta2) * (A.grad) ** 2
            A_update = lr_A * Gsum_A / (torch.sqrt(G_A)+ep)

            # A_validation[epoch] = np.corrcoef(A.numpy().T,A_groundTruth.T)[0,1]
            # phi_validation[epoch] = np.corrcoef(allPhi[:,k].numpy()[goodFrames],phi_groundTruth[goodFrames])[0,1]

            A = torch.maximum(A - A_update * mask, torch.zeros_like(A))
            output = cv2.connectedComponentsWithStats((A.numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)) > 1e-5).astype(np.int8),connectivity=4)
            A[output[1].reshape(A.shape) != np.argmax(output[2][1:,cv2.CC_STAT_AREA])+1] = 0

            A = (A / torch.sum(A)).requires_grad_(True)
            # A -= A_update * mask
            # A.grad.zero_()
            # print(A_update[A_update != 0])

            centroidRow = np.round(np.sum((np.sum(A.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)),axis=0) * np.arange(dmdPixelsPerRow)) \
                                          / np.sum(A.detach().numpy()))).astype(int)
            centroidCol = np.round(np.sum((np.sum(A.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow)),axis=1) * np.arange(dmdPixelsPerRow)) \
                                          / np.sum(A.detach().numpy()))).astype(int)

        mask = torch.zeros((dmdPixelsPerColumn, dmdPixelsPerRow))
        mask[centroidCol-maskSize:centroidCol+maskSize+1,centroidRow-maskSize:centroidRow+maskSize+1] = 1
        mask = mask.reshape((nPixels,1))


    Afinal = torch.cat((Afinal,A.detach()),dim=1)
    phiFinal = allPhi.detach()

    Yhat = reconstruct(Afinal,phiFinal,bleach).detach().numpy()
    phiVarFinal = torch.zeros_like(phiFinal)

    epMat = 1e-8 * torch.eye(Afinal.detach().shape[1])
    for t in range(numCycles):

        newR = r+np.round(motion[t,0])
        newC = c+np.round(motion[t,1])

        sparseHInds = np.concatenate((np.expand_dims(sparseMaskInds[:,1]-1,1),np.expand_dims(newR*(dmdPixelsPerRow)+newC,1)),1)

        validInds = np.nonzero((newR >= 0) & (newR < dmdPixelsPerColumn) & (newC >= 0) & (newC < dmdPixelsPerRow))

        H = torch.sparse_coo_tensor(sparseHInds[validInds[0],:].T,torch.ones(validInds[0].shape[0]),(numSuperPixels,nPixels))

        # H = torch.load(f"runtime_H_matrices/H{motInds[t]}.pt")

        tmp = torch.sparse.mm(H,Afinal.detach())
        tmpDot = tmp.T @ tmp
        tmpDotInv = torch.linalg.inv(tmpDot + epMat)

        phiVarFinal[t,:] = torch.diag(tmpDotInv @ tmp.T @ torch.diag(torch.from_numpy(Yhat[:,t])).float() @ tmp @ tmpDotInv.T) / (brightness[t]**2)

    # if k > 2:
    #     break
    #
    # plt.figure(figsize=(15,20))
    # plt.subplot(221)
    # plt.imshow(trueA.reshape((11,dmdPixelsPerRow,dmdPixelsPerColumn,200))[6,rowCenter-20:rowCenter+20,colCenter-20:colCenter+20,trueIdx].T)
    # plt.title('ground truth')
    # plt.colorbar(orientation='horizontal')
    # plt.subplot(222)
    # plt.imshow(Afinal[:,k].detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow))[colCenter-20:colCenter+20,rowCenter-20:rowCenter+20])
    # plt.title('estimate')
    # plt.colorbar(orientation='horizontal')
    #
    # ax1 = plt.subplot(413)
    # ax1.set_xlabel('frame')
    # ax1.set_ylabel('phi')
    # ax1.plot(truePhi[:,trueIdx] / np.max(truePhi[:200,trueIdx]) * np.max(phiFinal[:200,k].detach().numpy().T),color='black',alpha=0.6)
    # ax1.plot(phiFinal[:,k].detach().numpy().T,color='red')
    # ax1.fill_between(np.arange(numCycles),np.squeeze(phiFinal[:,k].detach().numpy().T-np.sqrt(phiVarFinal[:,k].numpy()).T),np.squeeze(phiFinal[:,k].detach().numpy().T+np.sqrt(phiVarFinal[:,k].numpy()).T),facecolor='r',alpha=.5)
    # ax1.tick_params(axis='y',labelcolor='red')
    # plt.legend(('ground truth','estimate'))
    #
    # plt.show()

# dir = "runtime_H_matrices/patches/"
# for f in os.listdir(dir):
#     os.remove(os.path.join(dir, f))

print('Found all sources! Saving results... ',end='')
np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/trueIdxs_profiling.npy',trueIdxs)
np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/phi_profiling.npy',phiFinal.detach().numpy())
np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/phiVar_profiling.npy',phiVarFinal.detach().numpy())
np.save('C:/Users/michael.xie/Documents/data/activityMotionSimulations/fullFOV/moreActivity/A_profiling.npy',Afinal.detach().numpy().reshape((dmdPixelsPerColumn, dmdPixelsPerRow,-1)))
print('FINISHED')