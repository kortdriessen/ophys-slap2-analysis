#Imports
import tifffile as tiff
import numpy as np
import pickle
import pandas as pd
from scipy.sparse.linalg import svds
import scipy.cluster.hierarchy as sch

#initatial parameters
nFrames = 10000
kClusters = 20
savePath = 'C:/temp/'

#Read Tiff
with tiff.TiffFile('C:/temp/scan_00003_20240816_115328_REGISTERED_DOWNSAMPLED-2x.tif') as tif:
    movie = tif.asarray()
    l,w,h = movie.shape

#Format Data For SVD
X_data = np.reshape(movie, (l,w*h)).T 
X_mean = np.nanmean(X_data,axis=1) #row averages
B = X_data - np.tile(X_mean, (X_data.shape[1],1)).T
Btruncated = B[:,:nFrames]
BnanMat =  np.isnan(Btruncated)
keepRows= np.nanmean(BnanMat, axis=1)<0.2
B_tilde = Btruncated[keepRows,:]
B_tilde = np.where(np.isnan(B_tilde),0,B_tilde)

#Run SVDS
U, S, Vt = svds(B_tilde, k=20)
print(
    'Shapes: \n',
    'U Shape ==> ', U.shape, '\n',
    'S Shape ==> ', S.shape, '\n',
    'Vt Shape ==> ', Vt.shape
)
'''
#         **Optional**  
#        View Dendrogram
#Compute Hierarchical Clustering Based off "EigenDendrites"
linkageMatrix = sch.linkage(Vt.T, method = 'complete', metric='euclidean')
plt.figure()
dendrogram = sch.dendrogram(
    linkageMatrix,
    # p=20,
    # truncate_mode='level',
    count_sort=True,
    orientation='right')
plt.xlabel('Samples')
plt.ylabel('Distance')  
plt.show()
kClusters = 20
clusterLabels = sch.fcluster(linkageMatrix, kClusters, criterion='maxclust')
'''
#Extract Clusters
clusterLabels = sch.fclusterdata(
    Vt.T,
    t=kClusters,
    criterion='maxclust',
    method='single')#using single because assuming we want outlier detection

#recover frames in data from clusters
clustersToOrganize = np.unique(clusterLabels)
imgDict= {}
for idx in clustersToOrganize:
    dataInCluster = Btruncated[:,clusterLabels == idx]
    avgImg = np.nanmean(dataInCluster, axis=1)
    imgDict[idx] = np.reshape(avgImg, (w,h))

#save dictionary for labeling
with open('C:/temp/imgDict.pkl','wb') as f:
    pickle.dump(imgDict, f)
