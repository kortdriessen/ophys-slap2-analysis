function [data, meta]= copyReadDeleteScanImageTiff(remotepath, localDir)
if nargin<2
    localDir = 'C:\temp\';
end
if ~exist(localDir, 'dir')
    mkdir(localDir);
end

assert(exist(remotepath, 'file'));
randName = [int2str(round(1e10*rand+0.1)) '.tif'];
localpath = [localDir randName];
copyfile(remotepath, localpath);
A = ScanImageTiffReader(localpath);
data = A.data();
if nargout>1
    meta = A.metadata;
end
delete(A);
delete(localpath);
end