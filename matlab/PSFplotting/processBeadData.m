function processBeadData

%select a .dat file
[fn, dr] = uigetfile('D:\Slap2Data\beads\*.tif*', 'Select a dataset to analyze');

% Read the metadata; determine Z planes and number of acquisitions per
% plane
A = ScanImageTiffReader([dr filesep fn]);
D = A.data;

fnmeta = dir([dr filesep fn(1:end-4) '*.meta']);
try
    load([dr filesep fnmeta.name], '-mat', 'stackZs')
    nPlanes = size(stackZs,1);
    nFramesPerPlane = size(D,3)./nPlanes;
    assert(nFramesPerPlane-round(nFramesPerPlane)<1e-10);
    dZ = diff(stackZs,1,1);
    assert(all(all(dZ==dZ(1))));
    zStep = dZ(1);
catch ME %there was no metadata file
    %ask user for number of images/plane
    valid = false;
    while ~valid
        nPlanes = input('How many planes are in this acquisition? >>');
        nFramesPerPlane = size(D,3)./nPlanes;
        valid = nFramesPerPlane-round(nFramesPerPlane)<1e-10;
        if ~valid
            warning('The number of planes does not divide the number of images!')
            continue
        end
        zStep = input('What is the Z step size (um)? >>');
        valid = zStep>0 && zStep<10 && isscalar(zStep);
        if ~valid
            warning('Put in a reasonable Z step!')
            continue
        end
    end
end

%average any acquisitions to obtain an X-by-Y-by-Z stack
im = squeeze(mean(reshape(D, size(D,1),size(D,2),nFramesPerPlane,nPlanes),3));

%characterizeBeads(im, 0.25, zStep, [dr filesep fn(1:end-4)]);
characterizeBeads(im, 0.2245, zStep, [dr filesep fn(1:end-4)]);
