function processBeadRefStack (imPath)

% Load one reference .tif file
if ~nargin
[fn,dr] = uigetfile('*.tif', 'Select a SLAP2 Reference Stack');
end

% Find the meta data of the file
A = ScanImageTiffReader([dr filesep fn]);
meta = A.descriptions();
metaS = jsondecode(meta{1});

% Get the XY pixel size
XYpixelSize = (metaS.dmdPixel2SampleTransform(1,1) + metaS.dmdPixel2SampleTransform(2,2))/2; % this doesn't work

% Find the Z pixel spacing in a robust way by finding the first and second
% distinct z values and then calculating the difference
newPlane = false;
pix =2;
channels = metaS.channel;
while ~newPlane
    metaS2 = jsondecode(meta{pix});
    newPlane = metaS2.z~=metaS.z;
    pix = pix+1;
    channels = unique([channels, metaS2.channel]);
end
dZ = abs(metaS2.z-metaS.z);

% Each reference stack contains intensities from multiple (usually 2) detector channels
nChan = length(channels);
% Select with channel you want to analyze and plot
analyzeChannel = 1;
chanInd = find(channels == analyzeChannel,1);
assert(~isempty(chanInd), 'The analysis channel wasn''t imaged');

% Extract intensity data and store as doubles
im = double(A.data());
im = im(:,:,analyzeChannel:nChan:end);

characterizeBeads(im, XYpixelSize, dZ, [dr filesep fn(1:end-4)]);