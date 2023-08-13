%alignBergamoStack
%Aligns a scanimage slow stack with repeats acquired on the Bergamo, generating a single
%high-quality image
%Stacks should be 1024 x 1024
%!!! The stack should be split into 1 file per Z-position. This reduces memory load


[fn, dr] = uigetfile('*.tif*', 'Select ONE FILE PER PLANE', 'multiselect', 'on');
if isnumeric(fn)
    return % user abort
end
if ~iscell(fn)
    fn = {fn};
end
tiffFn = fullfile(dr,fn);
tiffFn = sort_nat(tiffFn);

[IMs,meta] = readTiff(tiffFn{1});
try
    eval(meta);
catch ME
    disp(ME.message);
end

numChannels = length(SI.hChannels.channelSave);
numFramesPerSlice = SI.hStackManager.framesPerSlice;
zStep = SI.hStackManager.actualStackZStepSize;
numZs = length(fn);
pixelSizeUm = abs(diff(SI.hRoiManager.imagingFovUm(1:2,1)))./SI.hRoiManager.pixelsPerLine;

for Zix = numZs:-1:1
    disp(['Aligning plane: ' int2str(Zix) ' of ' int2str(numZs)])
    [IMs,~] = readTiff(tiffFn{Zix});

    % data order is XYCT
    IMs = permute(single(IMs), [2 1 3]); %transpose image; needed for bidi artifact detection by normcorr
    IMs = reshape(IMs,size(IMs,1),size(IMs,2),numChannels,numFramesPerSlice);

    % align
    refIM(:,:,:,Zix) = alignMultiChannel(IMs);
end

%straighten out the reference stack
refIM = straightenStack(refIM);

%save reference image, one stack for each channel
disp('Saving...')

outputPath = [tiffFn{1}(1:end-4) '_REF.ome.tif'];
metadata = createMinimalOMEXMLMetadata(IM);
pixelSizeObj = ome.units.quantity.Length(java.lang.Double(pixelSizeUm), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeX(pixelSizeObj, 0);
metadata.setPixelsPhysicalSizeY(pixelSizeObj, 0);
pixelSizeZObj = ome.units.quantity.Length(java.lang.Double(zStep), ome.units.UNITS.MICROMETER);
metadata.setPixelsPhysicalSizeZ(pixelSizeZObj, 0);
bfsave(IM, outputPath, 'BigTiff', true)

function refPlane = alignMultiChannel(IMs)
Y = squeeze(makeHighPass(sum(IMs,3))); % data order is XYCTZV
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',100,'us_fac',4,'init_batch',40, 'correct_bidir', true);
refPlane = [];
[~,shifts,~,options, col_shift] = normcorre(Y,options_rigid);
for chIx= size(IMs,3):-1:1
    tmp= apply_shifts(squeeze(IMs(:,:,chIx,:)),shifts,options,0,0,0,col_shift);
    refPlane(:,:,chIx) = mean(tmp,3);
end
end

function IMout = straightenStack(IM)
Y = makeHighPass(sum(IM,3)); % data order is XYCZ
Z2fft = fft2(Y(:,:,:,1));
IMout = nan*IM;
for Z1 = 1:(size(IM, 4)-1)
    Z1fft = Z2fft;
    Z2fft = fft2(Y(:,:,:,Z1+1));
    output = dftregistration_clipped(Z1fft,Z2fft,4,5);
    IMout(:,:,:,Z1+1) = imtranslate(IMout(:,:,:,Z1+1), [output(4) output(3)]);
end
end

% function refIM = alignMultiChannelbyPlane(IMs)
% nVol = size(IMs,6);
% nZ = size(IMs, 5); % data order is XYCTZV
% [p,n] = fileparts(tiffFn);
% 
% msg = sprintf(['Aligning file ' n '\n#Stacks: ' int2str(nVol) ' # planes: ' int2str(nZ)]);
% msg = most.idioms.latexEscape(msg);
% disp(msg)
% 
% Y = makeHighPass(sum(IMs,3)); % data order is XYCTZV
% 
% options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',4,'init_batch',30, 'correct_bidir', true);
% refIM = [];
% for Z = nZ:-1:1
%     [~,shifts,~,options, col_shift] = normcorre(squeeze(Y(:,:,1,:,Z)),options_rigid);
%     for chIx= 1:size(IMs,3)
%         tmp= apply_shifts(squeeze(IMs(:,:,chIx,:,Z)),shifts,options,0,0,0,col_shift);
%         refIM(:,:,chIx,1,Z) = mean(tmp,3);
%     end
% end
% end

function IMHP = makeHighPass(IM1)
IMHP = IM1; %sum the channels to speed alignment
IMHP(IMHP<0) = 0;
IMHP(isnan(IMHP)) = 0;
IMHP = imgaussfilt(IMHP, [0.5 0.5])-imgaussfilt(IMHP, [4 4]);
end

function [IMs,meta, ImDescriptions] = readTiff(tiffFn)
hTiff= ScanImageTiffReader(tiffFn);
try
    ImDescriptions = hTiff.descriptions;
    meta = hTiff.metadata();
    IMs = hTiff.data;
catch ME
    most.idioms.safeDeleteObj(hTiff);
    ME.rethrow();
end
most.idioms.safeDeleteObj(hTiff);
end

% function saveTif(IM,fileName, metadata)
% [filepath,fileName,ext] = fileparts(fileName);
% 
% if exist(filepath,'dir')
%     fileName = fullfile(filepath,[fileName,ext]);
% else
%     fileName = [fileName,ext];
% end
% 
% % [fileName,filePath] = uiputfile('.tif','Save ReferenceStack',fileName);
% % if isnumeric(fileName)
% %     return %User abort
% % end
% 
% fileName = fullfile(filePath,fileName);
% disp(['Saving file: ' fileName])
% 
% imageDescriptions = {};
% for zIdx = 1:size(IM,3)
%     metadata.ReferenceStackFileVersion = 1;
%     metadata.z= (zIdx-1)*metadata.zStep;
%     imageDescriptions{zIdx} = jsonencode(metadata);
% end
% 
% isTransposed = true;
% writeTiff(fileName, IM, isTransposed, imageDescriptions);
% end
