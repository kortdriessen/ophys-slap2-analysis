%alignBergamoStack
%Aligns a scanimage slow stack with repeats acquired on the Bergamo, generating a single
%high-quality image
%Stacks should be 1024 x 1024
%!!! The stack should be split into 1 file per Z-position. This reduces memory load


[fn, dr] = uigetfile('*.tif*', 'Select ONE FILE PER PLANE', 'multiselect', 'on');
% if isnumeric(fn)
%     return % user abort
% end
% if ~iscell(fn)
%     fn = {fn};
% end
tiffFn = fullfile(dr,fn);
tiffFn = sort_nat(tiffFn);

[IMs,meta] = readTiff(tiffFn{1});
metaLines = strsplit(meta, '\n');
for lineIx = 1:length(metaLines)
    try
        eval([metaLines{lineIx} ';']);
    catch
    end
end

numChannels = length(SI.hChannels.channelSave);
numFramesPerSlice = SI.hStackManager.framesPerSlice;
zStep = SI.hStackManager.actualStackZStepSize;
numZs = length(fn);
pixelSizeUm = abs(diff(SI.hRoiManager.imagingFovUm(1:2,1)))./SI.hRoiManager.pixelsPerLine;

[b,a] = butter(4, 0.05, 'high'); %butterworth filter for computing correlations
IM = []; IMc = []; IMsk = [];
for Zix = numZs:-1:1
    disp(['Aligning plane: ' int2str(Zix) ' of ' int2str(numZs)])
    [IMs,~] = readTiff(tiffFn{Zix});

    % data order is XYCT
    IMs = permute(single(IMs), [2 1 3]); %transpose image; needed for bidi artifact detection by normcorr
    IMs = reshape(IMs,size(IMs,1),size(IMs,2),numChannels,numFramesPerSlice);

    % align
    [IM(:,:,:,Zix), IMc(:,:,:,Zix), IMsk(:,:,:,Zix)] = alignMultiChannel(IMs, b, a);
end

%straighten out the reference stack
[IM, IMc, IMsk] = straightenStack(IM, IMc, IMsk);

%convolve correlation image in 3D and background subtract
IMc2 = squeeze(IMc(:,:,2,:)); %we use channel 2
IMc2 = IMc2 - imgaussfilt(IMc2, 29, 'padding', 'symmetric'); 
IMc2 = imgaussfilt3(IMc2, [4 4 1]);
IMc2 = max(IMc2,[],3);

IMsk2 = squeeze(IMsk(:,:,2,:)); %we use channel 2
IMsk2 = IMsk2 - imgaussfilt(IMsk2, 29, 'padding', 'symmetric'); 
IMsk2 = imgaussfilt3(IMsk2, [4 4 0.5]);
IMsk2 = max(IMsk2,[],3);

%convert to 16 bit
IM = max(0, IM- prctile(IM(:), 1));

%convert correlation image 
IMc2 = max(0, IMc2- prctile(IMc2(:), 1));
IMc2 = IMc2./prctile(IMc2(:), 99.9);
IMc2(~isfinite(IMc2)) = 0;
IMc2 = uint16(IMc2.*65000);
IMsk2 = max(0, IMsk2- prctile(IMsk2(:), 1));
IMsk2 = IMsk2./max(IMsk2, [], 'all');
IMsk2(~isfinite(IMsk2)) = 0;
IMsk2 = uint16(IMsk2.*65000);

%remap refIM to full dynamic range; off by default
normalizeIntensity = false;
if normalizeIntensity
    IM = IM./max(IM, [], 'all'); %./max(IMout, 'all')
    IM = uint16(IM.*65000); 
else
    IM = uint16(IM);
end

%save reference image, one stack for each channel
disp('Saving...')

outputPathCh1 = [tiffFn{1}(1:end-4) '-REF_Ch1.ome.tif'];
outputPathCh2 = [tiffFn{1}(1:end-4) '-REF_Ch2.ome.tif'];
outputPathCh1corr = [tiffFn{1}(1:end-4) '-CORR_Ch1.ome.tif'];
outputPathCh2corr = [tiffFn{1}(1:end-4) '-CORR_Ch2.ome.tif'];
outputPathCh1sk = [tiffFn{1}(1:end-4) '-SK_Ch1.ome.tif'];
outputPathCh2sk = [tiffFn{1}(1:end-4) '-SK_Ch2.ome.tif'];
outputPathCh2corrSummary = [tiffFn{1}(1:end-4) '-CORR_2D_Ch2.ome.tif'];

bfCheckJavaPath;
metadata1 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,1,:)));
metadata2 = createMinimalOMEXMLMetadata(squeeze(IM(:,:,2,:)));
pixelSizeObj = ome.units.quantity.Length(java.lang.Double(pixelSizeUm), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeX(pixelSizeObj,0); metadata2.setPixelsPhysicalSizeX(pixelSizeObj, 0);
metadata1.setPixelsPhysicalSizeY(pixelSizeObj, 0); metadata2.setPixelsPhysicalSizeY(pixelSizeObj, 0);
pixelSizeZObj = ome.units.quantity.Length(java.lang.Double(abs(zStep)), ome.units.UNITS.MICROMETER);
metadata1.setPixelsPhysicalSizeZ(pixelSizeZObj, 0); metadata2.setPixelsPhysicalSizeZ(pixelSizeZObj, 0);
bfsave(squeeze(IM(:,:,1,:)), outputPathCh1, 'BigTiff', true, 'metadata', metadata1);
bfsave(squeeze(IM(:,:,2,:)), outputPathCh2, 'BigTiff', true, 'metadata', metadata2);
bfsave(squeeze(IMc(:,:,1,:)), outputPathCh1corr, 'BigTiff', true);
bfsave(squeeze(IMc(:,:,2,:)), outputPathCh2corr, 'BigTiff', true);
bfsave(squeeze(IMsk(:,:,1,:)), outputPathCh1sk, 'BigTiff', true);
bfsave(squeeze(IMsk(:,:,2,:)), outputPathCh2sk, 'BigTiff', true);
bfsave(IMc2, outputPathCh2corrSummary);

function [refPlane, corrPlane, skewPlane] = alignMultiChannel(IMs, b,a)
Y = squeeze(makeHighPass(sum(IMs,3))); % data order is XYCTZV
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',20,'us_fac',4,'init_batch',40, 'correct_bidir', true);
refPlane = [];
[~,shifts,~,options, col_shift] = normcorre(Y,options_rigid);
for chIx= size(IMs,3):-1:1
    tmp= apply_shifts(squeeze(IMs(:,:,chIx,:)),shifts,options,0,0,0,col_shift);
   
    refPlane(:,:,chIx) = mean(tmp,3);

    %compute correlation image
    %downsample in space
    nDS = 1;
    for DSiter = 1:nDS
        tmp = tmp(1:2:2*floor(end/2), 1:2:2*floor(end/2),:) + tmp(2:2:2*floor(end/2), 1:2:2*floor(end/2),:) + tmp(1:2:2*floor(end/2), 2:2:2*floor(end/2),:) + tmp(2:2:2*floor(end/2), 2:2:2*floor(end/2),:);
    end
    %highpass in time
    HP = permute(filtfilt(b,a,permute(double(tmp), [3 1 2])), [2 3 1]);
    ss = sum(HP.^2,3, 'omitnan');
    vertC = sum(HP .* circshift(HP, [1 0 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0 ]));
    horzC = sum(HP .* circshift(HP, [0 1 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0 ]));
    C = mean(cat(3, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),3, 'omitnan');
    C(isnan(C)) = 0;
% 
    %compute skewness image
%     IMgamma = mean(tmp,3); %average downsampled image
%     IMgamma = max(0, IMgamma - prctile(IMgamma, 1,'all')); %estimate baseline
%     IMgamma = sqrt(IMgamma);
    sk = skewness(HP,0,3);
    
    corrPlane(:,:,chIx) = C;
    skewPlane(:,:,chIx) = sk;
end
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
