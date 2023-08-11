function alignBergamoStack(tiffFn)
%Aligns a scanimage slow stack acquired on the Bergamo
%uses the ScanImageTiffReader

if nargin<1 || isempty(tiffFn)
    [fn, dr] = uigetfile('*.tif*');
    if isnumeric(fn)
        return % user abort
    end
    tiffFn = fullfile(dr,fn);
end
[IMs,ImDescriptions] = readTiffNative(tiffFn);
numChannels = 2;
numFramesPerSlice = 100;
numZs = size(IMs,3)./numChannels./numFramesPerSlice;

framesPerCycle = numChannels*numFramesPerSlice*numZs;
nVol = size(IMs,3)/framesPerCycle;
if nVol ~= floor(nVol)
    warning('The number of saved TIF images is not a multiple of the stack size!')
end
nVol = floor(nVol);
totalNumFrames = nVol * framesPerCycle;
if totalNumFrames<size(IMs,3)
    IMs = IMs(:,:,1:totalNumFrames);
end

% data order is XYCTZV
IMs = permute(single(IMs), [2 1 3]);

%detect  and remove bidi artifact
IMs = reshape(IMs,size(IMs,1),size(IMs,2),numChannels,numFramesPerSlice,numZs,[]);

% align
refIM = alignMultiChannelbyPlane(IMs);

%save reference image, one stack for each channel
disp('Saving...')

metadata = struct;
metadata.zSpacing = 3; %um
metadata.pixelSize = 0.5; %um
for Ch = 1:numChannels
    metadata.channel = Ch;
    IMsave = squeeze(refIM(:,:,Ch,1,:,1));
    saveTif(IMsave, [tiffFn(1:end-4) '_REF_Ch' int2str(Ch) '.tif'], metadata)
end

keyboard

    function refIM = alignMultiChannelbyPlane(IMs)
        nVol = size(IMs,6);
        nZ = size(IMs, 5); % data order is XYCTZV
        [p,n] = fileparts(tiffFn);

        msg = sprintf(['Aligning file ' n '\n#Stacks: ' int2str(nVol) ' # planes: ' int2str(nZ)]);
        msg = most.idioms.latexEscape(msg);
        disp(msg)

        Y = makeHighPass(sum(IMs,3)); % data order is XYCTZV

        options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',15,'us_fac',4,'init_batch',30, 'correct_bidir', true);
        refIM = [];
        for Z = nZ:-1:1
            [~,shifts,~,options, col_shift] = normcorre(squeeze(Y(:,:,1,:,Z)),options_rigid);
            for chIx= 1:size(IMs,3)
                tmp= apply_shifts(squeeze(IMs(:,:,chIx,:,Z)),shifts,options,0,0,0,col_shift);
                refIM(:,:,chIx,1,Z) = mean(tmp,3);
            end
        end
    end

end
function IMHP = makeHighPass(IM1)
IMHP = IM1; %sum the channels to speed alignment
IMHP(IMHP<0) = 0;
IMHP(isnan(IMHP)) = 0;
IMHP = imgaussfilt(IMHP, [0.5 0.5])-imgaussfilt(IMHP, [4 4]);
end

function [IMs,ImDescriptions] = readTiff(tiffFn)
hTiff= ScanImageTiffReader(tiffFn);
try
    ImDescriptions = hTiff.descriptions;
    IMs = hTiff.data;
catch ME
    most.idioms.safeDeleteObj(hTiff);
    ME.rethrow();
end
most.idioms.safeDeleteObj(hTiff);
end

function [IMs,ImDescriptions] = readTiffNative(tiffFn)
IMs = {};
ImDescriptions = {};
w = warning();
warning('off','imageio:tiffmexutils:libtiffWarning');
hTiff = Tiff(tiffFn);
try
    while true
        IMs{end+1} = hTiff.read()';
        ImDescriptions{end+1} = hTiff.getTag('ImageDescription');
        try
            hTiff.nextDirectory();
        catch ME
            break;
        end
    end
catch ME
    warning(w);
    hTiff.delete();
    ME.rethrow();
end
warning(w);
hTiff.delete();
IMs = cat(3,IMs{:});
end

function saveTif(IM,fileName, metadata)
[filepath,fileName,ext] = fileparts(fileName);

if exist(filepath,'dir')
    fileName = fullfile(filepath,[fileName,ext]);
else
    fileName = [fileName,ext];
end

[fileName,filePath] = uiputfile('.tif','Save ReferenceStack',fileName);
if isnumeric(fileName)
    return %User abort
end
fileName = fullfile(filePath,fileName);

imageDescriptions = {};
for zIdx = 1:size(IM,3)
    metadata.ReferenceStackFileVersion = 1;
    metadata.z= zIdx;
    imageDescriptions{zIdx} = jsonencode(metadata);
end

isTransposed = true;
most.util.writeTiff(fileName, IM, isTransposed, imageDescriptions);
end
