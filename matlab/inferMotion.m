function inferMotion(varargin)

% get data file
if numel(varargin) < 1 || isempty(varargin{1})
    [fns_data, dr_data] = uigetfile({'*.meta'}, 'Select Data File', 'multiselect', 'on');
    dataFile = [dr_data fns_data];
else
    dataFile = varargin{1};
end

% load reference stack
if numel(varargin) < 2 || isempty(varargin{2})
    [fns, dr] = uigetfile('*.tif', 'Select Reference Stack', 'multiselect', 'on');
    try
        ReferenceStack_ = slap2.util.ReferenceStack.loadTif([dr fns]);
        refStack = permute(ReferenceStack_.data,[2 1 3]);
    catch
        refStack = tiffreadVolume([dr fns]);
    end

    if strcmpi(fns(end-6:end-5),'CH')
        channel = str2num(fns(end-4));
    else
        channel = 1;
    end
else
    refStackFlag = true;
    refStackDir = varargin{2};
    refStack = tiffreadVolume(refStackDir);
    if strcmpi(refStackDir(end-6:end-5),'CH')
        channel = str2num(refStackDir(end-4));
    else
        channel = 1;
    end
end 

% read in downsampling factor (default to 7x temporal ds)
if numel(varargin) < 3
    % ds = 7;
    ds = 3;
else
    ds = varargin{3};
end

% read in whether robust metrics should be used
if numel(varargin) < 4
    robust = false;
    b = 3;
else
    robust = varargin{4};
    b = varargin{5};
end

%% load SLAP2 data and metadata

hDataViewer = slap2.Slap2DataViewer(dataFile);
hLowLevelDataFile = hDataViewer.hDataFile.hDataFile;
fname = hLowLevelDataFile.filename;

numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
totalCycles = hLowLevelDataFile.numCycles;
numCycles = floor(totalCycles / 2^ds);

dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
numFastZs = length(hLowLevelDataFile.fastZs);


%% get list of superpixels and extract data

allSuperPixelIDs = [];

fastZ2RefZ = zeros(numFastZs,1);
for z = 1:numFastZs
    [~, ind] = min(abs(ReferenceStack_.zs - hLowLevelDataFile.fastZs(z)));
    fastZ2RefZ(z) = ind;
end

for lineSweepIdx = 1:numLinesPerCycle
    superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};

    if size(superPixIdxs,1) == 0; continue; end

    zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
    
    spIDs = superPixIdxs*100+zIdx; % add z plane to end of superpixel ID
    allSuperPixelIDs = [allSuperPixelIDs; spIDs]; % make list of unique superpixels across all Zs
end

[allSuperPixelIDs, ~, ic] = unique(allSuperPixelIDs);
spSampleCt = accumarray(ic, 1);

fprintf("%d superpixels detected\n", length(allSuperPixelIDs));

%% make sparse matrix with each superpixel's corresponding mask (roiMasks)

fprintf("Calculating ROI Masks... ");
tic;

% using sparse matrix
sparseMaskInds = [];
allPixelReplacementMaps = hLowLevelDataFile.metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;

for i = 1:length(allSuperPixelIDs)
    tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);

    sp = allSuperPixelIDs(i);
    zIdx = rem(sp,100);
    superPixIdx = (sp - zIdx) / 100;
    pixelReplacementMap = allPixelReplacementMaps{zIdx};

    open = double(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;
    openR = floor((open-1) /dmdPixelsPerRow)+1;
    openC = open - (openR-1) * dmdPixelsPerRow;
    openPixs = uint32(openR + (openC-1) * dmdPixelsPerColumn + double(zIdx-1) * dmdPixelsPerColumn * dmdPixelsPerRow);

    sparseMaskInds = [sparseMaskInds; openPixs, ones(size(openPixs))*i];
end
clear('tmpMask');

roiMasks = sparse(sparseMaskInds(:,1),sparseMaskInds(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs));
fprintf('done - took %f sec\n', toc);

%% make sure refStack is larger than imaged frame

bl = single(refStack)/100; %  reference image (X x Y x Z)
bl(bl < 0) = 0; % remove any negative values, should not happen
bl = padarray(bl,[floor((max(size(bl,1),dmdPixelsPerColumn)-size(bl,1))/2),...
                  floor((max(size(bl,2),dmdPixelsPerRow)-size(bl,2))/2),...
                  floor((max(size(bl,3),numFastZs)-size(bl,3))/2)],...
              mean(bl(:)),...
              'both');

if size(bl,1) < dmdPixelsPerColumn
    bl = padarray(bl,[1,0,0],mean(bl(:)),'pre');
end

if size(bl,2) < dmdPixelsPerRow
    bl = padarray(bl,[0,1,0],mean(bl(:)),'pre');
end

%% Calculate lookup table

xPre = 25; xPost = 25;
yPre = 25; yPost = 25;

zPre = 10;
zPost = 10;

xMotRange = xPre + xPost + 1;
yMotRange = yPre + yPost + 1;
zMotRange = zPre + zPost + 1;

% query user to select a existing lookup table
[lookupFile, lookupDir] = uigetfile('*.mat', 'Select Lookup Table');

% make lookup table if it doesn't exist
if sum(lookupFile == 0)

    likelihood_means = makeLookupTable(bl, sparseMaskInds, numFastZs, fastZ2RefZ,-yPre:yPost,-xPre:xPost,-zPre:zPost);
    
    if exist('fname')
        save([fname(1:end-5) '_LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
    else
        save([filepath '\LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
    end
else % load user selected lookup table
    fprintf("Loading lookup table... ")
    load([lookupDir lookupFile]);
    fprintf('done\n');
end

%%  Calculate log likelihoods and infer motion from MAP
log_means = log(likelihood_means + 1e-8);

% how much change is allowed in each dimension at each step
searchRadius = 4;

dsMotion = zeros(numCycles,3);
dsBrightness = zeros(numCycles,1);
dataMatrix = zeros(length(allSuperPixelIDs),numCycles);
expectedMatrix = zeros(length(allSuperPixelIDs),numCycles);

fprintf("Inferring motion... ")
tic

ySearch = 1:yMotRange;
xSearch = 1:xMotRange;
zSearch = 1:zMotRange;

for cycleIdx = 1:numCycles
    % load data
    data = zeros(length(allSuperPixelIDs),1);
    for c = 1:2^ds
        allLineData = hLowLevelDataFile.getLineData(1:numLinesPerCycle, ((cycleIdx-1)*2^ds+c)*ones(numLinesPerCycle,1));
        for lineSweepIdx = 1:numLinesPerCycle
    
            superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};
    
            if size(superPixIdxs,1) == 0; continue; end
    
            lineData = allLineData{lineSweepIdx};
            zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
    
            spID = superPixIdxs*100 + uint32(zIdx); % make superpixel index with Z plane
            [~,spIdx] = ismember(spID,allSuperPixelIDs);
    
            data(spIdx(spIdx>0)) = data(spIdx(spIdx>0)) + single(lineData(spIdx>0,channel));
        end
    end
    data = data ./ spSampleCt ./ 100;

    % calculate log likelihoods at all shifts
    [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, likelihood_means,log_means,ySearch,xSearch,zSearch,robust);

    [LL, I] = max(logLikelihoodTable(:));

    [My, Mx, Mz] = ind2sub(size(logLikelihoodTable),I);
    dsMotion(cycleIdx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)];

    dsBrightness(cycleIdx) = scalingFactorTable(My, Mx, Mz);
    dataMatrix(:,cycleIdx) = data;
    expectedMatrix(:,cycleIdx) = dsBrightness(cycleIdx) .* likelihood_means(ySearch(My), xSearch(Mx), zSearch(Mz),:);

    ySearch = max(1,dsMotion(cycleIdx,1) - searchRadius):min(xMotRange,dsMotion(cycleIdx,1) + searchRadius);
    xSearch = max(1,dsMotion(cycleIdx,2) - searchRadius):min(yMotRange,dsMotion(cycleIdx,2) + searchRadius);
    zSearch = max(1,dsMotion(cycleIdx,3) - searchRadius):min(zMotRange,dsMotion(cycleIdx,3) + searchRadius);
end

disp(['done - took ' num2str(toc/numCycles) ' sec per frame'])

%% Upsample if downsampled
frames = 1:totalCycles;
dsFrames = frames(1:2^ds:(2^ds*numCycles));

motion = interp1(dsFrames,dsMotion,frames);
motion = motion - [yPre+1 xPre+1 zPre+1];

brightness = interp1(dsFrames,dsBrightness,frames,'linear','extrap');

%% Plot results

figure(401); clf;
ax1 = subplot(3,1,1);
plot(motion(:,1)); hold on;
title('y motion')
legend('inferred','true','Location','best')

ax2 = subplot(3,1,2);
plot(motion(:,2)); hold on;
title('x motion')
legend('inferred','true','Location','best')

ax3 = subplot(3,1,3);
plot(motion(:,3)); hold on;
title('z motion')
legend('inferred','true','Location','best')

linkaxes([ax1,ax2])
linkaxes([ax1,ax2,ax3],'x')

[filepath, ~, ~] = fileparts(dataFile);
filepath = convertStringsToChars(filepath);

saveas(gcf,[filepath '\inferredMotion.fig']);

%% Save out data

inferMotionOut.motion = motion;
inferMotionOut.brightness = brightness';
inferMotionOut.dataMatrix = dataMatrix;
inferMotionOut.expectedMatrix = expectedMatrix;
inferMotionOut.sparseMaskInds = sparseMaskInds;

if robust
    save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_ROBUST_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
else
    save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
end