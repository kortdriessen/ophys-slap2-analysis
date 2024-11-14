function integrationRegistration(varargin)

% get data file
if numel(varargin) < 1 || isempty(varargin{1})
    [fns, dr] = uigetfile({'*.meta'}, 'Select Data File', 'multiselect', 'on');
    dataFile = [dr fns];
else
    dataFile = varargin{1};
end

% load reference stack
if numel(varargin) < 2 || isempty(varargin{2})
    [fns_ref, dr_ref] = uigetfile('*.tif', 'Select Reference Stack', 'multiselect', 'on');
    try
        ReferenceStack_ = slap2.gui.refstack.ReferenceStack.loadTif([dr_ref fns_ref]);
        numChannelsRefStack = numel(ReferenceStack_.data);
        for chIx = numChannelsRefStack:-1:1
            refStack(:,:,:,chIx) = permute(ReferenceStack_.data{chIx},[2 1 3]);
        end
    catch
        refStack = tiffreadVolume([dr_ref fns_ref]);
    end
else
    refStackFlag = true;
    refStackDir = varargin{2};
    refStack = tiffreadVolume(refStackDir);
end 

% % read in downsampling factor (default to 7x temporal ds)
% if numel(varargin) < 3
%     % ds = 7;
%     ds = 3;
% else
%     ds = varargin{3};
% end
% 
% % read in whether robust metrics should be used
% if numel(varargin) < 4
%     robust = false;
%     b = 3;
% else
%     robust = varargin{4};
%     b = varargin{5};
% end

%%
%load the trial Table, which sets correspondences between the two DMDs
load([dr filesep 'trialTable.mat'], 'trialTable');
alignHz = 80; %we will align data at this timescale, Hz

aData.alignHz = alignHz;
aData.maxshift = 50;
aData.clipShift = 5;%the maximum allowable shift per frame
aData.alpha = 0.005; %exponential time constant for template

%set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
nWorkers = 24;
if poolsize~=nWorkers
    delete(gcp('nocreate'));
    if nWorkers<15
        warning('You are using few parallel workers! adjust this in summarizeBCI.m');
    end
    disp(['Parallel workers:' int2str(nWorkers)])
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end

% make look up table first

[~,n] = fileparts(trialTable.DMD1filename{1});
n_base = regexprep(n,'-TRIAL[0-9]+$','','ignorecase');
metaDataFileName = fullfile(dr, [n_base '.meta']);
mustBeFile(metaDataFileName);
metaData = load(metaDataFileName, '-mat');

numLinesPerCycle = length(metaData.AcquisitionContainer.ParsePlan.acqParsePlan);

dmdPixelsPerColumn = metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = metaData.dmdPixelsPerRow;
numFastZs = length(metaData.AcquisitionContainer.ParsePlan.zs);

% get list of superpixels and extract data

allSuperPixelIDs = [];

fastZ2RefZ = zeros(numFastZs,1);
for z = 1:numFastZs
    [~, ind] = min(abs(ReferenceStack_.zs - metaData.AcquisitionContainer.ParsePlan.zs(z)));
    fastZ2RefZ(z) = ind;
end

for lineSweepIdx = 1:numLinesPerCycle
    superPixIdxs = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).superPixelID;

    if size(superPixIdxs,1) == 0; continue; end

    zIdx = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).sliceIdx(1) + 1;

    % only take integration mode pixels
    superPixIdxs(superPixIdxs <= dmdPixelsPerColumn*dmdPixelsPerRow) = [];

    spIDs = superPixIdxs*100+zIdx; % add z plane to end of superpixel ID
    allSuperPixelIDs = [allSuperPixelIDs; spIDs]; % make list of unique superpixels across all Zs
end

[allSuperPixelIDs, ~, ic] = unique(allSuperPixelIDs);
spSampleCt = accumarray(ic, 1);

fprintf("%d superpixels detected\n", length(allSuperPixelIDs));

% make sparse matrix with each superpixel's corresponding mask (roiMasks)

fprintf("Calculating ROI Masks... ");
tic;

% using sparse matrix
sparseMaskInds = [];
allPixelReplacementMaps = metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;

for i = 1:length(allSuperPixelIDs)
    tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);

    sp = allSuperPixelIDs(i);
    zIdx = rem(sp,100);
    superPixIdx = (sp - zIdx) / 100;
    pixelReplacementMap = allPixelReplacementMaps{zIdx};

    open = uint32(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;

    if isempty(open)
        open = superPixIdx+1;
    end

    openR = idivide(open-1, dmdPixelsPerRow, 'floor')+1;
    openC = open - (openR-1) * dmdPixelsPerRow;
    openPixs = uint32(openR + (openC-1) * dmdPixelsPerColumn + double(zIdx-1) * dmdPixelsPerColumn * dmdPixelsPerRow);

    sparseMaskInds = [sparseMaskInds; openPixs, ones(size(openPixs))*i];
end
clear('tmpMask');

roiMasks = sparse(sparseMaskInds(:,1),sparseMaskInds(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs));
fprintf('done - took %f sec\n', toc);

% make sure refStack is larger than imaged frame

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

% Calculate lookup table

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
    likelihood_means = likelihood_means .* repmat(reshape(spSampleCt,[1 1 1 length(allSuperPixelIDs)]),[size(likelihood_means,1:3) 1]);

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

%%
parfor f_ix = 1:length(trialTable.trueTrialIx)
    fnWrite = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD1'];
    alignIntegrationAsync(dr, trialTable.DMD1filename{f_ix}, fnWrite, trialTable.DMD1firstLine(f_ix), trialTable.DMD1lastLine(f_ix), aData, overwriteExisting);
    fnWrite = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD2'];
    alignIntegrationAsync(dr, trialTable.DMD2filename{f_ix}, fnWrite, trialTable.DMD2firstLine(f_ix), trialTable.DMD2lastLine(f_ix), aData, overwriteExisting);
end

disp('done integrationRegistration.')
end


function alignIntegrationAsync(dr, fn, fnW, firstLine, lastLine, aData, overwriteExisting)
disp(['Aligning: ' [dr filesep fn]])

fnwrite = [dr filesep fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
fnAdata = [dr filesep fnW '_ALIGNMENTDATA.mat'];

if ~overwriteExisting && exist(fnAdata, 'file') && exist(fnwrite, 'file')
    disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
    return
end

hSlap2DataFile = slap2.Slap2DataFile([dr filesep fn]);
hLowLevelDataFile = hSlap2DataFile.hDataFile;

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

    open = uint32(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;

    if isempty(open)
        open = superPixIdx+1;
    end

    openR = idivide(open-1, dmdPixelsPerRow, 'floor')+1;
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
    likelihood_means = likelihood_means .* repmat(reshape(spSampleCt,[1 1 1 length(allSuperPixelIDs)]),[size(likelihood_means,1:3) 1]);

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
    data = data ./ 100; %./ spSampleCt;

    % calculate log likelihoods at all shifts
    [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, likelihood_means,log_means,ySearch,xSearch,zSearch,robust);

    [LL, I] = max(logLikelihoodTable(:));

    [My, Mx, Mz] = ind2sub(size(logLikelihoodTable),I);

    if My>1 && My<size(logLikelihoodTable,1) && Mx>1 && Mx<size(logLikelihoodTable,2) && Mz>1 && Mz<size(logLikelihoodTable,3)
        %perform superresolution upsampling
        ratioY = min(1e6,(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My-1,Mx,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My+1,Mx,Mz)));
        dY = (1-ratioY)/(1+ratioY)/2;
        
        ratioX =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx-1,Mz))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx+1,Mz)));
        dX = (1-ratioX)/(1+ratioX)/2;

        ratioZ =min(1e6, (logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz-1))/(logLikelihoodTable(My,Mx,Mz) - logLikelihoodTable(My,Mx,Mz+1)));
        dZ = (1-ratioZ)/(1+ratioZ)/2;

        dsMotion(cycleIdx,:) = [ySearch(My)-dY; xSearch(Mx)-dX; zSearch(Mz)-dZ];
    
        % motion = shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
    else %the optimum is at an edge of search range; no superresolution
        dsMotion(cycleIdx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)];
    end

    dsBrightness(cycleIdx) = scalingFactorTable(My, Mx, Mz);
    dataMatrix(:,cycleIdx) = data;
    expectedMatrix(:,cycleIdx) = dsBrightness(cycleIdx) .* likelihood_means(ySearch(My), xSearch(Mx), zSearch(Mz),:);

    ySearch = max(1,round(dsMotion(cycleIdx,1)) - searchRadius):min(xMotRange,round(dsMotion(cycleIdx,1)) + searchRadius);
    xSearch = max(1,round(dsMotion(cycleIdx,2)) - searchRadius):min(yMotRange,round(dsMotion(cycleIdx,2)) + searchRadius);
    zSearch = max(1,round(dsMotion(cycleIdx,3)) - searchRadius):min(zMotRange,round(dsMotion(cycleIdx,3)) + searchRadius);
end

disp(['done - took ' num2str(toc/numCycles) ' sec per frame'])

%% Upsample if downsampled
frames = 1:totalCycles;
dsFrames = frames(1:2^ds:(2^ds*numCycles));

motion = interp1(dsFrames,dsMotion,frames,'linear','extrap');
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
inferMotionOut.fastZ2RefZ = fastZ2RefZ;

if robust
    save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_ROBUST_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
else
    save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
end
end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end