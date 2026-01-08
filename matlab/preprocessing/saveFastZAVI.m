function saveFastZAVI(fullPathToTrialTable)

if ~nargin
    [fn, dr] = uigetfile('*trialTable*.mat');
else
    [dr, fn, ext] = fileparts(fullPathToTrialTable); fn = [fn ext]; 
end

%PARAMETER SETTING
nWorkers = 6;

%load the trial Table, which sets correspondences between the two DMDs
load([dr filesep fn], 'trialTable');

%% set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 0;
else
    poolsize = p.NumWorkers;
end
nWorkers = min([nWorkers, numel(trialTable.filename), feature('numcores')]);
if poolsize~=nWorkers
    delete(gcp('nocreate'));
    if nWorkers<15
        warning('You are using few parallel workers!');
    end
    disp(['Parallel workers:' int2str(nWorkers)])
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end

%% get superpixel data
lookupTable = getSuperPixelData(dr,trialTable);

%% conduct alignment in parallel
nDMDs = size(trialTable.filename,1);

[dixs,fixs] = ndgrid(1:nDMDs,1:length(trialTable.trueTrialIx));
parfor p_ix = 1:numel(fixs)
    f_ix = fixs(p_ix); DMD_ix = dixs(p_ix);
    saveFastZTifAsync(dr, trialTable, lookupTable, f_ix, DMD_ix);
end

disp('done saving videos.')
end

function spData = getSuperPixelData(dr, trialTable)
nDMDs = size(trialTable.filename,1);
for DMDix = nDMDs:-1:1
    [~,n] = fileparts(trialTable.filename{DMDix,1});
    n_base = regexprep(n,'-TRIAL[0-9]+','','ignorecase');
    n_base = regexprep(n_base,'-CYCLE[0-9]+','','ignorecase');
    metaDataFileName = fullfile(dr, [n_base '.meta']);
    mustBeFile(metaDataFileName);
    metaData = load(metaDataFileName, '-mat');

    list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '_CONFIG2-REFERENCE*']);
    if isempty(list)
        list = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*-REFERENCE*']);
    end
    ReferenceStack_ = slap2.gui.refstack.ReferenceStack.loadTif(fullfile(list(1).folder, list(1).name));
    numChannelsRefStack = numel(ReferenceStack_.data);
    for chIx = numChannelsRefStack:-1:1
        refStack(:,:,:,chIx) = permute(ReferenceStack_.data{chIx},[2 1 3]);
    end

    if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
        numLinesPerCycle = length(metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs);

        zs_ix = horzcat(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{:});
        zs_ix = unique(zs_ix);
        zs = zeros(length(zs_ix),1);

        true_zs = cellfun(@(x) x.z, metaData.AcquisitionContainer.ROIs.rois);

        zPlanes = nan(1,numLinesPerCycle);
        for ix = 1:numLinesPerCycle
            if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix})
                zPlanes(ix) = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{ix}(1);
            end
        end
        for ix = 1:length(zs_ix)
            if all(isnan(metaData.AcquisitionContainer.AcquisitionPlan.zTrajectory))
                zs(ix) = metaData.remoteFocusPosition_um;
            else
                zPlane_um = mean(metaData.AcquisitionContainer.AcquisitionPlan.zTrajectory(zPlanes == zs_ix(ix)));
                [~, zPlane_trueIx] = min(abs(true_zs - zPlane_um));
                zs(ix) = true_zs(zPlane_trueIx);
            end
        end
    else
        zs = metaData.AcquisitionContainer.ParsePlan.zs;
        numLinesPerCycle = length(metaData.AcquisitionContainer.ParsePlan.acqParsePlan);
    end

    dmdPixelsPerColumn = metaData.dmdPixelsPerColumn;
    dmdPixelsPerRow = metaData.dmdPixelsPerRow;
    numFastZs = length(zs);

    % get list of superpixels and extract data

    allSuperPixelIDs{DMDix} = [];

    fastZ2RefZ{DMDix} = zeros(numFastZs,1);
    middleRefZ = (max(ReferenceStack_.zs) + min(ReferenceStack_.zs))/2;
    middleZ = (max(zs) + min(zs))/2;
    for z = 1:numFastZs
        [~, ind] = min(abs((ReferenceStack_.zs - middleRefZ) - (zs(z) - middleZ)));
        fastZ2RefZ{DMDix}(z) = ind;
    end

    for lineSweepIdx = 1:numLinesPerCycle
        if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
            superPixIdxs = metaData.AcquisitionContainer.AcquisitionPlan.superPixelIDs{lineSweepIdx}';
        else
            superPixIdxs = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).superPixelID;
        end

        if numel(superPixIdxs) == 0; continue; end
        
        if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
            zIdx = metaData.AcquisitionContainer.AcquisitionPlan.activeZs{lineSweepIdx}(1);
        else
            zIdx = metaData.AcquisitionContainer.ParsePlan.acqParsePlan(lineSweepIdx).sliceIdx(1)+1;
        end

        spIDs = superPixIdxs*100+zIdx; % add z plane to end of superpixel ID
        allSuperPixelIDs{DMDix} = [allSuperPixelIDs{DMDix}; spIDs]; % make list of unique superpixels across all Zs
    end

    [allSuperPixelIDs{DMDix}, ~, ic] = unique(allSuperPixelIDs{DMDix});
    spSampleCt = accumarray(ic, 1);

    fprintf("%d superpixels detected\n", length(allSuperPixelIDs{DMDix}));

    % make sparse matrix with each superpixel's corresponding mask (roiMasks)

    fprintf("Calculating ROI Masks... ");
    tic;

    % using sparse matrix
    sparseMaskInds{DMDix} = [];
    if ~isempty(metaData.AcquisitionContainer.AcquisitionPlan)
        allPixelReplacementMaps = metaData.AcquisitionContainer.AcquisitionPlan.pixelReplacementMaps;
    else
        allPixelReplacementMaps = metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
    end

    for i = 1:length(allSuperPixelIDs{DMDix})
        tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);

        sp = allSuperPixelIDs{DMDix}(i);
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

        sparseMaskInds{DMDix} = [sparseMaskInds{DMDix}; openPixs, ones(size(openPixs))*i];
    end
    clear('tmpMask');

    roiMasks = sparse(sparseMaskInds{DMDix}(:,1),sparseMaskInds{DMDix}(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs{DMDix}));
    fprintf('done DMD%d - took %f sec\n', DMDix, toc);
end

spData.allSuperPixelIDs = allSuperPixelIDs;
spData.sparseMaskInds = sparseMaskInds;
spData.fastZ2RefZ = fastZ2RefZ;

end


function fnwrite = saveFastZTifAsync(dr, trialTable, lookupTable, f_ix, DMD_ix)

if nargin < 6
    interpMov = true;
end

fn = trialTable.filename{DMD_ix,f_ix};
fnW = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD' int2str(DMD_ix)];
firstLine = trialTable.firstLine(DMD_ix,f_ix);
lastLine = trialTable.lastLine(DMD_ix, f_ix);

disp(['Processing: ' [dr filesep fn]])

hSlap2DataFile = slap2.Slap2DataFile([dr filesep fn]);
if isprop(hSlap2DataFile, 'hDataFile')
    hLowLevelDataFile = hSlap2DataFile.hDataFile;
else
    hLowLevelDataFile = hSlap2DataFile.hMultiDataFiles;
end
numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
totalCycles = hLowLevelDataFile.numCycles;
numChannels = hSlap2DataFile.numChannels;
channels = hLowLevelDataFile.metaData.channelsSave;

assert(length(channels) == numChannels, 'Saved channels does not match numChannels!');

linerateHz = 1/hLowLevelDataFile.metaData.linePeriod_s;
alignHz = linerateHz / numLinesPerCycle / 2;
dt = linerateHz/alignHz;

if ~exist([dr filesep 'qc_videos'],'dir')
    mkdir([dr filesep 'qc_videos'])
end

fnwrite = [dr filesep 'qc_videos' filesep fnW '_FASTZ-' int2str(alignHz) 'Hz.avi'];

fnwriteTmp = ['C:\temp' filesep int2str(round(rand(1)*10000)) '.avi'];

if numChannels==2
    fnwrite = [dr filesep 'qc_videos' filesep fnW '_FASTZ_CH1-' int2str(alignHz) 'Hz.avi'];
    fnwrite2 = [dr filesep 'qc_videos' filesep fnW '_FASTZ_CH2-' int2str(alignHz) 'Hz.avi'];
    fnwriteTmp2 = ['C:\temp' filesep int2str(round(rand(1)*10000)) '.avi'];
end

if exist(fnwrite, 'file')
    disp([fn ' is already saved as tif; skipping']);
    return
end

DSframes = ceil(firstLine:dt:lastLine);
nDSframes= length(DSframes); %number of downsampled frames

%% load metadata
dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
numFastZs = length(hLowLevelDataFile.fastZs);

%% prep tif

% pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
% fTIF = Fast_BigTiff_Write(fnwriteTmp,pixelscale,0);

video = VideoWriter(fnwriteTmp,'Motion JPEG AVI');
video.FrameRate = alignHz;

if numChannels==2
    video2 = VideoWriter(fnwriteTmp2,'Motion JPEG AVI');
    video2.FrameRate = alignHz;
end

% bandHeight = 9;

% Get unique superpixel IDs and their group memberships
[~, ~, groupIndices] = unique(lookupTable.sparseMaskInds{DMD_ix}(:, 2));

% Group the pixel indices by superpixel
splitPixels = accumarray(groupIndices, lookupTable.sparseMaskInds{DMD_ix}(:, 1), [], @(x) {x});

% Compute the median pixel index for each superpixel
medianIndices = cellfun(@(x) x(round(length(x)/2)), splitPixels);

% Convert linear indices to row and column coordinates
[spRows, spCols, spPlanes] = ind2sub([dmdPixelsPerColumn, dmdPixelsPerRow, numFastZs], medianIndices);

zcolors = turbo(max(lookupTable.fastZ2RefZ{DMD_ix})-min(lookupTable.fastZ2RefZ{DMD_ix})+1);

planeBorder = single(zeros(dmdPixelsPerColumn,dmdPixelsPerRow,3));
for planeIx = 1:numFastZs
    tmp = false(dmdPixelsPerColumn,dmdPixelsPerRow);
    activeInds = sub2ind([dmdPixelsPerColumn,dmdPixelsPerRow],spRows(spPlanes == planeIx),spCols(spPlanes == planeIx));
    tmp(activeInds) = true;
    planeMask = imerode(imdilate(tmp,ones(9,3),'same'),ones(1,3),'same');
    borderMask = single(imdilate(planeMask,ones(5,5),'same') & ~planeMask);
    borderMask = borderMask .* reshape(zcolors(lookupTable.fastZ2RefZ{DMD_ix}(planeIx)-min(lookupTable.fastZ2RefZ{DMD_ix})+1,:),[1 1 3]);
    planeBorder(borderMask > 0) = borderMask(borderMask > 0);
    clear tmp;
end

saveFailed = false;
fprintf("Saving each frame...\n")
tic

img_lims = [0 1];
open(video);
if numChannels==2
    img_lims2 = [0 1];
    open(video2);
end
for DSframeIx = 1:nDSframes
    timeWindow = max(1,floor(DSframes(DSframeIx)-2*dt)):min(ceil(DSframes(DSframeIx)+2*dt),hLowLevelDataFile.totalNumLines);

    if ~mod(DSframeIx, 10)
        disp([int2str(DSframeIx) ' of ' int2str(nDSframes) ' - ' num2str(toc/(DSframeIx-1)) ' sec/frame']);
    end

    lineIndices  = mod(timeWindow-1,numLinesPerCycle)+1;
    cycleIndices = floor((timeWindow-1) / numLinesPerCycle)+1;

    weights = exp(-abs((DSframes(DSframeIx) - timeWindow) / dt));

    % load data
    data = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),numChannels);
    spCt = zeros(length(lookupTable.allSuperPixelIDs{DMD_ix}),1);
    allLineData = hLowLevelDataFile.getLineData(lineIndices, cycleIndices, channels);
    for t = 1:length(allLineData)
        superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineIndices(t)};

        if numel(superPixIdxs) == 0; continue; end

        lineData = allLineData{t};
        badSPsMask = any(lineData < 0,2);
        
        if sum(badSPsMask) > 0
            superPixIdxs(badSPsMask) = [];
            lineData(find(badSPsMask),:) = [];
        end

        zIdx = hLowLevelDataFile.lineFastZIdxs(lineIndices(t));

        spID = superPixIdxs*100 + uint32(zIdx); % make superpixel index with Z plane
        [~,spIdx] = ismember(spID,lookupTable.allSuperPixelIDs{DMD_ix});
        
        weight = weights(t);
        data(spIdx(spIdx>0),:) = data(spIdx(spIdx>0),:) + single(lineData(spIdx>0,:) .* weight);
        spCt(spIdx(spIdx>0)) = spCt(spIdx(spIdx>0)) + weight;
    end
    data = data ./ spCt;
    data(spCt == 0,:) = nan;

    if mean(spCt == 0) > 0.5; continue; end

    A1 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
    A2 = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
    for zIdx = 1:numFastZs
        spIdxs = find(spPlanes == zIdx);
        tmp = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
        tmp(sub2ind(size(tmp),spRows(spIdxs),spCols(spIdxs))) = data(spIdxs,1);
        tmp = imdilate(tmp,ones(9,1));
        A1(:,:) = max(cat(3,A1,tmp),[],3,'omitmissing');
        if numChannels==2
            tmp = nan(dmdPixelsPerColumn,dmdPixelsPerRow);
            tmp(sub2ind(size(tmp),spRows(spIdxs),spCols(spIdxs))) = data(spIdxs,2);
            tmp = imdilate(tmp,ones(9,1));
            A2(:,:) = max(cat(3,A2,tmp),[],3,'omitmissing');
        end
    end
    if DSframeIx == 1
        img_lims(1) = 0;
        img_lims(2) = prctile(A1(:),99.9);
    end
    A1 = repmat(clip(A1,img_lims(1),img_lims(2)) / img_lims(2),[1 1 3]);
    A1(planeBorder > 0) = planeBorder(planeBorder > 0);
    writeVideo(video,A1);
    % fTIF.WriteIMG(uint16(A1));
    if numChannels==2
        if DSframeIx == 1
            img_lims2(1) = 0;
            img_lims2(2) = prctile(A2(:),99.9);
        end
        A2 = repmat(clip(A2,img_lims2(1),img_lims2(2)) / img_lims2(2),[1 1 3]);
        A2(planeBorder > 0) = planeBorder(planeBorder > 0);
        writeVideo(video2,A2);
        % fTIF.WriteIMG(uint16(A2));
    end
end

close(video);
if numChannels==2
    close(video2);
end

if saveFailed
    disp(['SAVE ERROR OCCURRED FOR FILE: ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
    return
end

copyfile(fnwriteTmp,fnwrite);
% delete(fnwriteTmp);
if numChannels==2
    copyfile(fnwriteTmp2,fnwrite2);
end

disp(['saving done - took ' num2str(toc) ' sec'])
end