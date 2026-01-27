function [E, footprints]= loadAndProcessVoltage(dr, fn, firstLine, lastLine, aData, masks, params)

hSlap2DataFile = slap2.Slap2DataFile([dr filesep fn]);
if isprop(hSlap2DataFile, 'hDataFile')
    hLowLevelDataFile = hSlap2DataFile.hDataFile;
else
    hLowLevelDataFile = hSlap2DataFile.hMultiDataFiles;
end
numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
% totalCycles = hLowLevelDataFile.numCycles;
% totalLines = totalCycles*numLinesPerCycle;
numChannels = hSlap2DataFile.numChannels;
channels = hLowLevelDataFile.metaData.channelsSave;
assert(length(channels) == numChannels, 'Saved channels does not match numChannels!');
dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
dmdNumPix = dmdPixelsPerRow*dmdPixelsPerColumn;
%metaData = hLowLevelDataFile.metaData;

linerateHz = 1/hLowLevelDataFile.metaData.linePeriod_s;
dt = ceil(2*linerateHz/params.analyzeHz);
sampleTimes_s = (firstLine/linerateHz):(1/params.analyzeHz):(lastLine/linerateHz);
lineIxs = round(sampleTimes_s.*linerateHz);
nSamps = numel(lineIxs);

%determine discard frames
nDSframes = numel(aData.DSframes);
tmp = aData.recNegErr(:)- medfilt1(aData.recNegErr(:), round(4*aData.alignHz)); %smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
tmp = -tmp./(min(-0.005, prctile(tmp,5))); %normalize to the median-to-5th prctile interval, or 0.5% of intensity, whichever is larger
thresh = params.motionThresh; %decrease thresh to be more stringent on motion correction
window = 2*ceil(0.1*aData.alignHz)+1;% a window in time to censor aronud motion events, ~0.2 seconds;
discard = imclose(imdilate(tmp>thresh, ones(window,1)) | (tmp>(thresh/2) & imdilate(tmp>thresh, ones(2*window+1,1))), ones(window,1));
E.discardFrames = logical(interp1(linspace(firstLine, lastLine, nDSframes),single(discard), lineIxs, 'nearest','extrap'));  %resample discardFrames to the sample times

%upsample motion and recNegErr
fieldNames = {'onlineXshift', 'onlineYshift', 'onlineZshift', 'recNegErr', 'motionDSc', 'motionDSr', 'aError'};
for fieldIx = 1:length(fieldNames)
    fieldName = fieldNames{fieldIx};
    E.upsampledMotion.(fieldName) = interp1(linspace(firstLine, lastLine, nDSframes),aData.(fieldName), lineIxs);
end

%Get superpixel traces
hTrace = slap2.util.datafile.trace.Trace(hLowLevelDataFile,1,1); %false(dmdPixelsPerRow, dmdPixelsPerColumn)
hTrace.setPixelIdxs([], true(dmdPixelsPerRow, dmdPixelsPerColumn)); %px);
E.superPixelIDs = hTrace.superPixelIds;
nPx = numel(hTrace.TracePixels);
F = nan(nSamps, numel(hTrace.TracePixels));
weight = nan(nSamps, numel(hTrace.TracePixels));
TracePixels = hTrace.TracePixels;
parfor idx = 1:nPx
    [~, fullTrace, fullWeight] = TracePixels(idx).process(dt,numLinesPerCycle);
    F(:,idx) = fullTrace(lineIxs);
    weight(:,idx) = fullWeight(lineIxs);
end
E.F = F;
E.weight = weight;

%make the superpixel footprints
if nargout>1
footprints = nan(dmdPixelsPerRow,dmdPixelsPerColumn, nPx);
pixelReplacementMap = hLowLevelDataFile.zPixelReplacementMaps{1};
intMask = pixelReplacementMap(:,2) >= dmdNumPix;
%pixelReplacementMap(intMask,1) = pixelReplacementMap(intMask,1) + dmdNumPix;
map = nan(dmdPixelsPerRow,dmdPixelsPerColumn);
map(pixelReplacementMap(intMask,1)) = pixelReplacementMap(intMask,2);
for idx = 1:nPx
    %tmp = nan(dmdPixelsPerRow,dmdPixelsPerColumn, 'single');
    %tmp(map==E.superPixelIDs(idx)) = 1;
    footprints(:,:,idx) = (map==E.superPixelIDs(idx));
end
end

nans = isnan(F)|isnan(weight);
weight(nans)=0;
E.global.F = sum(F.*weight,2, 'omitnan')./sum(weight,2);
clear F weight

%extract each ROI
nROIs = size(masks,3);
expecteddt = ceil(params.baselineWindow_s*linerateHz); %for single-pixel traces, just needs to be long enough to span a cycle
roiF = nan(nSamps, nROIs);
roiW = nan(nSamps, nROIs);
parfor idx = 1:nROIs
    hTrace = slap2.util.datafile.trace.Trace(hLowLevelDataFile,1,1);
    rasterPixels = masks(:,:,idx);
    integrationPixels = masks(:,:,idx);
    hTrace.setPixelIdxs(rasterPixels,integrationPixels);
    [roiTrace, roiWeight] = hTrace.process(dt,expecteddt);
    roiF(:,idx) = roiTrace(lineIxs);
    roiW(:,idx) = roiWeight(lineIxs);
end
E.ROIs.F = roiF;
E.ROIs.weight = roiW;
clear roiF roiWeight
end

% DSframes = ceil(firstLine:dt:lastLine);
% nDSframes= length(DSframes); %number of downsampled frames
% lineIndices = [];
% for lix = 1:numLinesPerCycle
%     if any(intersect(hLowLevelDataFile.lineSuperPixelIDs{lix},floor(lut.allSuperPixelIDs/100)))
%         lineIndices = [lineIndices; lix];
%     end
% end
% cycleIndices = (floor((firstLine-1) / numLinesPerCycle)+1) : (floor((lastLine-1) / numLinesPerCycle)+1);
% 
% [LI, CI] = ndgrid(lineIndices, cycleIndices);
% allLineData = hLowLevelDataFile.getLineData(LI(:), CI(:));
% 
% for DSframeIx = 1:nDSframes
%    timeWindow = max(1,floor(DSframes(DSframeIx)-2*dt)):min(ceil(DSframes(DSframeIx)+2*dt),hLowLevelDataFile.totalNumLines);
% 
%     if ~mod(DSframeIx, 1000)
%         disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
%     end
% lineIndices  = mod(timeWindow-1,numLinesPerCycle)+1;
% cycleIndices = floor((timeWindow-1) / numLinesPerCycle)+1;
% 
% % load data
% data = zeros(length(lut.allSuperPixelIDs),numChannels);
% spCt = zeros(length(lut.allSuperPixelIDs),1);
% allLineData = hLowLevelDataFile.getLineData(lineIndices, cycleIndices);
% 
% for t = 1:length(allLineData)
%     superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineIndices(t)};
% 
%     if numel(superPixIdxs) == 0; continue; end
% 
%     lineData = max(0,allLineData{t});
%     zIdx = hLowLevelDataFile.lineFastZIdxs(lineIndices(t))-1;
% 
%     spID = superPixIdxs*100 + uint32(zIdx); % make superpixel index with Z plane
%     [~,spIdx] = ismember(spID,lut.allSuperPixelIDs);
% 
%     data(spIdx(spIdx>0),:) = data(spIdx(spIdx>0),:) + single(lineData(spIdx>0,:));
%     spCt(spIdx(spIdx>0)) = spCt(spIdx(spIdx>0)) + 1;
% end
% data = data ./ 100; %./ spCt;
% data(spCt == 0,:) = nan;
