function summary = extractDendrites(dr_or_pathToTrialTable, paramsIn)
import ScanImageTiffReader.ScanImageTiffReader
zIdx = 1;

%PARAMETER SETTING
if nargin>1
    if ischar(paramsIn)  % Parse JSON String to Structure
        paramsIn = jsondecode(paramsIn);
    end
    params = setParams('extractDendrites', paramsIn);
else
    params = setParams('extractDendrites');
end
if ~nargin
    [trialTablefn, dr] = uigetfile('*trialTable*.mat');
else
    %parse dr
    %_or_pathToTrialTable
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trialTable.mat';
    else
        [dr trialTablefn ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

%load trial table
load([dr filesep trialTablefn], 'trialTable');

%use the trial table to make sure dat files have not been merged
numworkers = 16; %we won't make this a parameter because large is always good here, I think
pp = gcp('nocreate');
if isempty(pp) || pp.NumWorkers~=numworkers
    delete(pp);
    parpool(numworkers);
end

nDMDs =size(trialTable.filename,1);
rTot = 0;
for DMDix = 1:nDMDs
    %load image data
    firstValidTrial = find(~isempty(trialTable.filename(DMDix,:)), 1, 'first');
    hDataFile = slap2.util.DataFile([dr filesep trialTable.filename{DMDix,firstValidTrial}]);

    imagingRois = hDataFile.metaData.AcquisitionContainer.ROIs.rois;
    isIntegration = false(1,length(imagingRois));
    for rix = 1:length(imagingRois)
        isIntegration(rix) = strcmpi(imagingRois{rix}.imagingMode, 'Integrate');
    end
    integrationRois = imagingRois(isIntegration);
    nImagingROIs = numel(integrationRois);

    %GET ROI outlines
    outlines = {};
    maskImage{DMDix} = -1.*ones(hDataFile.metaData.dmdPixelsPerColumn, hDataFile.metaData.dmdPixelsPerRow);
    for rix = 1:nImagingROIs
        shape = integrationRois{rix}.shapeData;
        tmp = false(hDataFile.metaData.dmdPixelsPerColumn, hDataFile.metaData.dmdPixelsPerRow);
        tmp(sub2ind(size(tmp), shape(:,1), shape(:,2))) = true;
        summary.masks{DMDix}(:,:,rix) = tmp; 
        if DMDix==1
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rix;
        else
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rix+size(summary.masks{1},3);
        end
        %make an image where each pixel corresponds to the mask ID
        outlines = cat(1, outlines, bwboundaries(tmp));
    end

    %load reference image
    refFn = dir([dr filesep '**' filesep '*DMD' int2str(DMDix) '*REFERENCE*.tif']);
    assert(length(refFn)==1);
    A = ScanImageTiffReader([refFn.folder filesep refFn.name]);
    IDs = A.descriptions;
    for im_ix = length(IDs):-1:1
        js = jsondecode(IDs{im_ix});
        z(im_ix) = js.z;
        ch(im_ix) = js.channel;
    end
    nChan = numel(unique(ch));
    Zs = unique(z);
    stack = A.data();
    stack = reshape(stack, size(stack,1), size(stack,2), nChan, []);

    %get imaging Z plane
    metaZ = hDataFile.metaData.AcquisitionContainer.ParsePlan.zs;
    [~, bestZix] = min(abs(Zs-metaZ)); bestZix = bestZix(1);
    %bestZ = Zs(bestZix);

    %The image plane to trace on
    imPlane = permute(stack(:,:,:,bestZix),[2 1 3]);
    summary.refIM{DMDix} = imPlane;

    if params.manualROIs
        hROIs = drawROIs(imPlane, dr, refFn.name, outlines);
        waitfor(hROIs.hF);
        nAnalysisROIs(DMDix) = numel(hROIs.roiData);
        summary.masks{DMDix} = false(0);
        maskImage{DMDix}(:,:) = -1;
        for rix = 1:nAnalysisROIs(DMDix)
            rTot = rTot+1;
            summary.masks{DMDix}(:,:,rix) = hROIs.roiData{rix}.mask;
            %make an image where each pixel corresponds to the mask ID
            maskImage{DMDix}(summary.masks{DMDix}(:,:,rix)) = rTot; 
        end
    else
        nAnalysisROIs(DMDix) = nImagingROIs; 
    end

    figure('name', ['Reference Image for DMD' int2str(DMDix)]),
    imshow(summary.refIM{DMDix}(:,:,1),[]);
    figure('name', ['Masks for DMD' int2str(DMDix)]), 
    imshow(maskImage{DMDix},[]); colormap('jet')
    drawnow
end


%if continuous acquisition
if contains(trialTable.filename{1}, '-CYCLE-')
    ydata = nan(0, sum(nAnalysisROIs)); 
    r0 = 0;
    for DMDix = 1:nDMDs
            clear futures
            hDataFile = slap2.util.MultiDataFiles([dr filesep trialTable.filename{DMDix,1}]);
            for rix =  1:nAnalysisROIs(DMDix)
                hTrace = slap2.util.datafile.trace.Trace(hDataFile,zIdx,params.chIdx);
                pixelMask = summary.masks{DMDix}(:,:,rix);  % false(800,1280);
                rasterPixels = pixelMask;
                integrationPixels = pixelMask;
                hTrace.setPixelIdxs(rasterPixels,integrationPixels);
                futures(rix) = hTrace.processAsync(params.windowWidth_lines,params.expectedWindowWidth_lines);
            end
            for rix = 1:nAnalysisROIs(DMDix)
                tmp = futures(rix).fetchOutputs();
                ydata(1:length(tmp),r0+rix) = tmp;
            end
            r0 = r0 + nAnalysisROIs(DMDix);
    end
    summary.traces = ydata;
else
%for each trial
nTrials = size(trialTable.filename, 2);
for tix = 1:nTrials
    disp(['Processing trial: ' int2str(tix)])
    rTot = 0;
    clear futures
    for DMDix = 1:nDMDs
        hDataFile = slap2.util.DataFile([dr filesep trialTable.filename{DMDix,tix}]);
        for rix =  1:nAnalysisROIs(DMDix)
            rTot = rTot+1;
            hTrace = slap2.util.datafile.trace.Trace(hDataFile,zIdx,params.chIdx);
            pixelMask = summary.masks{DMDix}(:,:,rix);  % false(800,1280);
            rasterPixels = pixelMask;
            integrationPixels = pixelMask;
            hTrace.setPixelIdxs(rasterPixels,integrationPixels);
            futures(rTot) = hTrace.processAsync(params.windowWidth_lines,params.expectedWindowWidth_lines);
        end
    end
    waitfor(futures);
    rS = 0;
    ydata = nan(round(trialTable.lastLine(1, tix)-trialTable.firstLine(1, tix)+1), rTot);
    for DMDix = 1:nDMDs
        for rix = 1:nAnalysisROIs(DMDix)
            rS = rS+1;
            tmp = futures(rS).fetchOutputs();
            tmp = tmp(uint32(trialTable.firstLine(DMDix, tix)):min(length(tmp),uint32(trialTable.lastLine(DMDix, tix))));
            ydata(1:length(tmp),rS) = tmp;
        end
    end
    summary.traces{tix} = ydata;
    
    %plots
    sNorm = summary.traces{tix};
    sNorm = sNorm-mean(sNorm,1, 'omitnan');
    sNorm = sNorm./std(sNorm,0,1, 'omitnan');
    figure(99), imagesc(sNorm')
    title(['Trial' int2str(tix)])   
end
end
summary.nAnalysisROIs = nAnalysisROIs;
save([dr filesep 'dendriticVoltageSummary-' datestr(now, 'YYmmDD-HHMMSS') '.mat'], 'summary', '-v7.3');


