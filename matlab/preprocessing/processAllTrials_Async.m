function E = processAllTrials_Async(dr, fns, fls, els, selPix, sources, discardFrames, alignData, meanAligned, motOutput, roiData, validTrials, params)
curPool = gcp('nocreate');
if isempty(curPool) || ~strcmpi(class(curPool), 'parallel.ThreadPool') %  ~strcmpi(class(curPool), 'parallel.ProcessPool') %
    delete(curPool);
    parpool('Threads'); %use all available Threads
end
numDatasets = numel(fls);
E = cell(numDatasets,1); 
for i = 1:numel(validTrials)
    nLoad = validTrials(i);
    if isfield(params, 'batchIx') && isfield(params, 'nBatches')
        fprintf('[Batch %d/%d, Trial %d/%d]\n', params.batchIx, params.nBatches, i, numel(validTrials));
    else
        fprintf('[Trial %d/%d]\n', i, numel(validTrials));
    end
    CD = loadTrial(dr, fns{nLoad},fls(nLoad),els(nLoad),selPix,discardFrames{nLoad}, alignData{nLoad}, meanAligned(:,:,:,nLoad), motOutput(:,nLoad), roiData, params);
    E{nLoad}.ROIs = CD.ROIs; E{nLoad}.global = CD.global;
    E{nLoad}.discardFrames = CD.discardFrames;
    if i>1 %we process the previous trial after loading the next, to keep CPU usage up during loading
        [E{validTrials(i-1)}, B] = processResult(resultsFuture, E{validTrials(i-1)},params);
    end
    disp(['Processing dataset: ' fns{nLoad}])
    %Y = squeeze(CD.Yobs(:,1,:));
    Y = permute(CD.Yobs, [1 3 2]);
    resultsFuture = extractTrial(Y,CD.Finv, sources, any(selPix,3), params);
    clear CD;
end

[E{nLoad}, B] = processResult(resultsFuture,E{nLoad},params);
validateTrialResults(E, validTrials, params.numChannels);

end
% ---------------- Helper Functions ----------------
function plotE(E, Y,B,selPix) %debugging helper only
    %generate activity movie, baseline movie, and residual
    sel2D = any(selPix,3);
    Ht = reshape(E.footprints,numel(sel2D),[]);
    Ht = Ht(sel2D(:),:);

    Amov = max(0,Ht,'omitmissing')*max(0,E.dF.denoised,'omitmissing');
    Rmov = Y - Amov - B;

    render = zeros([size(sel2D) 500]);
    render(repmat(sel2D,1,1,500)) = B(:,501:1000);
end

function [E, B] = processResult(resultsFuture, E, params)
[H,B,Sanalysis,LSanalysis,F0analysis,SNR] = fetchOutputs(resultsFuture);
contract = channelContract(params.numChannels, params.activityChannel);
assert(size(Sanalysis,3) == params.numChannels && ...
    size(LSanalysis,3) == params.numChannels && ...
    size(F0analysis,3) == params.numChannels, ...
    'slap2:extractionChannelCountMismatch', ...
    'Extraction returned a channel count that differs from params.numChannels.');

discard = E.discardFrames;
E.footprints = H;
% E.baseline is intentionally not saved because it is very large.
E.dF.events = analysisToRawChannels(Sanalysis, ...
    contract.analysisChannelOrder, 3);
E.dF.events(:,discard,:) = nan;

denoisedAnalysis = nan(size(Sanalysis), 'like', Sanalysis);
denoisedAnalysis(:,:,1) = convn(Sanalysis(:,:,1), params.k, 'same');
for channelIx = 2:size(Sanalysis,3)
    denoisedAnalysis(:,:,channelIx) = convn(Sanalysis(:,:,channelIx), params.k2, 'same');
end
E.dF.denoised = analysisToRawChannels(denoisedAnalysis, ...
    contract.analysisChannelOrder, 3);
E.dF.denoised(:,discard,:) = nan;

E.dF.ls = analysisToRawChannels(LSanalysis, ...
    contract.analysisChannelOrder, 3);
E.dF.ls(:,discard,:) = nan;

E.F0 = analysisToRawChannels(F0analysis, ...
    contract.analysisChannelOrder, 3);
E.SNR = SNR;
end

function CD = loadTrial(dr, fn, startLine, endLine, selPix, discardFrames, alignData, meanIM, motOutput, roiData, params)
disp('Loading high-res data for file:')
disp([dr filesep fn])

numChannels = params.numChannels;
contract = channelContract(numChannels, params.activityChannel);
rawChannelOrder = contract.rawChannelOrder;
analysisChannelOrder = contract.analysisChannelOrder;

switch params.microscope
    case 'SLAP2'
        if params.includeIntegrationROIs
            warning('includeIntegration not implemented, using raster only!')
            spTypeFlag = 1; %use only raster superpixels
        else
            spTypeFlag = 1; %use only raster superpixels
        end

        %load the high time resolution data
        S2data = slap2.Slap2DataFile([dr filesep fn]);
        meta = loadMetadata([dr filesep fn]);
        
        linerateHz = 1/meta.linePeriod_s;
        dt = linerateHz/params.analyzeHz;
        frameLines = ceil(startLine:dt:endLine);
        nFrames= length(frameLines);
        selPx2D = any(selPix,3);

        %upsample motion
        motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
        motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);
        viewC = alignData.viewC(1,:);
        viewR = alignData.viewR(:,1);

        selPx2D = selPx2D(1:size(alignData.viewC,1), 1:size(alignData.viewC,2));
        nPx = numel(selPx2D);
        meanIM = meanIM(1:size(alignData.viewC,1), 1:size(alignData.viewC,2),:);
        assert(size(meanIM,3) == numChannels, 'slap2:staticChannelCountMismatch', ...
            'Registered mean image has %d channels; expected %d.', size(meanIM,3), numChannels);
        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        
        meanIMAnalysis = rawToAnalysisChannels(meanIM, analysisChannelOrder, 3);
        meanPxAnalysis = reshape(meanIMAnalysis, numel(labeled), numChannels);
        mLabeledAnalysis = meanPxAnalysis(labeled,:);

        blocksize = 600; %number of frames to load at a time; 100 frames ~= 1GB RAM usage
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumFAnalysis = sum(meanPxAnalysis,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            [Yraw, Fresh] = S2data.getImages(rawChannelOrder,frameLines(fIxs),ceil(dt),1,1);
            assert(size(Yraw,3) == numChannels, 'slap2:dynamicChannelCountMismatch', ...
                'Raw reader returned %d channels; expected %d.', size(Yraw,3), numChannels);
            Yanalysis = rawToAnalysisChannels(Yraw, analysisChannelOrder, 3);
            Yanalysis = Yanalysis(alignData.trimRows, alignData.trimCols,:,:);
            Fresh = Fresh(alignData.trimRows, alignData.trimCols,:);
            nFramesInBlock = size(Yanalysis,4);

            Y2 = nan(length(viewR),length(viewC),numChannels, nFramesInBlock);
            Finv = nan(length(viewR),length(viewC),nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Y2(:,:,:,frIx), Finv(:,:,frIx)] = interpFrames(Yanalysis(:,:,:,frIx),viewC+motionC(fIxs(frIx)), viewR+motionR(fIxs(frIx)), Fresh(:,:,frIx));
            end

            Y2 = reshape(Y2, nPx, numChannels, nFramesInBlock);
            Finv= reshape(Finv, nPx, nFramesInBlock);

            IMsel(:, :, fIxs) = Y2(selPx2D,:,:);
            Finvsel(:,fIxs) = Finv(selPx2D,:);

            %compute global ROI activity
            yLabeled = double(Y2(labeled(:),:,:)); nans= isnan(yLabeled);
            M = repmat(mLabeledAnalysis,1,1,nFramesInBlock); M(nans) = nan;
            CD.global.F(analysisChannelOrder,fIxs) = ...
                (sum(yLabeled,1, 'omitmissing')./sum(M,1, 'omitmissing')).*sumFAnalysis;

            %compute user ROI activity
            for rix = length(roiData):-1:1
                mask = roiData{rix}.mask;
                tmp1 = Y2(mask(:),:,:);  %the data over the ROI pixels
                Fpx{rix}(:,:,fIxs) = tmp1;
                tmp2 = repmat(meanPxAnalysis(mask(:),:), 1,1,numel(fIxs)); %the matching mean image in analysis order
                nans= isnan(tmp1) | isnan(tmp2);
                tmp1(nans) = 0;
                tmp2(nans) = 0;
                CD.ROIs.F(rix,analysisChannelOrder,fIxs) = ...
                    (sum(tmp1,1)./sum(tmp2,1)).*sum(meanPxAnalysis(mask(:),:),1,'omitmissing'); %saved in raw order
            end
        end

        %perform SVD on user ROIs to denoise
        CD.ROIs.Fsvd = nan(length(roiData), numChannels, nFrames);
        for rix = 1:length(roiData)
            if ~isempty(Fpx{rix}) && ~all(isnan(Fpx{rix}(:)))
                for cix = 1:numel(analysisChannelOrder)
                    Dtmp = squeeze(double(Fpx{rix}(:,cix,:)));
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,analysisChannelOrder(cix),:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end

    case 'bergamo'
        baseline = prctile(reshape(meanIM, [], numChannels), 10,1);

        if endsWith(fn, '.h5')
            desc = h5info([dr filesep fn]);
            IM = h5read([dr filesep fn], ['/', desc.Datasets.Name]);
        else
            IM = networkScanImageTiffReader([dr filesep fn]); %Raw movie
        end
        IM = double(IM);

        %rearrange into correct dimensions
        selPx2D = any(selPix,3);
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []);
        IM = IM-reshape(baseline, [1 1 numChannels 1]); %subtract baseline
        nPx = size(IM,1)*size(IM,2);
        nFrames = size(IM,4);

        %upsample motion
        viewC = (1:size(IM,2)) + motOutput(2);
        viewR = (1:size(IM,1))' + motOutput(1);
        
        meanIM = meanIM(1:size(IM,1), 1:size(IM,2),:);
        assert(size(meanIM,3) == numChannels, 'slap2:staticChannelCountMismatch', ...
            'Registered mean image has %d channels; expected %d.', size(meanIM,3), numChannels);
        scored = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        scored = (scored-prctile(scored(:),25,'all'))./(prctile(scored(:),99,'all')-prctile(scored(:),25,'all'));
        labeled = scored>0.1;
        meanIMAnalysis = rawToAnalysisChannels(meanIM, analysisChannelOrder, 3);
        meanPxAnalysis = reshape(meanIMAnalysis, numel(labeled), numChannels);
        mLabeledAnalysis = meanPxAnalysis(labeled,:);

        blocksize = 600; %number of frames to load at a time;
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumFAnalysis = sum(meanPxAnalysis,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        Fresh = ones(size(IM,1), size(IM,2), 'single'); %freshness currently unused for bergamo
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            nFramesInBlock = length(fIxs);
            
            Yraw = IM(:,:,:,fIxs); %canonical physical raw order
            assert(size(Yraw,3) == numChannels, 'slap2:dynamicChannelCountMismatch', ...
                'Movie has %d channels; expected %d.', size(Yraw,3), numChannels);
            Yanalysis = rawToAnalysisChannels(Yraw, analysisChannelOrder, 3);
            Y2 = nan(length(viewR),length(viewC),numChannels, nFramesInBlock);
            Finv = nan(length(viewR),length(viewC),nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Y2(:,:,:,frIx), Finv(:,:,frIx)] = interpFrames(Yanalysis(:,:,:,frIx),viewC, viewR, Fresh);
            end

            Y2 = reshape(Y2, nPx, numChannels, nFramesInBlock);
            Finv= reshape(Finv, nPx, nFramesInBlock);

            IMsel(:, :, fIxs) = Y2(selPx2D,:,:);
            Finvsel(:,fIxs) = Finv(selPx2D,:);

            %compute global ROI activity
            yLabeled = double(Y2(labeled(:),:,:)); nans= isnan(yLabeled);
            M = repmat(mLabeledAnalysis,1,1,nFramesInBlock); M(nans) = nan;
            CD.global.F(analysisChannelOrder,fIxs) = ...
                (sum(yLabeled,1, 'omitmissing')./sum(M,1, 'omitmissing')).*sumFAnalysis;

            %compute user ROI activity
            for rix = length(roiData):-1:1
                mask = roiData{rix}.mask;
                tmp1 = Y2(mask(:),:,:);  %the data over the ROI pixels
                Fpx{rix}(:,:,fIxs) = tmp1;
                tmp2 = repmat(meanPxAnalysis(mask(:),:), 1,1,numel(fIxs)); %the matching mean image in analysis order
                nans= isnan(tmp1) | isnan(tmp2);
                tmp1(nans) = 0;
                tmp2(nans) = 0;
                CD.ROIs.F(rix,analysisChannelOrder,fIxs) = ...
                    (sum(tmp1,1)./sum(tmp2,1)).*sum(meanPxAnalysis(mask(:),:),1,'omitmissing'); %saved in raw order
            end
        end

        %perform SVD on user ROIs to denoise
        CD.ROIs.Fsvd = nan(length(roiData), numChannels, nFrames);
        for rix = 1:length(roiData)
            if ~isempty(Fpx{rix}) && ~all(isnan(Fpx{rix}(:)))
                for cix = 1:numel(analysisChannelOrder)
                    Dtmp = squeeze(double(Fpx{rix}(:,cix,:)));
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,analysisChannelOrder(cix),:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end
end

discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,3)))>0; %upsample the discard frames
IMsel(:,:,discard) = nan;     %throw away movement frames as above
IMsel(IMsel<-2*std(IMsel,0,3, 'omitmissing')) = nan; %throw away any aberrantly negative data
CD.global.F(:,discard) = nan;
CD.ROIs.F(:, :, discard) = nan;
CD.ROIs.Fsvd(:,:,discard) = nan;

% Remove nans from IMsel and Finvsel
nans = isnan(IMsel);
nanFill = IMsel;
nanFill(all(nans,3)) = 0;
for channelIx = 1:numChannels
    while any(isnan(nanFill(:,channelIx,:)), 'all')
        nanFill(:,channelIx,:) = smoothdata(nanFill(:,channelIx,:),3,"movmean", ...
            params.analysisBaselineWindows_samps(channelIx), 'omitmissing');
    end
end
IMsel(nans) = nanFill(nans); 
Finvsel(squeeze(nans(:,1,:))) = 1000*mean(Finvsel,'all', 'omitmissing');
CD.Yobs = IMsel; 
CD.Finv = Finvsel;
CD.discardFrames = discard;
end
