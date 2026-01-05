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
    CD = loadTrial(dr, fns{nLoad},fls(nLoad),els(nLoad),selPix,discardFrames{nLoad}, alignData{nLoad}, meanAligned(:,:,:,nLoad), motOutput(:,nLoad), roiData, params);
    E{nLoad}.ROIs = CD.ROIs; E{nLoad}.global = CD.global;
    if i>1 %we process the previous trial after loading the next, to keep CPU usage up during loading
        [E{validTrials(i-1)}, B] = processResult(resultsFuture, E{validTrials(i-1)},params);
        doPlot = false;
        if doPlot
            plotE(E{validTrials(i-1)},Y,B,selPix); %the Y here is from the previous trial
        end
    end
    disp(['Processing dataset: ' fns{nLoad}])
    Y = squeeze(CD.Yobs(:,1,:));
    resultsFuture = extractTrial(Y,CD.Finv, sources, any(selPix,3), params);
    clear CD;
end
[E{nLoad}, B] = processResult(resultsFuture,E{nLoad},params);
if doPlot    
    plotE(E{nLoad},Y,B,selPix);
end

end
% ---------------- Helper Functions ----------------
function plotE(E, Y,B,selPix)
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
        [H,B,S,LS,F0,SNR] = fetchOutputs(resultsFuture);
        E.footprints = H;
        %E.baseline = single(B); %This is a huge variable and currently unused.
        E.dF.events = S;
        E.dF.denoised = convn(S,params.k,'same');
        E.dF.ls = LS;
        E.F0 = F0;
        E.SNR = SNR;
end

function CD = loadTrial(dr, fn, startLine, endLine, selPix, discardFrames, alignData, meanIM, motOutput, roiData, params)
disp('Loading high-res data for file:')
disp([dr filesep fn])

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
        numChannels = params.numChannels;

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
        
        orderedChannels = [params.activityChannel:numChannels, 1:params.activityChannel-1];

        selPx2D = selPx2D(1:size(alignData.viewC,1), 1:size(alignData.viewC,2));
        nPx = numel(selPx2D);
        meanIM = meanIM(1:size(alignData.viewC,1), 1:size(alignData.viewC,2),:);
        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        
        meanPx = reshape(meanIM, numel(labeled), numChannels);
        mLabeled = meanPx(labeled,:);

        blocksize = 600; %number of frames to load at a time; 100 frames ~= 1GB RAM usage
        nBlocks = ceil(nFrames./blocksize);
        blockEdges = round(linspace(1, nFrames+1, nBlocks+1));

        sumF = sum(meanPx,1,'omitmissing');
        IMsel = nan(sum(selPx2D(:)),numChannels, nFrames);
        Finvsel = nan(sum(selPx2D(:)),nFrames);
        for bix = nBlocks:-1:1
            fIxs  = blockEdges(bix):(blockEdges(bix+1)-1);
            [Y, Fresh] = S2data.getImages(orderedChannels,frameLines(fIxs),ceil(dt),1,1);
            Y = Y(alignData.trimRows, alignData.trimCols,:,:);
            Fresh = Fresh(alignData.trimRows, alignData.trimCols,:);
            nFramesInBlock = size(Y,4);

            Y2 = nan(length(viewR),length(viewC),numChannels, nFramesInBlock);
            Finv = nan(length(viewR),length(viewC),nFramesInBlock);
            parfor frIx = 1:nFramesInBlock
                [Y2(:,:,:,frIx), Finv(:,:,frIx)] = interpFrames(Y(:,:,:,frIx),viewC+motionC(fIxs(frIx)), viewR+motionR(fIxs(frIx)), Fresh(:,:,frIx));
            end

            Y2 = reshape(Y2, nPx, numChannels, nFramesInBlock);
            Finv= reshape(Finv, nPx, nFramesInBlock);

            IMsel(:, :, fIxs) = Y2(selPx2D,:,:);
            Finvsel(:,fIxs) = Finv(selPx2D,:);

            %compute global ROI activity
            yLabeled = double(Y2(labeled(:),:,:)); nans= isnan(yLabeled(:,1,:));
            M = repmat(mLabeled,1,1,nFramesInBlock); M(nans) = nan;
            CD.global.F(orderedChannels,fIxs) = (sum(yLabeled,1, 'omitmissing')./sum(M,1, 'omitmissing')).*sumF;

            %compute user ROI activity
            for rix = length(roiData):-1:1
                mask = roiData{rix}.mask;
                tmp1 = Y2(mask(:),:,:);  %the data over the ROI pixels
                Fpx{rix}(:,:,fIxs) = tmp1;
                tmp2 = repmat(meanPx(mask(:),:), 1,1,numel(fIxs)); %the mean image over the ROI pixels
                nans= isnan(tmp1) | isnan(tmp2);
                tmp1(nans) = 0;
                tmp2(nans) = 0;
                CD.ROIs.F(rix,:,fIxs) = (sum(tmp1,1)./sum(tmp2,1)).*sum(meanPx(mask(:),:),1,'omitmissing'); %dFF over the valid pixels, times the mean
            end
        end

        %perform SVD on user ROIs to denoise
        CD.ROIs.Fsvd = nan(length(roiData), numChannels, nFrames);
        for rix = 1:length(roiData)
            if ~isempty(Fpx{rix}) && ~all(isnan(Fpx{rix}(:)))
                for cix = 1:numel(orderedChannels)
                    Dtmp = squeeze(double(Fpx{rix}(:,cix,:)));
                    [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                    roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
                    [~,selPC] = max(roiLikeness);
                    CD.ROIs.Fsvd(rix,cix,:) = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                end
            end
        end

    case 'bergamo'
        error('bergamo not supported')
        activityChannel = params.activityChannel;
        numChannels = params.numChannels;
        IM = networkScanImageTiffReader([dr filesep fn]);
        IM = double(IM);

        selPx2D = any(selPix,3);

        %rearrange IM into correct dimensions
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []);

        if activityChannel>1
            IM = IM(:,:,[activityChannel:end, 1:activityChannel-1],:);
            disp('Reordering channels for analysis!')
        end
        IM1 = squeeze(IM(:,:,1,:));
        IMsel = interpArray(IM1, selPx2D, motOutput); %interpolate the movie at the shifted coordinates

        if numChannels==2
            IM2 =  squeeze(IM(:,:,2,:));
            clear IM;
            IM2sel = interpArray(IM2, selPx2D, motOutput); %interpolate the movie at the shifted coordinates
            clear IM2;
        else %1 channel
            clear IM;
        end

end

discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,2)))>0; %upsample the discard frames
IMsel(:,discard) = nan;     %throw away movement frames as above
CD.global.F(:,discard) = nan;
CD.ROIs.F(:, :, discard) = nan;
CD.ROIs.Fsvd(:,:,discard) = nan;

% Remove nans from IMsel and Finvsel
nans = isnan(IMsel);
nanFill = IMsel;
nanFill(all(nans,3)) = 0;
while any(isnan(nanFill), 'all')
    nanFill = smoothdata(nanFill,3,"movmean",params.baselineWindow_samps, 'omitmissing');
end
IMsel(nans) = nanFill(nans); 
Finvsel(squeeze(nans(:,1,:))) = 1000*mean(Finvsel,'all', 'omitmissing');
CD.Yobs = IMsel; 
CD.Finv = Finvsel;
CD.discardFrames = discard;
end