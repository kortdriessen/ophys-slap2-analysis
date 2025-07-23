function  [IM, meanIM, IMc, aData, peaks, discardFrames]= loadAndProcessTrialAsync(dr, fn, numChannels, params)
     %load the tiff
    IM = copyReadDeleteScanImageTiff([dr filesep fn]);
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;

    meanIM= mean(IM,4, 'omitnan');
    nanPx = repmat(mean(isnan(IM), [3 4])>params.nanThresh, 1,1,size(IM,3));
    meanIM(nanPx) = nan;

    IM = squeeze(IM(:,:,params.activityChannel,:));

    %load alignment data
    fnStemEnd = strfind(fn, '_REGISTERED') -1;
    load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
    aData.dsFac = 1; %SLAP2 data does not have downsampling per se
    params.dsFac = aData.dsFac;
    params.alignHz = aData.alignHz;
    nInitFrames = ceil(params.discardInitial_s * params.alignHz); %initial frames to discard

    %discard motion frames
    tmp = aData.recNegErr(:)- medfilt1(aData.recNegErr(:), round(4*params.alignHz)); %smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp = -tmp./(min(-0.005, prctile(tmp,5))); %normalize to the median-to-5th prctile interval, or 0.5% of intensity, whichever is larger
    thresh = params.motionThresh; %decrease thresh to be more stringent on motion correction
    window = 2*ceil(0.1*params.alignHz)+1;% a window in time to censor aronud motion events, ~0.2 seconds;
    discardFrames = imclose(imdilate(tmp>thresh, ones(window,1)) & (tmp>(thresh/2)), ones(window,1));
    discardFrames(1:nInitFrames) = true;
    IM(:,:,discardFrames) = nan;
    if isfield(aData, 'varFacDS')
        vIM = aData.varFacDS;
        aData = rmfield(aData,'varFacDS');
    else
        vIM = [];
    end

    try
        %[IMc, peaks] = localizeSourcesSLAP2(IM, [], params);
       [IMc, peaks] = localizeSources_vIM(IM, vIM, params);
    catch ME
        IM([1 end],:,:) = nan;
        IM(:,[1 end],:) = nan;
        %[IMc, peaks] = localizeSourcesSLAP2(IM, vIM, params);
        [IMc, peaks] = localizeSources_vIM(IM, vIM, params);
    end
end