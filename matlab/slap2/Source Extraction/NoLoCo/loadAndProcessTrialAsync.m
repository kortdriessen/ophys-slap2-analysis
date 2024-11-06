function  [rawIMs, meanIM, IMc, aData, peaks, discardFrames]= loadAndProcessTrialAsync(dr, fn, numChannels, params)
    
    nInitFrames = ceil(params.discardInitial_s/params.frametime); %initial frames to discard

    %load the tiff
    IM = copyReadDeleteScanImageTiff([dr filesep fn]);
    IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;

    meanIM= mean(IM,4, 'omitnan');
    nanPx = repmat(mean(isnan(IM), [3 4])>0.15, 1,1,size(IM,3));
    meanIM(nanPx) = nan;

    rawIMs = squeeze(IM(:,:,1,:));
    clear IM

    %load alignment data
    fnStemEnd = strfind(fn, '_REGISTERED') -1;
    load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
    aData.dsFac = 1; %SLAP2 data does not have downsampling per se
    params.dsFac = aData.dsFac;
    params.frametime = aData.frametime;

    %discard motion frames
    tmp = aData.recNegErr(:)- median(aData.recNegErr(:)); %smoothExp(aData.recNegErr(:),'movmedian', ceil(2/(aData.frametime*aData.dsFac))); %-smoothdata(aData.aRankCorrDS,2, 'movmedian', ceil(2/aData.frametime));
    tmp = tmp./(median(aData.recNegErr(:)) - prctile(aData.recNegErr(:),1)); 
    thresh = params.motionThresh; %decrease thresh to be more stringent on motion correction
    window = 2*ceil(0.1/params.frametime)+1;% a window in time to censor aronud motion events, ~0.2 seconds;
    discardFrames = imclose(imdilate(tmp>thresh, ones(window,1)) | (tmp>(thresh/2) & imdilate(tmp>thresh, ones(2*window+1,1))), ones(window,1));
    discardFrames(1:nInitFrames) = true;
    rawIMs(:,:,discardFrames) = nan;
    
    
    [IMc, peaks] = localizeSourcesSLAP2(rawIMs, aData, params);
end