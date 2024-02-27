function exptSummary = processTrialAsync(dr, fnRaw, startLine, endLine, W0, F0selDS, selPix, discardFrames, alignData, meanIM, motOutput, roiData, params)
    fn = fnRaw;
    disp('Loading high-res data for file:')
    disp([dr filesep fn])

    %load the high time resolution data
    S2data = slap2.Slap2DataFile([dr filesep fn]);
    meta = loadMetadata([dr filesep fn]);
    numChannels = S2data.numChannels;
    %numLines = S2data.totalNumLines;
    linerateHz = 1/meta.linePeriod_s;
    dt = linerateHz/params.analyzeHz;
    frameLines = ceil(startLine:dt:endLine);
    nFrames= length(frameLines);
    selPx2D = any(selPix,3);
    sz = size(selPx2D);

    %upsample motion
    motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
    motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);
    
    labeled = medfilt2(meanIM(:,:,1), [3 3]);
    labeled = labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
    for cix = numChannels:-1:1
        tmp = double(meanIM(:,:,cix));
        mLabeled{cix} = tmp(labeled);
    end

    Ysz = [length(alignData.trimRows) length(alignData.trimCols)];
    %interpolate the raw data at the shifted coordinates
    for fix = nFrames:-1:1
        for cix = 1:numChannels
            Y = S2data.getImage(cix, frameLines(fix), ceil(dt), 1);
            Y = Y(alignData.trimRows, alignData.trimCols);
            Y = interp2(1:Ysz(2), 1:Ysz(1), Y,alignData.viewC+motionC(fix), alignData.viewR+motionR(fix), 'linear', nan);
            if cix==1
                IMsel(:, fix) = Y(selPx2D);
            elseif cix==2
                IM2sel(:, fix) = Y(selPx2D);
            end

            %compute global ROI activity
            mIM = meanIM(:,:,cix);
            yLabeled = double(Y(labeled)); nans= isnan(yLabeled);
            FF = (sum(yLabeled(~nans))./sum(mLabeled{cix}(~nans))).*sum(mLabeled{cix});
            exptSummary.global.F(fix,cix) = FF;

            %compute user ROI activity
            for rix = 1:length(roiData)
                mask = roiData{rix}.mask;
                tmp1 = Y(mask); tmp2 = mIM(mask);
                nans= isnan(tmp1);
                FF = (sum(tmp1(~nans))./sum(tmp2(~nans))).*sum(tmp2);
                exptSummary.ROIs.F(rix, fix,cix) = FF;
            end

            % if fix<200 && cix==1
            %     sanitycheckTMP(:,:,201-fix) = Y;
            % end
        end
    end

    %IMsel = IMsel - min(0, min(mean(IMsel,2, 'omitnan'))); %ensure that the baseline is not overestimated
    
    discard = reshape(repmat(discardFrames(:), 1,params.dsFac)', 1,[]); %upsample the discard frames
    IMsel(:,discard) = nan;     %throw away movement frames as above

    [IMsel, F0sel, W1, selNans] = prepareNMFproblem(IMsel, W0, F0selDS, params);

    %NMF
    H0 = ones(size(W1,2), size(IMsel,2));
    for iter = 1:3 %perform nonnegative matrix division to initialize H0 with W0 constant
        H0 = H0.*(max(0,W1'*IMsel))./((W1'*W1)*(H0 + mean(H0(:)/100))); %confirm this is right; per Lee and Seung NIPS 2000 'Algorithms for nonnegative matrix factorization'
    end
    opts = statset('MaxIter', 10,  'Display', 'final');
    [W,H] = nnmf2(IMsel, size(W1,2),'algorithm', 'mult', 'w0', W1, 'h0', H0, 'options', opts); %!!nnmf2 has been modified to keep the ordering of the provided factors
    Wfull = nan([sz(1)*sz(2) size(W,2)]);
    Wfull(any(selPix,3),:) = W;
    Wfull = reshape(Wfull, sz(1),sz(2), []);
    
    %Frames where more than 25% of pixels in a source are nan
    nanFramesH = false(size(H));
    for sourceIx = 1:size(H,1)
        support = W(:,sourceIx)>0;
        nanFramesH(sourceIx,:) = mean(selNans(support,:)) > 0.25;
    end

    %deconvolve out matched filter we applied earlier
    kernel = [zeros(1,ceil(8*params.tau_full)) exp(-(0:ceil(8*params.tau_full))/params.tau_full)];
    fKernel = fliplr(kernel);
    kernel = kernel./sum(kernel);
    doubleKernel = conv(kernel, fKernel, 'same');
    doubleKernel = doubleKernel./sum(doubleKernel);
    
    %perform deconvolution, filling in NaNs with reconstructed values every
    %few iterations:
    J = deconvlucy({H},doubleKernel, 20);
    Js = deconvlucy({H},fKernel, 20);
    for iter = 1:5
            recon = convn(J{2}, doubleKernel, 'same');
            J{1}(nanFramesH) = recon(nanFramesH);
            J = deconvlucy(J,doubleKernel, 25);

            recon2 =  convn(Js{2}, fKernel, 'same');
            Js{1}(nanFramesH) = recon2(nanFramesH);
            Js = deconvlucy(Js,fKernel, 25);
    end
    H2 = J{2};
    H3 = convn(H2, kernel, 'same');
    H4 = Js{2};

    %errH = (H - convn(H2, doubleKernel, 'same')).^2;
    %errH = sqrt(convn(errH, doubleKernel, 'same')); % uncertainty at each point

    %compute F0
    F0mean = repmat(mean(F0sel,2, 'omitnan'),1,size(F0sel,2));
    F0sel(~isfinite(F0sel)) = max(0, F0mean(~isfinite(F0sel)));
    F0= (W./max(W,[],1))' *F0sel;
    F0(nanFramesH) = nan;
    
    %ensure F0 is positive
    Fmin = min(F0(:));
    desMin = prctile(F0(:),1) - Fmin;
    if Fmin<desMin
        F0 = F0 + desMin - Fmin;
    end

    %NaN out invalid data
    H(nanFramesH) = nan; %The match filtered signal
    H2(nanFramesH) = nan; %The detected events
    H3(nanFramesH) = nan; %The denoised activity

    %exptSummary.dFerr = sum(W,1)'.*errH;
    exptSummary.matchFilt(:,:,1) = sum(W,1)'.*H; %[source#, time, channel]
    exptSummary.events(:,:,1) = sum(W,1)'.*H2; %[source#, time, channel]
    exptSummary.denoised(:,:,1) = sum(W,1)'.*H3; %[source#, time, channel]
    exptSummary.dFraw(:,:,1) = sum(W,1)'.*H4; %[source#, time, channel]
    exptSummary.F0(:,:,1) = F0;
    exptSummary.footprints = Wfull;
    
    %exptSummary.sanityCheck = mean(sanitycheckTMP,3, 'omitnan');

    %compute channel 2 signals
    if numChannels==2
        F_2 = (W./max(W,[],1))' * IM2sel;
        F0_2 = computeF0(F_2', ceil(params.denoiseWindow_s*params.analyzeHz), ceil(params.baselineWindow_Ca_s*params.analyzeHz),1)';
        exptSummary.dFraw(:,:,2) = F_2;
        exptSummary.F0(:,:,2) = F0_2;
    end

    exptSummary.dFF = exptSummary.dFraw./exptSummary.F0;
end