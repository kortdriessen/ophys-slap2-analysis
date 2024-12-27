function results =  makeXValCovPlot(fns, HPfreqs, useBackground)

if nargin<1 || isempty(fns)
    fns = uipickfiles('FilterSpec','*_REGISTERED_RAW.tif');
end

if nargin < 2
    HPfreqs = 0;
end

if nargin < 3
    useBackground = false;
end

if ~iscell(fns)
    fns = {fns};
end

nTrials = 15;
nXVals = 20;
activityChannel = 2;

rng("default")

xval_splits = zeros(floor(nTrials/2),2,nXVals);

for idx = 1:nXVals
    xval_splits(:,:,idx) = reshape(randperm(nTrials,nTrials - mod(nTrials,2)),size(xval_splits,1:2));
end

results = table('Size',[0 8],'VariableTypes',{'single','single','single','single','uint16','string','uint16','uint16'},'VariableNames',["HP_freqs","distance","xval_cov","xval_cor","fileIx","file","pix1","pix2"]);

for fileIx = 1:length(fns)
    fn = fns{fileIx};
    disp(['Processing: ' fn])
    [dr,fn, ext] = fileparts(fn);
    fn = strcat(fn,ext);

    % get acquisition index
    [idxStart, idxEnd] = regexp(fn,'_\d*_');
    acqIdx = str2num(fn(idxStart+1:idxEnd-1));

    % get harp directory
    [idxStart, ~] = regexp(dr,'[\\/]scans[\\/]');
    harpDir = [dr(1:idxStart) 'harp'];

    photodiode = readmatrix(fullfile(harpDir,['photodiode_' num2str(acqIdx-1) '.csv']),'NumHeaderLines',1);
    orientations = readmatrix(fullfile(harpDir,['orientations_' num2str(acqIdx-1) '.csv']));
    
    opts = delimitedTextImportOptions("NumVariables", 2);
    
    % Specify range and delimiter
    opts.DataLines = [1, Inf];
    opts.Delimiter = ",";
    
    % Specify column names and types
    opts.VariableNames = ["Flyback", "Timestamp"];
    opts.VariableTypes = ["categorical", "double"];
    
    % Specify file level properties
    opts.ExtraColumnsRule = "ignore";
    opts.EmptyLineRule = "read";
    
    % Specify variable properties
    opts = setvaropts(opts, "Flyback", "EmptyFieldRule", "auto");
    
    % Import the data
    frameClk = readtable(fullfile(harpDir,['frameClk_' num2str(acqIdx-1) '.csv']), opts);
    
    clear numIdx opts
    
    frameClkSignal = table2array(frameClk(:,1))=='True';
    diffFrameClk = diff(frameClkSignal);
    
    frameClkTrue = frameClk([false; diffFrameClk == 1],2);
    frameClkTrue = table2array(frameClkTrue);
    disp([num2str(length(frameClkTrue)) ' flybacks detected']);
    
    uniqueOrientations = sort(unique(round(orientations * 360 / 2 / pi)));
    
    pd_mixture = fitgmdist(photodiode(:,2),2);
    [~, onCluster] = max(pd_mixture.mu);
    pdOn = cluster(pd_mixture,photodiode(:,2)) == onCluster;
    
    onTimes = photodiode(diff(pdOn) > 0,1);
    offTimes = [photodiode(diff(pdOn) < 0,1); photodiode(end,1)];


    % read alignment data
    [idxStart, ~] = regexp(fn,'_REGISTERED_RAW.tif');
    load(fullfile(fullfile(dr,[fn(1:idxStart) 'ALIGNMENTDATA.mat'])));

    [Ad, ~, ~] = networkScanImageTiffReader(fullfile(dr,fn));
    [nRow, nCol, ~] = size(Ad);

    Ad = reshape(Ad, size(Ad,1), size(Ad,2), aData.numChannels, []);
    disp([num2str(size(Ad,4)) ' frames read']);
    IM = reshape(Ad(:,:,min(end,activityChannel),:),[],size(Ad,4));
    clear Ad;
    
    frameClkTrue = frameClkTrue(1:min(length(frameClkTrue),size(IM,2)));
    clkErrs = find(diff(frameClkTrue)<0)+1;
    frameClkTrue(clkErrs) = (frameClkTrue(clkErrs-1)+frameClkTrue(clkErrs+1))/2;

    validPix = find(mean(isnan(IM),2) < 0.2);

    meanIM = mean(IM(validPix,:),2,'omitnan');
    meanIM = meanIM - min(meanIM);
    meanIM = meanIM / max(meanIM);

    thresh = graythresh(meanIM);

    numLabeled = sum(meanIM>thresh);
    thresh = prctile(meanIM,(100-numLabeled/numel(meanIM)*150));
    
    if useBackground
        signalMask = zeros(nRow,nCol);
        signalMask(validPix(meanIM > thresh)) = 1;
        backgroundMask = ~imdilate(signalMask,ones(9,9));
        backgroundMask(mean(isnan(IM),2) >= 0.2) = false;
        labeledPix = find(backgroundMask(:));
    else
        labeledPix = validPix(meanIM > thresh);
    end

    IMsel = IM(labeledPix,:);
    nanPix = isnan(IMsel);
    IMsel(nanPix) = 0;

    trialTraces = nan([size(IMsel,1), length(orientations) / length(uniqueOrientations), length(uniqueOrientations), ceil(2 / aData.frametime), length(HPfreqs)]);
    
    for freq_ix = 1:length(HPfreqs)
        if HPfreqs(freq_ix) == 0
            IMsel(nanPix) = nan;
            IMselHP = IMsel - movmedian(conv2(IMsel,ones(1,5)/5,'same')',30/aData.frametime,'omitnan')';
        else
            IMselHP = highpass(IMsel',HPfreqs(freq_ix),1/aData.frametime)';
        end
        IMselHP(nanPix) = nan;

        F = griddedInterpolant(frameClkTrue,IMselHP(:,1:length(frameClkTrue))','linear','none');
        
        for i = 1:length(orientations)
            startTime = onTimes(i);
            endTime = offTimes(i);
        
            trialTraces(:,floor((i-1) / length(uniqueOrientations))+1, uniqueOrientations == round(orientations(i) * 360 / 2 / pi),:,freq_ix) ...
                = F(linspace(startTime,startTime+2,size(trialTraces,4)))';
        end
    end

    warning('off','MATLAB:table:RowsAddedExistingVars')

    [rAll, cAll] = ind2sub([nRow, nCol], labeledPix);
    nPix = numel(labeledPix);

    pairs = [];
    for idx = 1:nPix
        pairs = [pairs; [repmat(idx, nPix-idx,1), (idx+1:nPix)']];
    end

    nPairs = size(pairs,1);

    % Pre-allocate result arrays
    pix1Arr = zeros(nPairs,1,'like',labeledPix);
    pix2Arr = zeros(nPairs,1,'like',labeledPix);
    fileIxArr = repmat(fileIx, [nPairs,1]); 
    fileArr = repmat(string(fn), nPairs, 1);
    HPfreqsArr = repmat(reshape(HPfreqs,1,[]),[nPairs,1]);
    distanceArr = zeros(nPairs,1);
    xval_covArr = nan(nPairs,length(HPfreqs));
    xval_corArr = nan(nPairs,length(HPfreqs));

    validPairs = false(nPairs,1);

    parfor p = 1:nPairs
        idx = pairs(p,1);
        jdx = pairs(p,2);
    
        r1 = rAll(idx);
        c1 = cAll(idx);
        r2 = rAll(jdx);
        c2 = cAll(jdx);
    
        distVal = sqrt((r1-r2)^2+(c1-c2)^2)/4;
    
        % Distance threshold check
        if distVal > 15
            continue; % Skip this pair
        end
    
        % Store results
        pix1Arr(p) = labeledPix(idx);
        pix2Arr(p) = labeledPix(jdx);
        distanceArr(p) = distVal;
    
        % Compute cross-validated covariance and correlation
        covs = nan(nXVals,length(HPfreqs));
        corrs = nan(nXVals,length(HPfreqs));
        for xvalIx = 1:nXVals
            t1_idx = xval_splits(:,1,xvalIx);
            t2_idx = xval_splits(:,2,xvalIx);
   
            for freq_ix = 1:length(HPfreqs)
                pix1resp = squeeze(mean(trialTraces(idx,t1_idx,:,:,freq_ix),2,'omitnan'));
                pix2resp = squeeze(mean(trialTraces(jdx,t2_idx,:,:,freq_ix),2,'omitnan'));
        
                % Flatten
                p1 = pix1resp(:);
                p2 = pix2resp(:);
        
                tmp = cov(p1,p2);
                covs(xvalIx,freq_ix) = tmp(1,2);
                corrs(xvalIx,freq_ix) = corr(p1,p2);
            end
        end
    
        xval_covArr(p,:) = mean(covs,1,'omitnan');
        xval_corArr(p,:) = mean(corrs,1,'omitnan');
    
        validPairs(p) = true;
    end

    % Filter out invalid pairs (distance > 15)
    pix1Arr = pix1Arr(validPairs);
    pix2Arr = pix2Arr(validPairs);
    fileIxArr = fileIxArr(validPairs);
    fileArr = fileArr(validPairs);
    HPfreqsArr = HPfreqsArr(validPairs,:);
    distanceArr = distanceArr(validPairs);
    xval_covArr = xval_covArr(validPairs,:);
    xval_corArr = xval_corArr(validPairs,:);

    results_movie = table(pix1Arr, pix2Arr, fileIxArr, fileArr, HPfreqsArr, distanceArr, xval_covArr, xval_corArr,...
    'VariableNames', {'pix1','pix2','fileIx','file','HP_freqs','distance','xval_cov','xval_cor'});

    results = [results; results_movie];


    % for idx = 1:numel(labeledPix)
    %     for jdx = (idx+1):numel(labeledPix)
    % 
    %         [r1,c1] = ind2sub([nRow,nCol],labeledPix(idx));
    %         [r2,c2] = ind2sub([nRow,nCol],labeledPix(jdx));
    % 
    %         if sqrt((r1-r2)^2+(c1-c2)^2)/4 > 15
    %             continue;
    %         end
    % 
    %         if mod(size(results,1),100) == 0
    %             fprintf('Processed %d pairs\n', size(results,1))
    %         end
    % 
    %         results.pix1(end+1) = labeledPix(idx);
    %         results.pix2(end) = labeledPix(jdx);
    %         results.fileIx(end) = fileIx;
    %         results.file(end) = fn;
    % 
    %         results.distance(end) = sqrt((r1-r2)^2+(c1-c2)^2)/4;
    % 
    %         covs = nan(nXVals,1);
    %         corrs = nan(nXVals,1);
    %         for xvalIx = 1:nXVals
    %             pix1resp = squeeze(mean(trialTraces(idx,xval_splits(:,1,xvalIx),:,:),2,'omitnan'));
    %             pix2resp = squeeze(mean(trialTraces(jdx,xval_splits(:,2,xvalIx),:,:),2,'omitnan'));
    %             tmp = cov(pix1resp(:),pix2resp(:));
    %             covs(xvalIx) = tmp(1,2);
    %             corrs(xvalIx) = corr(pix1resp(:),pix2resp(:));
    %         end
    %         results.xval_cov(end) = mean(covs);
    %         results.xval_cor(end) = mean(corrs);
    %     end
    % end
end