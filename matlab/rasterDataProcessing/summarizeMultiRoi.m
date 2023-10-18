function summarizeMultiRoi(nDMDs)
if ~nargin
    nDMDs = 2;
end

%select a set of aligned multiRoi recordings (trials) across two DMDs
[fns, dr] = uigetfile('*DMD1*REGISTERED*.tif', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

savedr = [dr filesep 'ExperimentSummary'];
if ~exist(savedr, 'dir')
    mkdir(savedr);
end
[fnsave, drsave] = uiputfile([savedr filesep 'Summary.mat'], 'Set filename for saving summary for this condition');

%confirm that all files exist for both DMDs
if nDMDs==2
for trialIx = length(fns):-1:1
        fn = fns{trialIx};
        dmdNumStringIndex = strfind(fn, 'DMD1')+3;
        fn(dmdNumStringIndex) = '2';
        assert(exist([dr filesep fn], 'file'), ['Missing aligned file:' fn]);
end
end

%load some metadata
fnStemEnd = strfind(fns{1}, '_REGISTERED') -1;
fnStem = fns{1}(1:fnStemEnd);
load([dr filesep fnStem '_ALIGNMENTDATA.mat'], 'aData');
numChannels = aData.numChannels;

%generate a concensus alignment across trials for further analysis
meanIM = cell(1,2);
for DMDix = nDMDs:-1:1
    meanIM{DMDix} = nan(1,1,2,1);
    actIM{DMDix} = nan(1,1,2,1);
    disp('Loading data and calculating correlation images...')
    for trialIx = length(fns):-1:1
        fn = fns{trialIx};
        
        dmdNumStringIndex = strfind(fn, 'DMD1')+3;
        fn(dmdNumStringIndex) = num2str(DMDix);

        %load the tiff
        A = ScanImageTiffReader([dr filesep fn]);
        IM = double(A.data);
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;
        meanIM{DMDix}(end:size(IM,1),:,:,:) = nan;
        meanIM{DMDix}(:, end:size(IM,2),:,:) = nan;
        meanIM{DMDix}(:,:,:,trialIx) = nan;
        meanIM{DMDix}(1:size(IM,1),1:size(IM,2),:,trialIx) = mean(IM,4, 'omitnan');

        %load alignment data
        fnStemEnd = strfind(fn, '_REGISTERED') -1;
        load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');

        %calculate correlation image
        [IMc, IMsk] = activityImage(IM, aData);
        actIM{DMDix}(end:size(IMc,1),:,:,:) = nan;
        actIM{DMDix}(:,end:size(IMc,2),:,:) = nan;
        actIM{DMDix}(:,:,:,trialIx) = nan;
        actIM{DMDix}(1:size(IMc,1),1:size(IMc,2),:,trialIx) = posNorm(IMc,IMsk);
    end

    %Make template
    disp('Making template for aligning across trials...')
    maxshift = 7;
    M = squeeze(sum(meanIM{DMDix}, 3));
    template = makeTemplateMultiRoi(M, maxshift); 
    
    %align all mean images to template
    disp('Aligning across trials...')
    Mpad = nan([size(template) size(M,3)]);
    Mpad(maxshift+(1:size(M,1)), maxshift+(1:size(M,2)),:) = M;
    for trialIx = length(fns):-1:1
        disp(['trial: ' fn])
        mot1 = xcorr2_nans(Mpad(:,:,trialIx), template, [0 ; 0], maxshift);
        [motOutput(:,trialIx, DMDix), corrCoeff(trialIx, DMDix)] = xcorr2_nans(Mpad(:,:,trialIx), template, round(mot1'), maxshift);
        [rr,cc] = ndgrid(1:size(meanIM{DMDix},1), 1:size(meanIM{DMDix},2));
        for chIx = 1:size(meanIM{DMDix},3)
            meanAligned{DMDix}(:,:,chIx,trialIx) = interp2(meanIM{DMDix}(:,:,chIx,trialIx), cc+motOutput(2,trialIx, DMDix), rr+motOutput(1,trialIx, DMDix));
            actAligned{DMDix}(:,:,chIx,trialIx) = interp2(actIM{DMDix}(:,:,chIx,trialIx), cc+motOutput(2,trialIx, DMDix), rr+motOutput(1,trialIx, DMDix));
        end
    end

    %identify outliers in alignment quality
    cc = corrCoeff(:, DMDix);
    corrThresh(DMDix) = min(min(0.99,median(cc)), mean(cc(cc>median(cc)) - 4*std(cc(cc>median(cc)))));
    selTrials(:, DMDix) = cc>corrThresh(DMDix);

    keyboard %for the second DMD aligned, we should use a template corresponding to
    %the selected trials from the first DMD!!
end

validTrials = all(selTrials,2);

%prepare file for saving
exptSummary.fns = fns;
exptSummary.dr = dr;
exptSummary.validTrials = validTrials;
%global images
exptSummary.meanIM = cellfun(@(x)(mean(x(:,:,:,validTrials),4, 'omitnan')), meanAligned, 'UniformOutput', false);
exptSummary.actIM = cellfun(@(x)(sqrt(mean(x(:,:,:,validTrials).^2,4, 'omitnan'))), actAligned, 'UniformOutput', false);
%per-trial images
exptSummary.perTrialMeanIMs = meanIM;
exptSummary.perTrialActIms = actIM;
exptSummary.perTrialAlignmentOffsets = motOutput; %the alignment vector for each trial

%save
save([drsave filesep fnsave], 'exptSummary');

%make some plots
%plot the peak correlation
figure('name', fnsave)
plot(corrCoeff(:, 1), 'r'); hold on, plot([1, length(fns)], corrThresh(1)*[1 1], ':r')
plot(corrCoeff(:, 2), 'b'); hold on, plot([1, length(fns)], corrThresh(2)*[1 1], ':b')
xlabel('Trials');
ylabel('alignment quality')

%plot the combined correlation image
Clevel = 2;
figure,
IMnorm = exptSummary.meanIM{1};
IMnorm = max(0, IMnorm./max(IMnorm, [], [1 2]));
subplot(4,1,1)
imshow(cat(3, exptSummary.actIM{1}(:,:,1)./Clevel, sqrt(IMnorm(:,:,1)), sqrt(IMnorm(:,:,1))));
subplot(4,1,2)
imshow(cat(3, exptSummary.actIM{1}(:,:,2)./Clevel, sqrt(IMnorm(:,:,2)), sqrt(IMnorm(:,:,2))));
IMnorm = exptSummary.meanIM{2};
IMnorm = max(0, IMnorm./max(IMnorm, [], [1 2]));
subplot(4,1,3)
imshow(cat(3, exptSummary.actIM{2}(:,:,1)./Clevel, sqrt(IMnorm(:,:,1)), sqrt(IMnorm(:,:,1))));
subplot(4,1,4)
imshow(cat(3, exptSummary.actIM{2}(:,:,2)./Clevel, sqrt(IMnorm(:,:,2)), sqrt(IMnorm(:,:,2))));

disp('Done summarizeMultiRoi')
end



function [IMc, IMsk]= activityImage(IM, aData)
            %compute correlation and skewness images
            [bb,aa] = butter(4, 0.05, 'high');
            IMgamma = sqrt(mean(IM,4, 'omitnan'));

            %select valid frames for cross-correlation analysis
            err = aData.aError;
            errStd = sqrt(estimatenoise(err));
            valid = err < (ordfilt2(err, 40, ones(1,200), 'symmetric')+5*errStd);

            noNan = double(IM(:,:,:,valid));
            noNan = noNan - mean(noNan,4, 'omitnan');
            nanInds = isnan(noNan);
            noNan(nanInds)=0;
            HP = permute(filtfilt(bb,aa,permute(noNan, [4 1 2 3])), [2 3 4 1]); clear noNan;
            HP(nanInds) = nan;

            %discard pixels with few measurements
            nanOut = sum(~isnan(HP(:,:,1,:)),4)<(size(HP,4)/3);
            HP(repmat(nanOut, [1 1 size(HP,3) size(HP,4)])) = nan;

            ss = sum(HP.^2,4, 'omitnan');
            vertC = sum(HP .* circshift(HP, [1 0 0 0]),4, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0 0]));
            horzC = sum(HP .* circshift(HP, [0 1 0 0]),4, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0 0]));

            C = mean(cat(4, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),4, 'omitnan');
            IMc = (C-median(C, [1 2], 'omitnan'))./std(C, 0,[1 2],'omitnan');
            
            HP = imgaussfilt(HP, [0.6 0.6]);
            sk = skewness(HP,0,4).*IMgamma;
            sknan = isnan(sk);
            sk(sknan)=1;
            sk = sk - imgaussfilt(sk, 10);
            sk(sknan) = nan;
            IMsk = sk./std(sk, 0,[1 2],'omitnan');
            %IMsk = (sk-median(sk(:), 'omitnan'))./std(sk, 0,'all','omitnan');

            %check if the bias in IMact is due to skewness or correlation
            %bias probably has to do with imaging rate at each pixel;
            %normalize to rate?
end

function IMout = posNorm(IM1, IM2)
    %computes the 2-norm of two images but only considers positive values
    IMout = sqrt(max(IM1,0).^2 + max(IM2,0).^2);
    IMout(isnan(IM1) & isnan(IM2)) = nan;
end