function summarizeMultiRoi_Peaks(nDMDs)
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
    actIM{DMDix} = nan(1,1,1,1);
    disp('Loading data and calculating correlation images...')
    for trialIx = length(fns):-1:1
        fn = fns{trialIx};
        
        dmdNumStringIndex = strfind(fn, 'DMD1')+3;
        fn(dmdNumStringIndex) = num2str(DMDix);

        %load the tiff
        A = ScanImageTiffReader([dr filesep fn]);
        IM = double(A.data);
        if size(IM,3)<100
            error(['The file:' fn 'is very short. You should probably not include it?']);
        end
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []); %deinterleave;
        meanIM{DMDix}(end:size(IM,1),:,:,:) = nan;
        meanIM{DMDix}(:, end:size(IM,2),:,:) = nan;
        meanIM{DMDix}(:,:,:,trialIx) = nan;
        meanIM{DMDix}(1:size(IM,1),1:size(IM,2),:,trialIx) = mean(IM,4, 'omitnan');

        %load alignment data
        fnStemEnd = strfind(fn, '_REGISTERED') -1;
        load([dr filesep fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');

        %discard motion frames
        tmp = aData.aError-smoothdata(aData.aError,2, 'movmedian', ceil(2/aData.frametime));
        discardFrames = tmp>(3*std(tmp));
        IMf = squeeze(IM(:,:,1,:));
        IMf(:,:,discardFrames) = nan;

        [IMc, peaks{DMDix, trialIx}] = localizeFlashes(IMf);

        %calculate correlation image
        %[IMc, IMsk] = activityImage(IM, aData);
        actIM{DMDix}(end:size(IMc,1),:,:,:) = nan;
        actIM{DMDix}(:,end:size(IMc,2),:,:) = nan;
        actIM{DMDix}(:,:,:,trialIx) = nan;
        actIM{DMDix}(1:size(IMc,1),1:size(IMc,2),:,trialIx) = IMc; %posNorm(IMc,IMsk);
    end
end

%validTrials = all(selTrials,2);
validTrials = true(1, length(fns));

%prepare file for saving
exptSummary.fns = fns;
exptSummary.dr = dr;
exptSummary.validTrials = validTrials;
%global images
%exptSummary.meanIM = cellfun(@(x)(mean(x(:,:,:,validTrials),4, 'omitnan')), meanAligned, 'UniformOutput', false);
exptSummary.meanIM = cellfun(@(x)(mean(x(:,:,:,validTrials),4, 'omitnan')), meanIM, 'UniformOutput', false);
%exptSummary.actIM = cellfun(@(x)(sqrt(mean(x(:,:,:,validTrials).^2,4, 'omitnan'))), actAligned, 'UniformOutput', false);
exptSummary.actIM = cellfun(@(x)(sqrt(mean(x(:,:,:,validTrials).^2,4, 'omitnan'))), actIM, 'UniformOutput', false);

%per-trial images
exptSummary.peaks = peaks;
exptSummary.perTrialMeanIMs = meanIM;
exptSummary.perTrialActIms = actIM;
exptSummary.perTrialAlignmentOffsets = zeros(2,length(fns)); 
%exptSummary.perTrialAlignmentOffsets = motOutput; %the alignment vector for each trial

%save
save([drsave filesep fnsave], 'exptSummary');

%plot the combined correlation image
Clevel = 1000;
for DMDix = 1:2
figure,
IMnorm = exptSummary.meanIM{DMDix};
IMnorm = max(0, IMnorm./prctile(IMnorm(:), 99));
imshow(cat(3, sqrt(exptSummary.actIM{DMDix}(:,:,1))./Clevel, sqrt(IMnorm(:,:,1)), sqrt(IMnorm(:,:,1))));
colors = spring(size(peaks,2));
for trialIx = 1:size(peaks,2)
hold on
scatter(peaks{DMDix,trialIx}.col, peaks{DMDix,trialIx}.row, 5*peaks{DMDix,trialIx}.val.^2, 'MarkerEdgeColor',colors(trialIx,:));
end
end

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
            try
            HP = permute(filtfilt(bb,aa,permute(noNan, [4 1 2 3])), [2 3 4 1]); clear noNan;
            HP(nanInds) = nan;
            catch
                error(['A file is likely too short!']);
            end

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