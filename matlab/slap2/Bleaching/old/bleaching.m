function [bleachingDataPerTrial, yFit] = bleachingAnalysis(dr, fn)
%pick expSummary file
if nargin<1
    [fn, dr] = uigetfile;
end
%grab first and last trial
experimentalData = load([dr filesep fn]);
data = experimentalData.exptSummary.E;
nTrials = size(data,1); %.F0(:,:,1)
%find first quality trial (bad trials show up as empty cells)
validTrials = (cellfun(@(x) isstruct(x),data(:)));
%convert trial F0 to 1D
bleachingDataPerTrial = zeros(1, nTrials);
firstTrialFlag = 0;
for trial = 1:nTrials
    if validTrials(trial) == 1
        firstTrialFlag = firstTrialFlag + 1;
        bleachingData = data{trial}.F0(:);
        meanBleached = mean(bleachingData, 'omitnan');
        if firstTrialFlag == 1
            figure(1)
            hist(bleachingData, 300)
            title('Raw Pixel Intensity of First and Last Trials')
            hold on
        end
    else
        bleachingData = nan;
        meanBleached = nan;
    end
    bleachingDataPerTrial(trial) = meanBleached;
end
%plot hist with means
hist(bleachingData,300)
hold off
%compare means
figure(2)
plot(bleachingDataPerTrial)
hold on
%linear fit to this
lengthOfFit = length(bleachingDataPerTrial(~isnan(bleachingDataPerTrial))) ;
[p,S] =  polyfit(1:lengthOfFit,bleachingDataPerTrial(~isnan(bleachingDataPerTrial)),1);
[yFit, ~] = polyval(p, 1:lengthOfFit,S);
title('Bleaching Decay Over Trials')
plot(yFit)
hold off
end