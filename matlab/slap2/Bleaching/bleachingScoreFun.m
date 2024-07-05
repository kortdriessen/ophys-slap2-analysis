function [output, nTrials] = bleachingScoreFun(experimentalData)

DMDix = 1;

%grab first and last trial
data = experimentalData.E;
nTrials = size(data,1);

%find first quality trial (bad trials show up as empty cells)
validTrials = (cellfun(@(x) isstruct(x),data(:, DMDix)));
bleachingDataPerTrial = zeros(1, nTrials);

firstValidTrial = find(validTrials,1,'first');
lastValidTrial = find(validTrials,1,'last');
bleachingData(:,1) = mean(data{firstValidTrial,DMDix}.F0(:,:,1),2, 'omitmissing');
bleachingData(:,2) = mean(data{lastValidTrial,DMDix}.F0(:,:,1),2, 'omitmissing');

meanF0perTrial = nan(1,nTrials);
meanF0perTrial(validTrials)= cellfun(@(x)(mean(x.F0(:,:,1),'all', 'omitmissing')), data(validTrials,DMDix));

% OUTPUT REGRESSION DATA
[b,bint] = regress(meanF0perTrial(validTrials)',[find(validTrials) ones(sum(validTrials),1)]);
output = b;

%OPTIONAL PLOTTING
%figure;
% plot([1 nTrials]', b'* [1 nTrials; 1 1])
% if bint(1)>0 & bint(2)>0
%     msgbox('the slope is significantly rising! Your ROIs are getting brighter, likely due to damage');
% elseif bint(1)<0 & bint(2)<0
%     msgbox('the slope is significantly falling! Your ROIs are getting dimmer, likely due to bleaching');
% end


end