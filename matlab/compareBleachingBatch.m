% dirs = uipickfiles(); %picking out session files 
%% Run code
bleaching = true; 

xx = cell(length(dirs),1);
yy = cell(length(dirs),1);
labels = cell(length(dirs),1);

for dirNum = 1:length(dirs)
    currentDr = dirs{dirNum};
    sessionDetails = dir(currentDr);
    if bleaching
        sessionID = 'bleachTest';
        variant=sessionDetails(3).name;
        bleachResults.variantName = variant;
    else
        sessionID = ['*Neuron*'];
        variant=sessionDetails(3).name;
    end
    disp(['Processing ' variant])

    %GET DATA
    s = dir([currentDr filesep sessionID filesep 'ExperimentSummary' filesep '*.mat']);
    exptSummary = load([s.folder filesep s.name]).exptSummary ;

    %ASSES DECAY
    [decaySlope, nTrials] = bleachingScoreFun(exptSummary);
    yPlot = decaySlope'* [1 nTrials; 1 1];
    xPlot = [1 nTrials]';
    xx{dirNum} = xPlot;
    yy{dirNum} = yPlot; 
    if variant(7) == 'v'
        thisLabel = variant(7:11);
    else
        thisLabel = 'iGlu3';
    end
    labels{dirNum} = thisLabel;

    % CREATE MOVIE FOR EACH DMD AND SAVE
    for DMDix = 1:2
        meanIMs = squeeze(exptSummary.perTrialMeanIMs{DMDix}(:,:,1, :)); %x,y,color,trial
        actIMs = squeeze(exptSummary.perTrialActIms{DMDix}(:,:,1, :)); %x,y,color,trial

        %normalize
        meanIMs = normalizeLocal(meanIMs, 21, 0.9, false);
        actIMs = actIMs./prctile(actIMs(~isnan(actIMs)),99);
        actIMs = actIMs - prctile(actIMs(~isnan(actIMs)), 66);
        actIMs = normalizeLocal(actIMs, 21, 1);
        RGB = cat(4, max(0,actIMs), (max(0,meanIMs)), (max(0,meanIMs))); %removing sqrt to better show bleaching (temp)
        saveName = [s.folder filesep 'sessionMovie_DMD' int2str(DMDix) '.tif'];

        validSlices = find(any(RGB, [1 2 4]));
        started = false;
        for sliceIx = 1:length(validSlices)
            slice = validSlices(sliceIx);
            writingIM = insertText(squeeze(RGB(:,:,slice,:)), [5 5], sprintf('Trial %i', slice), FontSize=20, TextColor = "white");
            if ~started
                imwrite(writingIM, saveName);
                started = true;
            else
                imwrite(writingIM, saveName, 'WriteMode', 'append');
            end
        end
    end
        disp(['Done processing ' variant])

end

%% ANALYZE BLEACHING WITH BOXPLOTS
variantUnderScrutiny = unique(labels);
slopesForVariants = cell(length(variantUnderScrutiny),1);
for vN = 1:length(variantUnderScrutiny)
    v = variantUnderScrutiny{vN};
    vSlopes = yy(cellfun(@(x) strcmp(x, v), labels));
    if length(vSlopes) >1 
        combinedDecays = cell2mat(vSlopes);
        initial = (combinedDecays(:,1));
        final = (combinedDecays(:,2));
        combinedDecays = (initial-final)./30;
    else
        combinedDecays = cell2mat(vSlopes);
        initial = combinedDecays(1);
        final = combinedDecays(2);
        combinedDecays = (initial - final)/30;
    end
    slopesForVariants{vN} = combinedDecays;
end
combinedData = cell2mat(slopesForVariants);
grouping = arrayfun( ...
    @(x) repmat(x, size(slopesForVariants{x}, 1), 1), ...
    1:numel(slopesForVariants), ...
    'UniformOutput',false ...
    );
grouping = cell2mat(grouping');
xx = [1:length(variantUnderScrutiny)];


%% PLOT DECAYS
figure('Name', 'F0 versus trial number')
hold on
boxplot(combinedData, grouping, 'Labels', variantUnderScrutiny)
% errorbar(xx, slopesForVariants, err, '.')
title('F0 Decay Across 30 Trials')
xlabel('Variant')
ylabel('F0 Decay')
ylim([-10 10])
xticklabels(variantUnderScrutiny)
grid on
hold off