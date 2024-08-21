% dirs = uipickfiles(); %picking out session files 
%% Run code
bleaching = false; 

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
    % for DMDix = 1:2
    meanIMs_1 = squeeze(exptSummary.perTrialMeanIMs{1}(:,:,1, :)); %x,y,color,trial
    meanIMs_2 = squeeze(exptSummary.perTrialMeanIMs{2}(:,:,1, :)); %x,y,color,trial
    actIMs_1 = squeeze(exptSummary.perTrialActIms{1}(:,:,1, :)); %x,y,color,trial
    actIMs_2 = squeeze(exptSummary.perTrialActIms{2}(:,:,1, :)); %x,y,color,trial

    %normalize
    meanIMs_1 = normalizeLocal(meanIMs_1, 21, 0.9, false);
    meanIMs_2 = normalizeLocal(meanIMs_2, 21, 0.9, false);


    actIMs_1 = actIMs_1./prctile(actIMs_1(~isnan(actIMs_1)),99);
    actIMs_1 = actIMs_1 - prctile(actIMs_1(~isnan(actIMs_1)), 66);
    actIMs_1 = normalizeLocal(actIMs_1, 21, 1);

    actIMs_2 = actIMs_2./prctile(actIMs_2(~isnan(actIMs_2)),99);
    actIMs_2 = actIMs_2 - prctile(actIMs_2(~isnan(actIMs_2)), 66);
    actIMs_2 = normalizeLocal(actIMs_2, 21, 1);

    RGB1 = cat(4, max(0,actIMs_1), (max(0,meanIMs_1)), (max(0,meanIMs_1))); %removed sqrt transform
    RGB2 = cat(4, max(0,actIMs_2), (max(0,meanIMs_2)), (max(0,meanIMs_2))); %removed sqrt transform
    
    
    saveName = [s.folder filesep 'sessionMovie.tif'];
    % saveName2 = [s.folder filesep 'sessionMovie_DMD' int2str(2) '.tif'];
    
    %rescaling to fit both dmds in same figure
    rescaled = max( size(squeeze(RGB1(:,:,1,:))) , size(squeeze(RGB2(:,:,1,:))) ); %1 is h, 2, is w, 3 is depth (n/a)
    mins = min( size(squeeze(RGB1(:,:,1,:))) , size(squeeze(RGB2(:,:,1,:))) ); %1 is h, 2, is w, 3 is depth (n/a)
    validSlices = find(any(RGB1, [1 2 4]));
    % validSlices2 = find(any(RGB2, [1 2 4]));
    started = false;
    for sliceIx = 1:length(validSlices)
        slice = validSlices(sliceIx);
        dmd1IMG = squeeze(RGB1(:,:,slice,:)); dmd2IMG = squeeze(RGB2(:,:,slice,:));
        dmd1Padded = NaN(rescaled); dmd2Padded = NaN(rescaled);
        dmd1Padded(1:size(dmd1IMG,1), 1:size(dmd1IMG,2), :) = dmd1IMG;
        dmd2Padded(1:size(dmd2IMG,1), 1:size(dmd2IMG,2), :) = dmd2IMG;
        concatenatedImage = [dmd1Padded; dmd2Padded];
        % writingIM = insertText(squeeze(RGB(:,:,slice,:)), [5 5], sprintf('Trial %i', slice), FontSize=20, TextColor = "white");
        writingIM = insertText(concatenatedImage, [5 5], sprintf('Trial %i', slice), FontSize=20, TextColor = "white");
        if ~started
            imwrite(writingIM, saveName);
            started = true;
        else
            imwrite(writingIM, saveName, 'WriteMode', 'append');
        end
    end
% end
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

