dirs = uipickfiles();
%% Run code
bleaching = false;
analysisResults = cell(length(dirs), 1);
tic;
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
s = dir([currentDr filesep sessionID filesep 'ExperimentSummary' filesep '*.mat']);
exptSummary = load([s.folder filesep s.name]).exptSummary ;
currentOutput = bleachingAnalysis(s.folder, s.name);
bleachResults.bint = currentOutput;
bleachResults.visualShowcase = cell(1,2);
bleachResults.dr = s.folder;

for DMDix = 1:2
    meanIMs = squeeze(exptSummary.perTrialMeanIMs{DMDix}(:,:,1, :)); %x,y,color,trial
    actIMs = squeeze(exptSummary.perTrialActIms{DMDix}(:,:,1, :)); %x,y,color,trial

    %normalize
    meanIMs = normalizeLocal(meanIMs, 21, 0.9, false);
    %meanIMs = meanIMs./prctile(meanIMs(:),99);
    actIMs = actIMs./prctile(actIMs(~isnan(actIMs)),99);
    actIMs = actIMs - prctile(actIMs(~isnan(actIMs)), 66);
    actIMs = normalizeLocal(actIMs, 21, 1);
    RGB = cat(4, max(0,actIMs), sqrt(max(0,meanIMs)), sqrt(max(0,meanIMs)));
    bleachResults.visualShowcase{DMDix} = RGB;
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
    analysisResults{dirNum} = bleachResults;
    disp(['Done processing ' variant])
end
