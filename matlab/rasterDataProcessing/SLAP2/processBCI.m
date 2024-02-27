function processBCI

dr = uigetdir;

%generate the trial table
if ~exist([dr filesep 'trialTable.mat'], 'file')
    buildTrialTable(dr);
end

%align files
multiRoiRegBCI(dr)

%summarize
summarizeBCI(dr);