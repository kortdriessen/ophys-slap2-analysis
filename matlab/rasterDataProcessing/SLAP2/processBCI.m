function processBCI

dr = uigetdir;

%generate the trial table
trialTable = buildTrialTable(dr)

%align files
multiRoiRegBCI(dr)