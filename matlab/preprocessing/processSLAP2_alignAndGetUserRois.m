function processSLAP2_alignAndGetUserRois(dr)
%PROCESSSLAP2_ALIGNANDGETUSERROIS Run SLAP2 alignment, then draw user ROIs.

disp(['Script started at: ' char(datetime("now"))])
if ~nargin
    dr = uigetdir; % neuron folder where scans are, not project folder
end

aParams = setParams('multiRoiRegSLAP2');
fullPathToTrialTable = [dr filesep 'trialTable.mat'];

% generate the trial table
if ~exist(fullPathToTrialTable, 'file')
    buildTrialTableSLAP2(dr);
end

% log header timestamps
log_header_timestamps(dr)

% align files
multiRoiRegSLAP2(fullPathToTrialTable, aParams)

% draw user ROIs from the aligned TIFFs, then stop
get_user_rois(dr);
end
