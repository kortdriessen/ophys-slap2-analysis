function processSLAP2(dr)
disp(['Script started at: ' char(datetime("now"))])
if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end
aParams = setParams('multiRoiRegSLAP2');
sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams.batchSize = 100;

%sParams.nParallelWorkers = 8;
sParams = setParams('summarize_LoCo', sParams, true);

fullPathToTrialTable = [dr filesep 'trialTable.mat'];
%generate the trial table
if ~exist(fullPathToTrialTable, 'file')
    buildTrialTableSLAP2(dr);
end

%log header timestamps
log_header_timestamps(dr)

%align files
multiRoiRegSLAP2(fullPathToTrialTable,aParams)

%summarize
summarize_LoCo_ISOLATED(dr, sParams);