function processVoltage
if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end

%dr = 'Z:\scratch\ophys\Adrian\dendritic vm\ASAP7y-JEDI2P comparison\slap2_786811_2025-02-25_11-20-21\cell1';
aParams = setParams('multiRoiRegSLAP2');

sParams.manualRois = true;
sParams.nParallelWorkers = 8;
sParams.analyzeHz = 1000;
sParams.baselineWindow_s = 1;
sParams = setParams('summarize_Voltage', sParams, true);

fullPathToTrialTable = [dr filesep 'trialTable.mat'];
%generate the trial table
if ~exist(fullPathToTrialTable, 'file')
    buildTrialTableSLAP2(dr);
end

%align files
multiRoiRegSLAP2(fullPathToTrialTable,aParams)

%summarize
summarize_Voltage(dr, sParams);
end