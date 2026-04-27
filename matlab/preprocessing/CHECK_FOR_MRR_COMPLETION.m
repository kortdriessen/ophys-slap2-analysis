function tf = CHECK_FOR_MRR_COMPLETION(fullPathToTrialTable)
% Returns true iff multiRoiRegSLAP2 finished cleanly on the given trialTable,
% as evidenced by trialTable.alignParams.endTime being set.

tf = false;
if ~exist(fullPathToTrialTable, 'file')
    return
end
S = load(fullPathToTrialTable, 'trialTable');
tf = isfield(S, 'trialTable') && isfield(S.trialTable, 'alignParams') ...
    && isfield(S.trialTable.alignParams, 'endTime') && ~isempty(S.trialTable.alignParams.endTime);
end
