function BATCH_LoCo_Pipeline()
% does multiRoiRegSLAP2 for all directories in the raw mirror that contain 'acq_'.

root = '/data/raw_mirror/miaplacidus';

% Walk directory tree and collect directories named 'acq_*' whose parent is named 'loc_*'
items = dir(fullfile(root, '**'));
items = items([items.isdir] & ~ismember({items.name}, {'.', '..'}));

directoriesToProcess = {};
for i = 1:length(items)
    [~, parentName] = fileparts(items(i).folder);
    if contains(items(i).name, 'acq_') && contains(parentName, 'loc_')
        directoriesToProcess{end+1} = fullfile(items(i).folder, items(i).name);
    end
end


sParams.operator = 'KD';
sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams = setParams('summarize_LoCo', sParams);
sParams.batchSize = 100;
sParams.nParallelWorkers = 8;

for i = 1:length(directoriesToProcess)
    dr = directoriesToProcess{i};
    fullPathToTrialTable = [dr filesep 'trialTable.mat'];
    
    mrr_done = CHECK_FOR_MRR_COMPLETION(fullPathToTrialTable);
    if ~mrr_done
        fprintf('MRR NOT done for %s, skipping.\n', dr);
        continue;
    end
    
    % check if any LoCo error file already exists for this directory
    errFiles = dir(fullfile(dr, 'LoCo_ERROR_*.txt'));
    if ~isempty(errFiles)
        fprintf('Previous LoCo error detected for %s, skipping. Found error file: %s\n', dr, errFiles(1).name);
        continue;
    end

    % check if an ExperimentSummary directory already exists in this directory, and if so, skip
    if exist(fullfile(dr, 'ExperimentSummary'), 'dir')
        fprintf('ExperimentSummary already exists for %s, skipping.\n', dr);
        continue;
    end


    try
        summarize_LoCo_ISOLATED(dr, sParams);
    catch ME
        timestamp = char(datetime('now','Format','yyyyMMdd-HHmmss'));
        errFile = fullfile(dr, ['LoCo_ERROR_' timestamp '.txt']);
        fid = fopen(errFile, 'w');
        fprintf(fid, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
        fprintf('ERROR in summarize_LoCo_ISOLATED for %s — logged to %s\n', dr, errFile);
    end
end
