function BATCH_MRR()
% does multiRoiRegSLAP2 for all directories in the raw mirror that contain 'acq_'.

root = '/data/raw_mirror';

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


aParams.operator = 'KD';
aParams.nWorkers = 64;
aParams = setParams('multiRoiRegSLAP2', aParams);

for i = 1:length(directoriesToProcess)
    dr = directoriesToProcess{i};
    fullPathToTrialTable = [dr filesep 'trialTable.mat'];
    
    if ~exist(fullPathToTrialTable, 'file')
        buildTrialTableSLAP2(dr);
    end
    mrr_done = CHECK_FOR_MRR_COMPLETION(fullPathToTrialTable);
    if mrr_done
        fprintf('MRR already done for %s, skipping.\n', dr);
        continue;
    end
    try
        multiRoiRegSLAP2(fullPathToTrialTable, aParams)
    catch ME
        timestamp = char(datetime('now','Format','yyyyMMdd-HHmmss'));
        errFile = fullfile(dr, ['MRR_ERROR_' timestamp '.txt']);
        fid = fopen(errFile, 'w');
        fprintf(fid, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
        fclose(fid);
        fprintf('ERROR in multiRoiRegSLAP2 for %s — logged to %s\n', dr, errFile);
    end
end
