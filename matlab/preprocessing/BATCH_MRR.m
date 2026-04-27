function BATCH_MRR()
% does multiRoiRegSLAP2 for all directories in the raw mirror that contain 'acq_'.

root = '/data/raw_mirror';

% Walk directory tree and collect directories whose path contains 'acq_'
items = dir(fullfile(root, '**'));
items = items([items.isdir] & ~ismember({items.name}, {'.', '..'}));

directoriesToProcess = {};
if contains(root, 'acq_')
    directoriesToProcess{end+1} = root;
end
for i = 1:length(items)
    fullPath = fullfile(items(i).folder, items(i).name);
    if contains(fullPath, 'acq_')
        directoriesToProcess{end+1} = fullPath;
    end
end

for i = 1:length(directoriesToProcess)
    dr = directoriesToProcess{i};
    fullPathToTrialTable = [dr filesep 'trialTable.mat'];
    aParams.operator = 'KD';
    aParams = setParams('multiRoiRegSLAP2', aParams);
    if ~exist(fullPathToTrialTable, 'file')
        buildTrialTableSLAP2(dr);
    end
    multiRoiRegSLAP2(fullPathToTrialTable, aParams)
end
