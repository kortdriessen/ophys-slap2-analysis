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
    directory = directoriesToProcess{i};
    buildTrialTableSLAP2(directory);
end
