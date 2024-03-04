dirs = uipickfiles();

dataDirs = {};
for i = 1:length(dirs)
    dataDirs = [dataDirs, findAllScanDirs(dirs{i})];
end


parfor idx = 1:length(dataDirs)
    [dFolder, dName] = fileparts(dataDirs{idx});
    disp(['Running directory ' dName])
    stripRegistrationBergamo_saveinplace([],fullfile(dFolder,dName,[dName '.tif']));
    
    tmp = dir(fullfile(dFolder,dName,'*_DOWNSAMPLED*'));
    downsampledFile = tmp.name;
    summarizeBergamo_Peaks(fullfile(dFolder,dName),downsampledFile);
end

%%
function scanDirs = findAllScanDirs(dirName)

scanDirs = {};

[~, dirNameOnly] = fileparts(dirName);
if (length(dirNameOnly) > 5) && all(dirNameOnly(1:5) == 'scan_')
    scanDirs = [scanDirs, dirName];
    return;
end

subDirs = dir(dirName);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name},{'.','..'}));

for idx = 1:length(subDirs)
    subDirPath = fullfile(dirName, subDirs(idx).name);
    scanDirs = [scanDirs, findAllScanDirs(subDirPath)];
end

end