
doRegistration = true;
if doRegistration
    aParams.denoise20Hz = true;
    aParams.overwriteExisting = false;
    aParams.ds_time = 1;
    aParams.nWorkers = 1; %disable the parallel pool within stripRegBergamo
    aParams = setParams('stripRegBergamo', aParams, true);
end

sParams.microscope = "bergamo";
sParams.drawUserRois = false;
sParams.nParallelWorkers = 1;
sParams.nanThresh = 0.5;
sParams.activityChannel = 2;
sParams = setParams('summarize_NoLoCo', sParams, true);

dirs = uipickfiles();
dataDirs = {};
for i = 1:length(dirs)
    dataDirs = [dataDirs, findScanDirsRecursive(dirs{i})];
end

%initialize parallel pool with n workers
%parpool(10);

if doRegistration
    parfor idx = 1:length(dataDirs)
        [dFolder, dName] = fileparts(dataDirs{idx});
        disp(['batch processing directory ' dName])
        try
            %built a trial table
            tmp = dir([dataDirs{idx} filesep '*.tif'])
            %ignore registered or figure files

            fns = {tmp.name};
            fns = fns(~contains(fns, 'REGISTERED') & ~contains(fns, 'FIGURE'));
            stripRegBergamo(dataDirs{idx}, fns, aParams)
            summarize_NoLoCo(dataDirs{idx}, sParams);
    catch ME
        fprintf("%s: %s\n",dName, ME.message)
        end
    end
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

%%
function scanTiffs = findAllScanTiffs(dirName)

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

function scanDirs = findScanDirsRecursive(dirName)

scanDirs = {};

subDirs = dir(dirName);
subDirs = subDirs([subDirs.isdir] & ~ismember({subDirs.name},{'.','..'}));

for idx = 1:length(subDirs)
    subDirPath = fullfile(dirName, subDirs(idx).name);
    scanDirs = [scanDirs, findScanDirsRecursive(subDirPath)];
end

dd = dir([dirName filesep 'scan*.tif']);
containsTifs = ~isempty(dd);
if containsTifs
     scanDirs = [scanDirs dirName];
     return
end
end



