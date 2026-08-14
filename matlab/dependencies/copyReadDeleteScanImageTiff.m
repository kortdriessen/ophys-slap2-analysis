function [data, meta]= copyReadDeleteScanImageTiff(remotepath, localDir)
envLocalDir = getenv('SLAP2_TIFF_SCRATCH_DIR');
if (nargin < 2 || isempty(localDir)) && ~isempty(envLocalDir)
    localDir = envLocalDir;
end

if isempty(remotepath) %PASS EMPTY MATRIX TO INITIALIZE
    if nargin >= 2 && ~isempty(localDir) || ~isempty(envLocalDir)
        if ~exist(localDir, 'dir')
            mkdir(localDir);
        end
    else
        thisDir = fileparts(which('copyReadDeleteScanImageTiff'));
        if exist([thisDir filesep 'pathToFastDrive.mat'], 'file')
            return
        end
        try
            localDir = '/data/temp';
            if ~exist(localDir, 'dir')
                warning('User canceled setting local directory. defaulting to C:/temp on Windows else /scratch')
                if ispc
                    localDir = '/data/temp';
                else
                    localDir = '/data/temp';
                end
            end
        catch
            if ispc
                localDir = '/data/temp';
            else
                localDir = '/data/temp';
            end
        end
        save([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
    end
    return
end

if ~exist('localDir', 'var') || isempty(localDir)
    %read user settings for local directory
    thisDir = fileparts(which('copyReadDeleteScanImageTiff'));
    if exist([thisDir filesep 'pathToFastDrive.mat'], 'file')
        load([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
    else
        try
            localDir = uigetdir('', 'Select a local directory to save temporary files, on your FASTEST local drive!');
            if exist(localDir, 'dir')
               save([thisDir filesep 'pathToFastDrive.mat'], 'localDir');
            else
                warning('User canceled setting local directory. defaulting to C:/temp on Windows else /scratch')
                if ispc
                    localDir = '/data/temp';
                else
                    localDir = '/data/temp';
                end
            end
        catch
            if ispc
                localDir = '/data/temp';
            else
                localDir = '/data/temp';
            end
        end
    end
end
if ~exist(localDir, 'dir')
    mkdir(localDir);
end
disp(localDir)
assert(exist(remotepath, 'file'));
randName = [int2str(round(1e10*rand+0.1)) '.tif'];
localpath = [localDir  filesep   randName];
disp(localpath)
cleanupLocalCopy = onCleanup(@()deleteIfExists(localpath));


success = false; retries = 0;
while ~success
    try
        copyfile(remotepath, localpath);
        success = true;
    catch ME
        if retries<10 %we sometimes get system resource errors
            pause(5);
            retries = retries+1;
        else
            rethrow(ME)
        end
    end
end

[data, desc, meta] = networkScanImageTiffReader(localpath);
end

function deleteIfExists(path)
if exist(path, 'file')
    delete(path);
end
end
