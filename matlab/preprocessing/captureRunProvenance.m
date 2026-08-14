function provenance = captureRunProvenance(acquisitionDirectory, params)
%CAPTURERUNPROVENANCE Capture reproducibility metadata for a summary run.

scriptPath = mfilename('fullpath');
ophysRoot = fileparts(fileparts(fileparts(scriptPath)));
slapClassPath = which('slap2.Slap2DataFile');
packageMarker = [filesep '+slap2' filesep];
markerIx = strfind(slapClassPath, packageMarker);
assert(~isempty(markerIx), 'slap2:cannotLocateSlap2Repository', ...
    'Could not infer the SLAP2 repository from slap2.Slap2DataFile.');
slap2Root = slapClassPath(1:markerIx(1)-1);

provenance.runIdentifier = getParamOrDefault(params, ...
    'validationRunId', char(java.util.UUID.randomUUID));
provenance.designation = getParamOrDefault(params, ...
    'validationDesignation', 'unspecified');
[~, defaultAcquisitionId] = fileparts(acquisitionDirectory);
provenance.acquisitionIdentifier = getParamOrDefault(params, ...
    'validationAcquisitionId', defaultAcquisitionId);
provenance.executionCommand = getParamOrDefault(params, ...
    'validationExecutionCommand', 'not supplied');
provenance.stagedAcquisitionDirectory = acquisitionDirectory;
provenance.ophysGit = gitState(ophysRoot);
provenance.slap2Git = gitState(slap2Root);
provenance.matlabVersion = version;
provenance.matlabRelease = version('-release');
provenance.startTimestamp = params.startTime;
provenance.endTimestamp = '';
provenance.tiffScratchDirectory = getenv('SLAP2_TIFF_SCRATCH_DIR');
provenance.matlabTempDirectory = tempdir;
[~, hostname] = system('hostname');
provenance.hostname = strtrim(hostname);
[~, osDescription] = system('uname -a');
provenance.operatingSystem = strtrim(osDescription);
end

function value = getParamOrDefault(params, fieldName, defaultValue)
if isfield(params, fieldName) && ~isempty(params.(fieldName))
    value = params.(fieldName);
else
    value = defaultValue;
end
if isstring(value)
    value = char(value);
end
end

function state = gitState(repositoryPath)
state.repositoryPath = repositoryPath;
[headStatus, head] = system(sprintf('git -C "%s" rev-parse HEAD', repositoryPath));
assert(headStatus == 0, 'slap2:gitProvenanceFailed', ...
    'Could not resolve Git HEAD for %s.', repositoryPath);
state.commit = strtrim(head);
[dirtyStatus, dirtyOutput] = system(sprintf( ...
    'git -C "%s" status --porcelain --untracked-files=all', repositoryPath));
assert(dirtyStatus == 0, 'slap2:gitProvenanceFailed', ...
    'Could not resolve Git status for %s.', repositoryPath);
state.dirty = ~isempty(strtrim(dirtyOutput));
state.statusPorcelain = strtrim(dirtyOutput);
end
