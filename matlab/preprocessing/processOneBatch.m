function processOneBatch(inputFile, batchFile, outputFile)
%PROCESSONEBATCH  Run a batch of trials in a subprocess for memory isolation.
%   processOneBatch(inputFile, batchFile, outputFile)
%
%   Called via: matlab -batch "processOneBatch('in.mat','batch.mat','out.mat')"
%
%   When this function returns and the subprocess exits, ALL memory
%   (MATLAB internal allocator + glibc arenas) is returned to the OS.
%   This is the only approach that guarantees bounded memory usage
%   regardless of recording length.

fprintf('=== SUBPROCESS START (PID %d) ===\n', feature('getpid'));

% Restore parent's MATLAB path so all functions are accessible
tmp = load(inputFile, 'matlabPath');
path(tmp.matlabPath);

% Load shared input data
data = load(inputFile, 'dr', 'fns', 'fls', 'els', 'selPix', 'sources', ...
    'discardFrames', 'alignData', 'meanAligned', 'motOutput', ...
    'roiData', 'params');

% Load batch-specific data
batch = load(batchFile, 'batchTrials');

fprintf('Processing %d trials...\n', numel(batch.batchTrials));

% Call the ORIGINAL unmodified processAllTrials_Async (13 params)
E = processAllTrials_Async(data.dr, data.fns, data.fls, data.els, ...
    data.selPix, data.sources, data.discardFrames, data.alignData, ...
    data.meanAligned, data.motOutput, data.roiData, ...
    batch.batchTrials, data.params);

% Save results
save(outputFile, 'E', '-v7.3');

fprintf('=== SUBPROCESS COMPLETE ===\n');
end
