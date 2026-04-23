function trial_table_fix(trialTablePath)
%TRIAL_TABLE_FIX One-off repair for the known E4T145 lastLine overflow.
%
% Usage:
%   trial_table_fix('/data/bellatrix/exp_6/loc_K/all_acqs/trialTable.mat')

arguments
    trialTablePath (1,:) char
end

assert(isfile(trialTablePath), 'trialTable file not found: %s', trialTablePath);

S = load(trialTablePath, 'trialTable');
assert(isfield(S, 'trialTable'), ...
    'File does not contain a variable named ''trialTable'': %s', trialTablePath);

trialTable = S.trialTable;

targetTrial = 145;
expectedEpoch = 4;
expectedBadLastLine = 2099161;
fixedLastLine = 2099160;

expectedDMD1 = 'acq-4_20260401_142927_DMD1-CYCLE-000000.dat';
expectedDMD2 = 'acq-4_20260401_142927_DMD2-CYCLE-000000.dat';

assert(isfield(trialTable, 'lastLine'), 'trialTable is missing field ''lastLine''.');
assert(isfield(trialTable, 'epoch'), 'trialTable is missing field ''epoch''.');
assert(isfield(trialTable, 'filename'), 'trialTable is missing field ''filename''.');

assert(size(trialTable.lastLine, 2) >= targetTrial, ...
    'trialTable.lastLine has only %d trials; expected at least %d.', ...
    size(trialTable.lastLine, 2), targetTrial);

assert(numel(trialTable.epoch) >= targetTrial, ...
    'trialTable.epoch has only %d entries; expected at least %d.', ...
    numel(trialTable.epoch), targetTrial);

assert(size(trialTable.filename, 2) >= targetTrial, ...
    'trialTable.filename has only %d trials; expected at least %d.', ...
    size(trialTable.filename, 2), targetTrial);

assert(trialTable.epoch(targetTrial) == expectedEpoch, ...
    'Trial %d is epoch %d, expected epoch %d.', ...
    targetTrial, trialTable.epoch(targetTrial), expectedEpoch);

assert(strcmp(trialTable.filename{1, targetTrial}, expectedDMD1), ...
    'Trial %d DMD1 file mismatch. Got: %s', ...
    targetTrial, trialTable.filename{1, targetTrial});

assert(strcmp(trialTable.filename{2, targetTrial}, expectedDMD2), ...
    'Trial %d DMD2 file mismatch. Got: %s', ...
    targetTrial, trialTable.filename{2, targetTrial});

oldLastLine = trialTable.lastLine(:, targetTrial);

if all(oldLastLine == fixedLastLine)
    fprintf('trialTable already patched: trialTable.lastLine(:,%d) = %s\n', ...
        targetTrial, mat2str(oldLastLine.'));
    return
end

assert(all(oldLastLine == expectedBadLastLine), ...
    ['Unexpected lastLine values for trial %d. ' ...
     'Expected all entries to equal %d, got %s.'], ...
    targetTrial, expectedBadLastLine, mat2str(oldLastLine.'));

trialTable.lastLine(:, targetTrial) = fixedLastLine;

save(trialTablePath, 'trialTable');

fprintf('Patched %s\n', trialTablePath);
fprintf('trialTable.lastLine(:,%d): %s -> %s\n', ...
    targetTrial, mat2str(oldLastLine.'), ...
    mat2str(trialTable.lastLine(:, targetTrial).'));
end
