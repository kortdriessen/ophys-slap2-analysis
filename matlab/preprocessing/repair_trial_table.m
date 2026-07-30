function info = repair_trial_table(trialTablePath, opts)
%REPAIR_TRIAL_TABLE  Clamp trial-table lastLine values that overrun the data.
%
% =========================================================================
% WHAT THIS UTILITY IS FOR
% =========================================================================
% It repairs a single, specific, recurring defect in a `trialTable.mat`
% that causes multiRoiRegSLAP2 (the MultiRoiReg / MRR step) to crash at the
% very end of an otherwise-successful registration with:
%
%     Error using slap2.util.MultiDataFiles/getLineDataIdxs (line 99)
%     Assertion failed.
%     Error in slap2.util.MultiDataFiles/getLineHeader (line 136)
%     Error in getOnlineMotion (line 19)
%     Error in multiRoiRegSLAP2>alignAsync (line 389)
%     Error in multiRoiRegSLAP2 (line 44)   [the parfor over trials]
%
% The failed assertion is `assert(all(cycleIdx <= obj.numCycles))` inside
% slap2.util.MultiDataFiles.getLineDataIdxs.
%
% -------------------------------------------------------------------------
% ROOT CAUSE (exactly what is wrong)
% -------------------------------------------------------------------------
% buildTrialTableSLAP2.m, in continuous-acquisition mode, carves each
% continuous recording into ~100k-line trials like this:
%
%     nLinesTot   = min(numLinesDMD1(fileIx), numLinesDMD2(fileIx)); % true #lines
%     trialEdges  = linspace(1, nLinesTot + 1, nTrialsInFile + 1);   % <-- the bug
%     ...
%     trialTable.lastLine(:, trueTrialIx) = trialEdges(trialIx + 1);
%
% Because the final edge is `nLinesTot + 1`, the LAST trial's `lastLine` is
% set to exactly ONE LINE PAST the end of the real data for whichever DMD
% has the FEWER total lines (the `min` DMD). The two DMDs are written with
% the SAME shared edge value, so the DMD with the larger line count still
% has margin and is fine; only the limiting DMD overflows.
%
% This overflow only materializes when the acquisition's final cycle-file is
% truncated mid-cycle (i.e. the recording was stopped manually rather than
% ending on a clean cycle boundary). The partial final cycle is correctly
% discarded by `floor()` in loadFileHeaderV2, so the DMD's true line count is
%
%     totalNumLines = numCycles * linesPerCycle
%
% and the trial table's `nLinesTot + 1` lands one line beyond it.
%
% -------------------------------------------------------------------------
% WHY IT SLIPS THROUGH REGISTRATION BUT KILLS getOnlineMotion
% -------------------------------------------------------------------------
% During the registration loop, frames are read through
% Slap2DataFile.getImageBatched, which CLAMPS the read window:
%
%     timeWindow = max(1, time-dtRead) : min(totalNumLines, time+dtRead);
%
% so an out-of-range line index silently reads the last valid line and no
% error is raised (all per-trial TIFFs get written). But AFTER the loop,
% alignAsync calls getOnlineMotion, which indexes line headers DIRECTLY with
% NO clamp:
%
%     cycleIdxs = ceil(lineIds / linesPerCycle);   % lineIds == DSframes
%     lineHeader = hDataFile.getLineHeader(resIdxs(idx), cycleIdxs(idx));
%
% For the overflowing final line, ceil(lastLine / linesPerCycle) == numCycles+1,
% which exceeds numCycles and trips the assertion. That asymmetry is why the
% run dies right at the finish line of a trial rather than during it.
%
% -------------------------------------------------------------------------
% WHAT THIS FUNCTION DOES
% -------------------------------------------------------------------------
% For every (DMD, trial) cell in the table it looks up the trial's data file,
% asks the authoritative slap2 loader for that file's true `totalNumLines`,
% and clamps any `lastLine` that exceeds it down to exactly `totalNumLines`
% (the largest valid line index). Setting lastLine == totalNumLines makes
%   ceil(lastLine / linesPerCycle) == numCycles  --> the assertion passes.
%
%   * Per-DMD: each DMD is clamped against ITS OWN line count, so the limiting
%     DMD is fixed while the other DMD's already-valid value is left alone.
%     (summarize_LoCo consumes firstLine/lastLine strictly per-DMD and never
%     requires the two DMDs to be equal, so a 1-line asymmetry is harmless.)
%   * Idempotent: a table with no overflow is reported clean and left
%     untouched; re-running after a fix is a no-op.
%   * Safe: before writing, it makes a one-time `*.prerepair.bak` copy of the
%     original table (never overwriting an existing backup).
%
% -------------------------------------------------------------------------
% WHAT THIS FUNCTION DELIBERATELY DOES NOT DO
% -------------------------------------------------------------------------
%   * It does NOT fix the source bug in buildTrialTableSLAP2.m. It only
%     repairs already-built tables. Fixing the generator is a separate change.
%   * It does NOT modify `firstLine`. A `firstLine` past the end of the data
%     is NOT the known off-by-one; if one is found, this function ABORTS with
%     an error rather than guessing, because that signals a different and more
%     serious problem (e.g. a mis-built table or wrong data files).
%   * It is NOT a general trial-table validator. Its single job is the
%     lastLine-overflow clamp described above.
%
% -------------------------------------------------------------------------
% WORKED EXAMPLES (real cases this was built from)
% -------------------------------------------------------------------------
%   miaplacidus/exp_1/loc_D/acq_1 : trial T539, DMD1
%       totalNumLines = 53,812,845 ; table lastLine = 53,812,846 (overflow 1)
%       -> clamped to 53,812,845
%   bellatrix/exp_6/loc_K (the original one-off this replaces) : trial E4T145
%       lastLine 2,099,161 -> 2,099,160 (overflow 1)
%
% -------------------------------------------------------------------------
% REQUIREMENTS
% -------------------------------------------------------------------------
% Must be run with the `slap2` toolbox on the MATLAB path (i.e. the same
% environment used to run the preprocessing pipeline / BATCH_MRR), because it
% uses slap2.Slap2DataFile to read each file's true line count. Opening the
% data files for a long continuous acquisition can take a few seconds; results
% are cached per unique filename so each file is opened only once.
%
% -------------------------------------------------------------------------
% USAGE
% -------------------------------------------------------------------------
%   % Preview without writing anything:
%   repair_trial_table('/data/raw_mirror/miaplacidus/exp_1/loc_D/acq_1/trialTable.mat', dryRun=true)
%
%   % Apply the fix (writes trialTable.mat, makes a .prerepair.bak first):
%   info = repair_trial_table('/data/raw_mirror/miaplacidus/exp_1/loc_D/acq_1/trialTable.mat');
%
% After applying, delete the acquisition's MRR_ERROR_*.txt (BATCH_MRR skips
% any directory that still has one) and re-run MRR. With overwriteExisting=false
% (the default) the re-run skips already-aligned trials and only finishes the
% unaligned ones plus the second DMD.
%
% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
%   info : struct with fields
%       .trialTablePath, .nDMD, .nTrials, .dryRun, .changed (logical),
%       .changes (struct array; one entry per clamped cell with dmd, trial,
%                 filename, oldLastLine, newLastLine, overflow, linesPerCycle)
%
% HISTORY
%   Supersedes the hard-coded `trial_table_fix.m` (which patched only the
%   bellatrix E4T145 / 2,099,161 case by literal value). This version detects
%   and clamps the overflow generically for any acquisition.
% =========================================================================

arguments
    trialTablePath (1,:) char
    opts.dryRun (1,1) logical = false
    opts.backup (1,1) logical = true
end

assert(isfile(trialTablePath), 'trialTable file not found: %s', trialTablePath);

S = load(trialTablePath, 'trialTable');
assert(isfield(S, 'trialTable'), ...
    'File does not contain a variable named ''trialTable'': %s', trialTablePath);
trialTable = S.trialTable;

requiredFields = ["filename", "firstLine", "lastLine"];
for f = requiredFields
    assert(isfield(trialTable, f), 'trialTable is missing required field ''%s''.', f);
end

dr = fileparts(trialTablePath);

[nDMD, nTrials] = size(trialTable.lastLine);
assert(isequal(size(trialTable.firstLine), [nDMD nTrials]), ...
    'trialTable.firstLine %s and trialTable.lastLine %s have mismatched sizes.', ...
    mat2str(size(trialTable.firstLine)), mat2str(size(trialTable.lastLine)));
assert(isequal(size(trialTable.filename), [nDMD nTrials]), ...
    'trialTable.filename %s and trialTable.lastLine %s have mismatched sizes.', ...
    mat2str(size(trialTable.filename)), mat2str(size(trialTable.lastLine)));

% Cache the authoritative line count (and linesPerCycle) per unique filename.
boundCache = containers.Map('KeyType', 'char', 'ValueType', 'double');
lpcCache   = containers.Map('KeyType', 'char', 'ValueType', 'double');

% Collect every cell whose lastLine runs past the real end of its data file.
changes = struct('dmd', {}, 'trial', {}, 'filename', {}, ...
                 'oldLastLine', {}, 'newLastLine', {}, 'overflow', {}, ...
                 'linesPerCycle', {});

for d = 1:nDMD
    for t = 1:nTrials
        fnCell = trialTable.filename{d, t};
        if isempty(fnCell)
            continue
        end

        [bound, lpc] = lineCountFor(dr, fnCell, boundCache, lpcCache);

        % A firstLine past the end is NOT the known off-by-one; refuse to guess.
        if trialTable.firstLine(d, t) > bound
            error('repair_trial_table:firstLineOverflow', ...
                ['Trial T%d (DMD%d) has firstLine=%g > totalNumLines=%d for file %s.\n' ...
                 'This is not the known lastLine off-by-one and may indicate a ' ...
                 'mis-built trial table or the wrong data files. Aborting without changes.'], ...
                t, d, trialTable.firstLine(d, t), bound, fnCell);
        end

        ll = trialTable.lastLine(d, t);
        if ll > bound
            changes(end + 1) = struct( ...
                'dmd', d, 'trial', t, 'filename', fnCell, ...
                'oldLastLine', ll, 'newLastLine', bound, ...
                'overflow', ll - bound, 'linesPerCycle', lpc); %#ok<AGROW>
        end
    end
end

info = struct();
info.trialTablePath = trialTablePath;
info.nDMD = nDMD;
info.nTrials = nTrials;
info.dryRun = opts.dryRun;
info.changes = changes;
info.changed = ~isempty(changes);

if isempty(changes)
    fprintf('[repair_trial_table] No lastLine overflow found - trial table already valid.\n  %s\n', ...
        trialTablePath);
    return
end

fprintf('[repair_trial_table] Found %d lastLine overflow(s) in:\n  %s\n', ...
    numel(changes), trialTablePath);
for k = 1:numel(changes)
    c = changes(k);
    flag = '';
    if c.overflow >= c.linesPerCycle
        flag = '   <-- WARNING: overflow >= one full cycle; verify this is the expected off-by-one';
    end
    fprintf('  DMD%d  T%d  (%s):  lastLine %g -> %d   [overflow %g line(s)]%s\n', ...
        c.dmd, c.trial, c.filename, c.oldLastLine, c.newLastLine, c.overflow, flag);
end

if opts.dryRun
    fprintf('[repair_trial_table] DRY RUN - nothing written. Re-run with dryRun=false to apply.\n');
    return
end

for k = 1:numel(changes)
    c = changes(k);
    trialTable.lastLine(c.dmd, c.trial) = c.newLastLine;
end

if opts.backup
    bakPath = [trialTablePath '.prerepair.bak'];
    if ~isfile(bakPath)
        copyfile(trialTablePath, bakPath);
        fprintf('[repair_trial_table] Backed up original to: %s\n', bakPath);
    else
        fprintf('[repair_trial_table] Backup already exists, left untouched: %s\n', bakPath);
    end
end

save(trialTablePath, 'trialTable');
fprintf('[repair_trial_table] Applied %d fix(es) and saved: %s\n', numel(changes), trialTablePath);

end


function [bound, lpc] = lineCountFor(dr, fnLocal, boundCache, lpcCache)
%LINECOUNTFOR  True totalNumLines + linesPerCycle for one data file (cached).
% containers.Map is a handle, so writing into the caches here updates the
% caller's copies too.
if isKey(boundCache, fnLocal)
    bound = boundCache(fnLocal);
    lpc   = lpcCache(fnLocal);
    return
end

fullPath = fullfile(dr, fnLocal);
assert(isfile(fullPath), ...
    'Referenced data file not found (needed to verify its line count): %s', fullPath);

try
    hDat = slap2.Slap2DataFile(fullPath);
catch ME
    error('repair_trial_table:dataFileOpenFailed', ...
        ['Failed to open %s via slap2.Slap2DataFile.\n' ...
         'This utility must run with the slap2 toolbox on the MATLAB path ' ...
         '(the same environment as the preprocessing pipeline).\n' ...
         'Underlying error: %s'], fullPath, ME.message);
end

hMulti = hDat.hMultiDataFiles;
bound  = double(hDat.totalNumLines);
lpc    = double(hMulti.header.linesPerCycle);
% Slap2DataFile defines no delete(); delete the contained MultiDataFiles so its
% DataFiles close their memmaps + MEX streams (delete cascades there).
delete(hMulti);

boundCache(fnLocal) = bound;
lpcCache(fnLocal)   = lpc;
end
