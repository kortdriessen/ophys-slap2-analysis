function [trialTable, keepTrials] = verifyFiles(dr, params)

switch params.microscope
    case 'SLAP2'
load([dr filesep 'trialTable.mat'], 'trialTable');
nTrials = length(trialTable.trueTrialIx);
keepTrials = true(nDMDs, nTrials);
for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        trialStr = ['E' int2str(trialTable.epoch(trialIx)) 'T' int2str(trialIx) 'DMD' int2str(DMDix)];
        tiffFn = [trialStr '_REGISTERED_DOWNSAMPLED-*Hz.tif'];
        if isempty(dir([dr filesep tiffFn]))
            disp(['Missing tiff file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        end
        alignFn =  [trialStr '_ALIGNMENTDATA.mat'];
        if ~exist([dr filesep alignFn], 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        end
        fnRaw{trialIx} = trialTable.(strcat('DMD', int2str(DMDix), 'filename')){trialIx};
        if ~exist([dr filesep fnRaw{trialIx}], 'file')
            disp(['Missing Dat file:' alignFn])
            keepTrials(DMDix, trialIx) = false;
        end
    end
end
if ~all(keepTrials)
    disp(['Files were missing for ' int2str(sum(~keepTrials(:))) ' [DMDs x trials]; likely failed alignments. Proceeding without them.']);
end
if ~any(keepTrials)
    error('All trials were rejected due to missing alignment files!');
end
    case 'bergamo'




end