function [trialTable, keepTrials] = verifyFiles(fn,dr, params)

load([dr filesep fn], 'trialTable');
nTrials = length(trialTable.trueTrialIx);
nDMDs = size(trialTable.filename,1);
keepTrials = true(nDMDs, nTrials);
for trialIx = nTrials:-1:1
    for DMDix = 1:nDMDs
        [~, tiffFn, ext] = fileparts(trialTable.fnRegDS{DMDix,trialIx}); 
        if isempty(dir([dr filesep tiffFn '.tif'])) & isempty(dir([dr filesep tiffFn '.h5']))
            disp(['Missing tiff or h5 file:' tiffFn])
            keepTrials(DMDix,trialIx) = false;
        else
            trialTable.fnRegDS{DMDix,trialIx} = [tiffFn, ext];
        end

        [~, alignFn] = fileparts(trialTable.fnAdata{DMDix,trialIx}); 
        if ~exist([dr filesep alignFn '.mat'], 'file')
            disp(['Missing alignData file:' alignFn])
            keepTrials(DMDix,trialIx) = false;
        else
            trialTable.fnAdata{DMDix,trialIx} = [alignFn '.mat'];
        end

        sourceFn = trialTable.filename{DMDix,trialIx};
        if ~exist([dr filesep sourceFn], 'file')
            disp(['Missing source data file:' sourceFn])
            keepTrials(DMDix, trialIx) = false;
        end
        [~,~,ext] = fileparts(sourceFn);
        if strcmpi(ext, '.dat')
            trialTable.fnRaw{DMDix,trialIx} = sourceFn;
        else
            rawFn = trialTable.fnRaw{DMDix,trialIx};
            if ~exist([dr filesep rawFn], 'file')
                disp(['Missing raw data file:' rawFn])
                keepTrials(DMDix, trialIx) = false;
            end
        end
    end
end
if ~all(keepTrials)
    disp(['Files were missing for ' int2str(sum(~keepTrials(:))) ' recordings; likely failed alignments. Proceeding without them.']);
end
if ~any(keepTrials)
    error('All trials were rejected due to missing alignment files!');
end
end