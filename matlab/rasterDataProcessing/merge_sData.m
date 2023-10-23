function merged = merge_sData
stop = false;
dr = [];
epoch = 0;
while ~stop
    epoch = epoch+1;
    %(quest,dlgtitle,btn1,btn2,opts)
    stop = strcmpi('done selecting', questdlg(['Select Analyzed sData files (all DMDs) for Epoch ' int2str(epoch)], 'Merging sData', 'Select Files', 'Done Selecting', 'Select Files'));
    if stop
        break
    end
    [fns, dr] = uigetfile('*.mat', 'multiselect', 'on');
    if ~iscell(fns)
        if fns==0
            %user hit cancel, let them try again
            epoch = epoch-1;
            continue
        else
            fns = {fns}; %selected one file
        end
    end

    if length(fns)>1 %ensure DMDs are in order
        if contains(lower(fns{1}), 'dmd1')
            if ~contains(lower(fns{2}), 'dmd2')
                error('Filenames must include "DMD1" and "DMD2"')
            end
        elseif contains(lower(fns{2}), 'dmd1')
            if contains(lower(fns{1}), 'dmd2')
                fns = fns([2 1]);
            else
                error('Filenames must include "DMD1" and "DMD2"')
            end
        else
            error('Filenames must include "DMD1" and "DMD2"');
        end
        load([dr filesep fns{1}], 'sData');
        merged(epoch,1) = sData;
        load([dr filesep fns{2}], 'sData');
        merged(epoch,2) = sData;
    else
        load([dr filesep fns{1}], 'sData');
        merged(epoch,1) = sData;
    end
    mergedFilenames{epoch} = fns;
end

save([dr filesep 'sData_merged.mat'], 'merged', 'mergedFilenames')
end