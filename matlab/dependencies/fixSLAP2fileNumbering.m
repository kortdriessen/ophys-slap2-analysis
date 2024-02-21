function fixSLAP2fileNumbering
%Addresses a bug in SLAP2 trial numbering as of Feb 2024, where trial
%numbers sometimes fail to increment. THis makes some files extra long,
% and subsequent trial numbers get out of sync

%We renumber the files and move files that don't have a matching partner
%to the QC folder. Some files remain extra-long. We will have to trim them
%in a future step (?);

%get a list of dat files in a given folder
dr = uigetdir;
files = dir([dr filesep '*.dat']);

%get time written from dat file timestamp

%for each DMD1 recording, find the corresponding DMD2 recording, if it
%exists
isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));

pat =  ['TRIAL' + digitsPattern(6)];
for fileIx = find(isDMD1)
    fn  = files(fileIx).name;
    DMD1trialNumStr=  extract(fn, pat);
    
    timeDiff = abs([files.datenum] - files(fileIx).datenum);
    minTime = min(timeDiff(~isDMD1));
    if minTime<1e-4
        DMD2fileIx = find((timeDiff<=minTime) & ~isDMD1,1,'first');
        DMD2trialNumStr=  extract(files(DMD2fileIx).name, pat);
        allDMD2files = dir([dr filesep '*DMD2*' DMD2trialNumStr '*']);

        renameTrialString = ['TRIAL' num2str(randi(999999), '%06u')];
        for fIx2 = 1:length(allDMD2files)
            oldName =  allDMD2files(fIx2).name;
            newName = replace(oldName,DMD2trialNumStr, DMD1trialNumStr); 
            if exist(newName, 'file')
                    rename = replace(newName,DMD2trialNumStr, renameTrialString); 
                    movefile([dr filesep newName], [dr filesep rename]);
            end
            movefile([dr filesep oldName], [dr filesep newName]);
        end
    else %no matching DMD2 file found
        %move the DMD1 files to QC folder
        if ~exist([dr filesep 'QC'], 'dir')
            mkdir([dr filesep 'QC']);
        end
        allDMD1files = dir([dr filesep '*DMD1*' DMD1trialNumStr '*']);
        for fIx2 = 1:length(allDMD1files)
            oldName =  allDMD1files(fIx2).name;
            movefile([dr filesep oldName], [dr filesep 'QC' filesep oldName]);
        end
    end
end


%move all unaccounted-for DMD2 files to QC
files = dir([dr filesep '*.dat']);
isDMD2 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD2'));
for fileIx = find(isDMD2)
    fn  = files(fileIx).name;
    DMD2trialNumStr=  extract(fn, pat);

    timeDiff = abs([files.datenum] - files(fileIx).datenum);
    minTime = min(timeDiff(~isDMD2));
    if minTime>=1e-4
        if ~exist([dr filesep 'QC'], 'dir')
            mkdir([dr filesep 'QC']);
        end
        allDMD1files = dir([dr filesep '*DMD1*' DMD2trialNumStr '*']);
        for fIx2 = 1:length(allDMD1files)
            oldName =  allDMD1files(fIx2).name;
            movefile([dr filesep oldName], [dr filesep 'QC' filesep oldName]);
        end
    end
end
    



