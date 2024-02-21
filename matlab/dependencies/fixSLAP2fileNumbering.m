function trialTable = fixSLAP2fileNumbering(dr)
%Addresses a bug in SLAP2 trial numbering as of Feb 2024, where trial
%numbers sometimes fail to increment. THis makes some files extra long,
% and subsequent trial numbers get out of sync

%build a tabl
trialTable.DMD1filename = {};
trialTable.DMD1firstLine = [];
trialTable.DMD1lastLine = [];
trialTable.DMD2filename = {};
trialTable.DMD2firstLine = [];
trialTable.DMD2lastLine = [];
%DMD1filename;DMD1firstLine;DMD1lastLine;DMD2filename;DMD2firstLine;DMD2lastLine

%get a list of dat files in a given folder
if ~nargin
    dr = uigetdir;
end
files = dir([dr filesep '*.dat']);

%get number of lines from dat file timestamp
for fIx = 1:length(files)
    hDat = slap2.Slap2DataFile([dr filesep files(fIx).name]);
    numLines(fIx) = hDat.totalNumLines;
end
endTimes = [files.datenum];

isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));


% %for each DMD1 recording, find the corresponding DMD2 recording, if it
% %exists
% isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));
% 
% pat =  ['TRIAL' + digitsPattern(6)];
% for fileIx = find(isDMD1)
%     fn  = files(fileIx).name;
%     DMD1trialNumStr=  extract(fn, pat);
%     epochStr = extractBefore(fn,'_DMD');
% 
%     timeDiff = abs([ - files(fileIx).datenum);
%     numLinesDiff = abs(numLines -numLines(fileIx));
%     [minTime, ix1] = min(timeDiff(~isDMD1));
%     [minLines, ix2] = min(numLinesDiff(~isDMD1));
%     if minTime<1e-4
%         if ix1~=ix2
%             keyboard
%         end
%         DMD2fileIx = find((timeDiff<=minTime) & ~isDMD1,1,'first');
%         DMD2trialNumStr=  extract(files(DMD2fileIx).name, pat);
%         if strcmpi(DMD2trialNumStr,DMD1trialNumStr)
%             continue %go to next file
%         end
%         allDMD2files = dir([dr filesep epochStr '_DMD2*' DMD2trialNumStr{1} '*']);
% 
%         renameTrialString = ['TRIAL' num2str(randi(999999), '%06u')];
%         for fIx2 = 1:length(allDMD2files)
%             oldName =  allDMD2files(fIx2).name;
%             newName = replace(oldName,DMD2trialNumStr{1}, DMD1trialNumStr{1}); 
%             if exist(newName, 'file')
%                     rename = replace(newName,DMD2trialNumStr{1}, renameTrialString); 
%                     movefile([dr filesep newName], [dr filesep rename]);
%             end
%             movefile([dr filesep oldName], [dr filesep newName]);
%         end
%     else %no matching DMD2 file found
%         %move the DMD1 files to QC folder
%         if ~exist([dr filesep 'QC'], 'dir')
%             mkdir([dr filesep 'QC']);
%         end
%         allDMD1files = dir([dr filesep epochStr '_DMD1*' DMD1trialNumStr{1} '*']);
%         for fIx2 = 1:length(allDMD1files)
%             oldName =  allDMD1files(fIx2).name;
%             movefile([dr filesep oldName], [dr filesep 'QC' filesep oldName]);
%         end
%     end
% end
% 
% 
% %move all unaccounted-for DMD2 files to QC
% files = dir([dr filesep '*.dat']);
% isDMD2 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD2'));
% for fileIx = find(isDMD2)
%     fn  = files(fileIx).name;
%     DMD2trialNumStr=  extract(fn, pat);
% 
%     timeDiff = abs([files.datenum] - files(fileIx).datenum);
%     minTime = min(timeDiff(~isDMD2));
%     if minTime>=1e-4
%         if ~exist([dr filesep 'QC'], 'dir')
%             mkdir([dr filesep 'QC']);
%         end
%         allDMD1files = dir([dr filesep '*DMD1*' DMD2trialNumStr '*']);
%         for fIx2 = 1:length(allDMD1files)
%             oldName =  allDMD1files(fIx2).name;
%             movefile([dr filesep oldName], [dr filesep 'QC' filesep oldName]);
%         end
%     end
% end
    



