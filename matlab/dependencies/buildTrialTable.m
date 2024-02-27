function trialTable = buildTrialTable(dr)
%Addresses a bug in SLAP2 trial numbering as of Feb 2024, where trial
%numbers sometimes fail to increment. THis makes some files extra long,
% and subsequent trial numbers get out of sync

%parameters
lineDiffThresh = 2000; %difference threshold for calling two recordings the same length, in lines. ~0.2 seconds

trialTable.DMD1filename = {};
trialTable.DMD1firstLine = [];
trialTable.DMD1lastLine = [];
trialTable.DMD2filename = {};
trialTable.DMD2firstLine = [];
trialTable.DMD2lastLine = [];
trialTable.trialEndTimeFromPC = [];
trialTable.trialStartTimeInferred = [];

%get a list of dat files in a given folder
if ~nargin
    dr = uigetdir;
end
unpickedfiles = dir([dr filesep '*.dat']);

epoch = 0;
while ~isempty(unpickedfiles)
    epoch = epoch+1;
    [indx,tf] = listdlg('ListString',{unpickedfiles.name}, 'PromptString',['Select files for EPOCH ' int2str(epoch)]);
    if ~tf
        return
    end
    epochfiles{epoch} = unpickedfiles(indx);
    unpickedfiles(indx) = [];
end

trueTrialIx = 0;
for eIx = 1:epoch %for each epoch
    %get number of lines from dat file timestamp
    files = epochfiles{eIx};
    disp(['Loading metadata from ' int2str(length(epochfiles{eIx})) ' DAT files...']);
    for fIx = length(files):-1:1
        hDat = slap2.Slap2DataFile([dr filesep files(fIx).name]);
        numLines(fIx) = hDat.totalNumLines;
    end
    linePeriod_s = hDat.hDataFile.metaData.linePeriod_s;
    disp('done loading ')
    isDMD1 = cellfun(@(x)(~isempty(x)), strfind({files.name}, 'DMD1'));

    DMD1files = files(isDMD1);
    [~, sortorder] = sort([DMD1files.datenum], 'ascend');
    DMD1files = DMD1files(sortorder);
    numLinesDMD1 = numLines(isDMD1); numLinesDMD1 = numLinesDMD1(sortorder);

    DMD2files = files(~isDMD1);
    [~, sortorder] = sort([DMD2files.datenum], 'ascend');
    DMD2files = DMD2files(sortorder);
    numLinesDMD2 = numLines(~isDMD1); numLinesDMD2 = numLinesDMD2(sortorder);

    lastDMD1fIx = 0;
    lastDMD2fIx = 0;
    accumLines = [0 0];
    while lastDMD2fIx<length(DMD2files) && lastDMD1fIx<length(DMD1files)

        %ocnfirm that the trials match up AT SOME POINT SOON
        cumLines1 = cumsum(numLinesDMD1(lastDMD1fIx+1: min(end, lastDMD1fIx+3)))-accumLines(1);
        cumLines2 = cumsum(numLinesDMD2(lastDMD2fIx+1: min(end, lastDMD2fIx+3)))-accumLines(2);
        assert(min(abs(cumLines1 - cumLines2'), [],'all')<lineDiffThresh, 'Error lining up trials across DMDs!')

        nLines = min(cumLines1(1), cumLines2(1));
        trueTrialIx = trueTrialIx+1;
        trialTable.DMD1filename{trueTrialIx} = DMD1files(lastDMD1fIx+1).name;
        trialTable.DMD2filename{trueTrialIx} = DMD2files(lastDMD2fIx+1).name;
        trialTable.DMD1firstLine(trueTrialIx) = accumLines(1)+1;
        trialTable.DMD2firstLine(trueTrialIx) = accumLines(2)+1;
        trialTable.DMD1lastLine(trueTrialIx) = accumLines(1)+nLines;
        trialTable.DMD2lastLine(trueTrialIx) = accumLines(2)+nLines;
        trialTable.trueTrialIx(trueTrialIx) = trueTrialIx;
        trialTable.epoch(trueTrialIx) = eIx;

        if abs(cumLines1(1)-cumLines2(1))<lineDiffThresh
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD2 file
            accumLines = [0 0]; %reset accumulated lines

        elseif cumLines1(1) < cumLines2(1)
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD1files(lastDMD1fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD1fIx = lastDMD1fIx+1; %we finished this DMD1 file
            accumLines(1) = 0;
            accumLines(2) = accumLines(2) + nLines;
        elseif cumLines2(1) < cumLines1(1)
            trialTable.trialEndTimeFromPC(trueTrialIx) = DMD2files(lastDMD2fIx+1).datenum;
            trialTable.trialStartTimeInferred(trueTrialIx) = trialTable.trialEndTimeFromPC(trueTrialIx) - datenum(duration(0,0,nLines.*linePeriod_s));

            lastDMD2fIx = lastDMD2fIx+1; %we finished this DMD1 file
            accumLines(1) = accumLines(1) + nLines;
            accumLines(2) = 0;
        end
    end
end

save([dr filesep 'trialTable'], 'trialTable');
end


