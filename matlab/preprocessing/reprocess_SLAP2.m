function reprocess_SLAP2(dr, redoAlignment, sParams, aParams, noGUI)

if nargin<5
    noGUI = false; %for batch reprocessing
end

%find the trial table file
ttfn = [dr filesep 'trialTable.mat'];

%load the existing trial table
load(ttfn, 'trialTable')

if nargin<4 || isempty(aParams)
    aParams = setParams('multiRoiRegSLAP2', trialTable.alignParams, ~noGUI);
end
if nargin<3 || isempty(sParams)
    sParams.microscope = "SLAP2";
    sParams.drawUserRois = true;
    sParams.nParallelWorkers = 16; %relevant only for the initial reading segment
    sParams = setParams('summarize_LoCo', sParams, true);
end

if redoAlignment
    %delete all the alignmentdata and aligned data files
    for ix = 1:numel(trialTable.fnRegDS)
        delete(trialTable.fnRegDS{ix});
    end
    for ix = 1:numel(trialTable.fnAdata)
        delete(trialTable.fnAdata{ix});
    end
    trialTable = rmfield(trialTable, {'fnRegDS', 'fnAdata'});
    movefile(ttfn, [dr filesep 'ttOld-' datestr(now, 'YYYYMMDD') '.mat'])
    save(ttfn, 'trialTable');

    multiRoiRegSLAP2(ttfn,aParams)
end
summarize_LoCo(ttfn, sParams);