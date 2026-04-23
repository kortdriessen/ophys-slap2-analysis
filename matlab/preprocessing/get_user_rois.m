function ROIs = get_user_rois(dr_or_pathToTrialTable)
%GET_USER_ROIS Draw and save user ROIs from aligned SLAP2 TIFFs.
%   ROIs = GET_USER_ROIS(dr)
%   ROIs = GET_USER_ROIS(pathToTrialTable)
%
%   Uses the first valid registered TIFF for each DMD, matching the
%   existing summarize_* ROI workflow. If ANNOTATIONS.mat already exists,
%   the saved ROIs are loaded and returned without reopening the GUI.

if ~nargin
    [trialTablefn, dr] = uigetfile('*.mat', 'Select a trialTable file', '*trialTable*.mat');
else
    if exist(dr_or_pathToTrialTable, 'dir')
        dr = dr_or_pathToTrialTable;
        trialTablefn = 'trialTable.mat';
    else
        [dr, trialTablefn, ext] = fileparts(dr_or_pathToTrialTable);
        trialTablefn = [trialTablefn ext];
    end
end

copyReadDeleteScanImageTiff([]); % make sure the TIFF helper is initialized
[trialTable, keepTrials] = verifyFiles(trialTablefn, dr, struct());
nDMDs = size(trialTable.filename, 1);
fnAnn = [dr filesep 'ANNOTATIONS.mat'];

if exist(fnAnn, 'file')
    load(fnAnn, 'ROIs');
    disp(['Using existing user ROI annotations: ' fnAnn])
    return
end

for DMDix = 1:nDMDs
    firstValidTrial = find(keepTrials(DMDix,:), 1, "first");
    if isempty(firstValidTrial)
        error('No valid aligned trial was found for DMD %d; cannot draw user ROIs.', DMDix);
    end

    [~, fn, ext] = fileparts(trialTable.fnRegDS{DMDix, firstValidTrial});
    IM = copyReadDeleteScanImageTiff([dr filesep fn ext]);
    IM = squeeze(mean(IM, [3 4], 'omitnan'));
    hROIs(DMDix) = drawROIs(sqrt(max(0, IM)), dr, fn); %#ok<AGROW>
    ROIs(DMDix).dr = dr; %#ok<AGROW>
    ROIs(DMDix).fn = fn;
end

for DMDix = 1:nDMDs
    waitfor(hROIs(DMDix).hF);
    ROIs(DMDix).roiData = hROIs(DMDix).roiData;
end

save(fnAnn, 'ROIs');
disp(['Saved user ROI annotations: ' fnAnn])
end
