function review_user_rois(dr_or_pathToTrialTable)
%REVIEW_USER_ROIS Open the SLAP2 ROI editor for an existing ANNOTATIONS.mat.
%   review_user_rois(dr)
%   review_user_rois(pathToTrialTable)
%
%   Opens one drawROIs_REVIEW window per DMD on the first valid registered
%   TIFF (matching the windows that get_user_rois / summarize_LoCo would
%   use). If an ANNOTATIONS.mat already lives in the working directory it
%   is auto-loaded into the windows so you can review it immediately.
%
%   Workflow:
%     - "Load ANNOTATIONS" lets you pick any pipeline-style ANNOTATIONS.mat
%       and replace the displayed ROIs with its contents (matched by DMD
%       index).
%     - Edit, move, rename (right-click -> Set Label), or delete ROIs as
%       needed in either window.
%     - "Resave All" backs up any existing ANNOTATIONS.mat in the working
%       directory (timestamped), then writes a new one capturing every
%       currently-open review window. The output is structured identically
%       to what summarize_LoCo / get_user_rois would have written, so it
%       drops in to the rest of the pipeline transparently.
%     - Closing the windows without pressing "Resave All" leaves all files
%       untouched.

if ~nargin
    [trialTablefn, dr] = uigetfile('*.mat', 'Select a trialTable file', '*trialTable*.mat');
    if isequal(trialTablefn, 0)
        return
    end
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

hROIs = drawROIs_REVIEW.empty;
for DMDix = 1:nDMDs
    firstValidTrial = find(keepTrials(DMDix,:), 1, "first");
    if isempty(firstValidTrial)
        error('No valid aligned trial was found for DMD %d; cannot open review GUI.', DMDix);
    end
    [~, fn, ext] = fileparts(trialTable.fnRegDS{DMDix, firstValidTrial});
    IM = copyReadDeleteScanImageTiff([dr filesep fn ext]);
    IM = squeeze(mean(IM, [3 4], 'omitnan'));
    hROIs(DMDix) = drawROIs_REVIEW(sqrt(max(0, IM)), dr, fn, DMDix); %#ok<AGROW>
end

% Auto-load the working-directory ANNOTATIONS.mat if present, so the user
% can immediately see the saved ROIs without having to click Load.
fnAnn = [dr filesep 'ANNOTATIONS.mat'];
if exist(fnAnn, 'file')
    populateFromFile(hROIs, fnAnn);
    disp(['Auto-loaded existing annotations: ' fnAnn])
end

for DMDix = 1:nDMDs
    hROIs(DMDix).hButtonLoadAnn.ButtonPushedFcn = @(~,~) onLoadAnnotations(hROIs);
    hROIs(DMDix).hButtonResave.ButtonPushedFcn  = @(~,~) onResaveAll(hROIs, dr);
end

end

% -------------------------------------------------------------------------

function onLoadAnnotations(hROIs)
%ONLOADANNOTATIONS Pick an ANNOTATIONS.mat and replace ROIs in all windows.
startDir = '';
for DMDix = 1:numel(hROIs)
    if isvalid(hROIs(DMDix)) && isvalid(hROIs(DMDix).hF)
        startDir = hROIs(DMDix).dr;
        break
    end
end
if isempty(startDir)
    return
end
[fn, dr] = uigetfile('*.mat', 'Select an ANNOTATIONS.mat file', ...
    [startDir filesep 'ANNOTATIONS.mat']);
if isequal(fn, 0)
    return
end
populateFromFile(hROIs, [dr filesep fn]);
end

function populateFromFile(hROIs, fullPath)
S = load(fullPath, 'ROIs');
if ~isfield(S, 'ROIs')
    warning('File does not contain a variable named ''ROIs'': %s', fullPath);
    return
end
loadedROIs = S.ROIs;
if numel(loadedROIs) ~= numel(hROIs)
    warning('Loaded ANNOTATIONS has %d DMD entries, current session has %d. Matching by index.', ...
        numel(loadedROIs), numel(hROIs));
end
for DMDix = 1:numel(hROIs)
    if ~isvalid(hROIs(DMDix)) || ~isvalid(hROIs(DMDix).hF)
        continue
    end
    if DMDix > numel(loadedROIs) || isempty(loadedROIs(DMDix).roiData)
        continue
    end
    hROIs(DMDix).loadROIsFromCell(loadedROIs(DMDix).roiData);
end
end

% -------------------------------------------------------------------------

function onResaveAll(hROIs, dr)
%ONRESAVEALL Save all open windows' ROIs to dr/ANNOTATIONS.mat.
%   Output is structured identically to what summarize_LoCo /
%   get_user_rois would have written: a struct array `ROIs` with fields
%   `dr`, `fn`, and `roiData` (cell array of structs from
%   `rmfield(get(roi), {Parent, Children, ContextMenu})` plus
%   `mask = createMask(roi)`).

fnAnn = [dr filesep 'ANNOTATIONS.mat'];

ROIs = struct('dr', {}, 'fn', {}, 'roiData', {});
anyCaptured = false;
for DMDix = 1:numel(hROIs)
    if ~isvalid(hROIs(DMDix)) || ~isvalid(hROIs(DMDix).hF)
        warning('Window for DMD %d is closed; it will not be included in the resave.', DMDix);
        continue
    end
    ROIs(DMDix).dr      = hROIs(DMDix).dr;
    ROIs(DMDix).fn      = hROIs(DMDix).fn;
    ROIs(DMDix).roiData = captureROIData(hROIs(DMDix).hROIs);
    anyCaptured = true;
end

if ~anyCaptured
    warning('No open review windows; nothing to save.');
    return
end

if exist(fnAnn, 'file')
    backup = [dr filesep 'ANNOTATIONS_backup_' ...
        char(datetime('now','Format','yyyyMMdd_HHmmss')) '.mat'];
    movefile(fnAnn, backup);
    disp(['Backed up previous annotations to: ' backup])
end

save(fnAnn, 'ROIs');
disp(['Saved annotations: ' fnAnn])
end

function roiData = captureROIData(hROIObjs)
%CAPTUREROIDATA Mirror drawROIs.close to produce the exact pipeline format.
%   For each ROI handle: take all observable properties via `get`, drop the
%   {Parent, Children, ContextMenu} fields (matching drawROIs.close), and
%   attach the binary mask via `createMask`. Failures are surfaced as
%   warnings rather than swallowed silently, so a regression here is
%   visible instead of producing an empty save like the original close().
roiData = {};
for rix = 1:length(hROIObjs)
    try
        S = rmfield(get(hROIObjs(rix)), {'Parent', 'Children', 'ContextMenu'});
        S.mask = createMask(hROIObjs(rix));
        roiData = [roiData {S}]; %#ok<AGROW>
    catch ME
        warning('Could not capture ROI %d (%s); it will not be saved.', rix, ME.message);
    end
end
end
