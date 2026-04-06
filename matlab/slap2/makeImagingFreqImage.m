function freqIM = makeImagingFreqImage(metaFilePath)
% makeImagingFreqImage Create and display imaging frequency image from meta file.
%   freqIM = makeImagingFreqImage(metaFilePath) loads the specified .meta MAT
%   file (or prompts to select one if not provided), computes the average
%   imaging rate per DMD pixel, displays the image, and returns the matrix.
%
%   Input:
%     metaFilePath (optional) - full path to the .meta MAT file.
%
%   Output:
%     freqIM - matrix of imaging frequencies (Hz) for each DMD pixel.

if nargin < 1 || isempty(metaFilePath)
    [fname, fpath] = uigetfile({'*.meta','SLAP2 meta files (*.meta)'; '*.*','All Files'}, 'Select meta file');
    if isequal(fname,0)
        error('No meta file selected.');
    end
    metaFilePath = fullfile(fpath, fname);
end

s = load(metaFilePath, '-mat');

% Expect variables in file: AcquisitionContainer, dmdPixelsPerRow, dmdPixelsPerColumn
if isfield(s,'AcquisitionContainer')
    AcquisitionContainer = s.AcquisitionContainer;
else
    error('Loaded file does not contain AcquisitionContainer.');
end
if isfield(s,'dmdPixelsPerRow') && isfield(s,'dmdPixelsPerColumn')
    dmdPixelsPerRow = s.dmdPixelsPerRow;
    dmdPixelsPerColumn = s.dmdPixelsPerColumn;
else
    error('Loaded file must contain dmdPixelsPerRow and dmdPixelsPerColumn.');
end

linesPerCycle = AcquisitionContainer.ParsePlan.linesPerCycle;

freqIM = zeros(dmdPixelsPerRow, dmdPixelsPerColumn*2);
for i = 1:linesPerCycle
    freqIM(AcquisitionContainer.ParsePlan.acqParsePlan(i).superPixelID) = ...
        freqIM(AcquisitionContainer.ParsePlan.acqParsePlan(i).superPixelID) + 1;
end

freqIM = freqIM / double(linesPerCycle) * AcquisitionContainer.ParsePlan.lineRateHz;

% Apply pixel replacement map if present
if isfield(AcquisitionContainer.ParsePlan, 'pixelReplacementMaps') && ...
        ~isempty(AcquisitionContainer.ParsePlan.pixelReplacementMaps) && ...
        size(AcquisitionContainer.ParsePlan.pixelReplacementMaps{1},2) >= 2
    repl = AcquisitionContainer.ParsePlan.pixelReplacementMaps{1};
    freqIM(repl(:,1)) = freqIM(repl(:,2));
end

freqIM = freqIM(:,1:dmdPixelsPerColumn);

figure; imagesc(freqIM'); colormap jet; axis image off;
cb = colorbar();
ylabel(cb,'Avg. imaging rate (Hz)');

end