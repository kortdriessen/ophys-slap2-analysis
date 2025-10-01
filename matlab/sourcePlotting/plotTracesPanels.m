function hFig = plotTracesPanels(traces, varargin)
% plotTracesPanels  Plot stacked traces by ROI across multiple panels.
%
% hFig = plotTracesPanels(traces, 'Channel',1, 'ROIs',1:size(traces,1), ...
%     'PerPanel',25, 'NumCols',6, 'Scale',1e4, 'Offset',15, ...
%     'Color',[0 0.447 0.741], 'Fs',[], 'Title','')
%
% INPUT
%   traces : [nROIs x nFrames] or [nROIs x nFrames x nChannels] numeric array
%
% Name-Value options (all optional)
%   'Channel'  : which channel to plot (default: 1, only if traces is 3D)
%   'ROIs'     : vector of ROI indices to include (default: 1:nROIs)
%   'PerPanel' : number of ROIs per panel (default: 25)
%   'NumCols'  : number of columns in the panel grid (default: 6)
%   'Scale'    : divide traces by this value (default: 1e4)
%   'Offset'   : vertical spacing between stacked traces (default: 15)
%   'Color'    : [r g b] for the trace (default: MATLAB blue)
%   'Fs'       : sampling rate in Hz; if provided, x-axis is time (s); else frames
%   'Title'    : figure title string
%
% OUTPUT
%   hFig : figure handle

% ------------ Parse inputs
validateattributes(traces, {'numeric'}, {'nonempty','ndims',2,3}, mfilename, 'traces', 1);
[nROI, nFrames, nChan] = size(traces);

if ndims(traces) == 3
    nChan = size(traces, 3);
end

p = inputParser;
addParameter(p,'Channel',1,@(x)isnumeric(x)&&isscalar(x)&&x>=1&&x<=nChan);
addParameter(p,'ROIs',1:nROI,@(x)isnumeric(x)&&isvector(x)&&all(x>=1)&all(x<=nROI));
addParameter(p,'PerPanel',25,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'NumCols',6,@(x)isnumeric(x)&&isscalar(x)&&x>=1);
addParameter(p,'Scale',1e4,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Offset',15,@(x)isnumeric(x)&&isscalar(x)&&x>0);
addParameter(p,'Color',[0 0.447 0.741],@(x)isnumeric(x)&&numel(x)==3);
addParameter(p,'Fs',[],@(x)(isnumeric(x)&&isscalar(x)&&x>0) || isempty(x));
addParameter(p,'Title','',@(x)ischar(x) || isstring(x));
parse(p,varargin{:});
opt = p.Results;

roiIdx  = opt.ROIs(:).';
nSel    = numel(roiIdx);
nPanels = ceil(nSel / opt.PerPanel);
nRows   = ceil(nPanels / opt.NumCols);

% X-axis: frames or time
if ~isempty(opt.Fs)
    x    = (0:nFrames-1) / opt.Fs;
    xlab = 'time (s)';
else
    x    = 1:nFrames;
    xlab = 'frame';
end

% ------------ Figure + layout
hFig = figure('Color','w');
tlo  = tiledlayout(nRows, opt.NumCols, 'TileSpacing','compact','Padding','compact');

titleStr = string(opt.Title);
if strlength(titleStr)==0
    titleStr = sprintf('Source Traces (Channel %d)', opt.Channel);
end
title(tlo, titleStr, 'FontWeight','bold');

% Fixed y-limits per panel so scales match across panels
yTop    = opt.Offset * 1.3;            % headroom above top trace
yBottom = -opt.Offset * opt.PerPanel;  % space for PerPanel traces

for pidx = 1:nPanels
    a = (pidx-1)*opt.PerPanel + 1;
    b = min(pidx*opt.PerPanel, nSel);
    chunk = roiIdx(a:b);
    nChunk = numel(chunk);

    ax = nexttile; hold(ax,'on');
    for k = 1:nChunk
        r = chunk(k);
        if ndims(traces) == 3
            y = squeeze(traces(r,:,opt.Channel)) ./ opt.Scale - opt.Offset*(k-1);
        else
            y = squeeze(traces(r,:)) ./ opt.Scale - opt.Offset*(k-1);
        end
        plot(ax, x, y, 'Color', opt.Color, 'LineWidth', 0.8);
    end

    yt = -opt.Offset*(nChunk-1) : opt.Offset : 0;
    set(ax, 'YLim', [yBottom yTop], ...
            'YTick', yt, ...
            'YTickLabel', string(fliplr(chunk)));

    xlim(ax, [x(1) x(end)]);
    xlabel(ax, xlab);
    ylabel(ax, 'Source #');
    box(ax,'off');
end
end
