function optsOut = optionsGUI(opts, tooltips, funcName)
%OPTIONSGUI Simple editable options dialog with on-screen, scrollable layout.
%
%   optsOut = optionsGUI(opts, tooltips, funcName)
%
% This version keeps the original Save/Load/OK behavior, but fixes a common
% Windows display issue where the parameter window is created partly off-screen
% when many options are present. The dialog is now centered, constrained to the
% monitor height, resizable, and scrollable when the option list is taller than
% the visible window.

if nargin > 2
    caller = funcName;
else
    caller = dbstack;
    if length(caller) > 1
        caller = caller(2).name;
    else
        caller = 'Unknown Function';
    end
end

if nargin < 2 || isempty(tooltips)
    tooltips = struct();
end

optsOut = opts;

[optNames, ~] = sort(fieldnames(opts)); %#ok<ASGLU>
N = length(optNames);

% Layout constants.
titlesX = 8;
titlesW = 180;
etX = titlesX + titlesW + 20;
etW = 170;
resetW = 12;
H = 18;
H0 = 8;
rowStep = H + H0;
HOK = 40;
margin = 8;
sliderW = 18;
buttonH = 22;
buttonW = 55;

contentH = max(N * rowStep + 10, 40);
figW = titlesW + etW + resetW + 70 + sliderW;

% Keep the figure within the current screen. On Windows, fixed large figures
% can otherwise be created with the title bar off-screen and become impossible
% to move.
screen = get(groot, 'ScreenSize');  % [left bottom width height]
maxFigH = max(220, screen(4) - 140);
figH = min(contentH + HOK + 2 * margin, maxFigH);
figH = max(figH, 220);
figX = screen(1) + round((screen(3) - figW) / 2);
figY = screen(2) + round((screen(4) - figH) / 2);
figY = max(figY, screen(2) + 40);

scrollOffset = 0;
maxScrollOffset = 0;

handles = struct();
handles.F = figure( ...
    'Name', caller, ...
    'Units', 'pixels', ...
    'Position', [figX figY figW figH], ...
    'Toolbar', 'none', ...
    'Menubar', 'none', ...
    'Resize', 'on', ...
    'NumberTitle', 'off', ...
    'WindowStyle', 'normal', ...
    'ResizeFcn', @resizeGUI);

% Move on-screen after creation as an additional guard against multi-monitor
% or Windows display-scaling quirks.
drawnow;
movegui(handles.F, 'onscreen');

handles.panel = uipanel( ...
    'Parent', handles.F, ...
    'Units', 'pixels', ...
    'BorderType', 'none');

handles.scroll = uicontrol( ...
    'Units', 'pixels', ...
    'Parent', handles.F, ...
    'Style', 'slider', ...
    'Min', 0, ...
    'Max', 1, ...
    'Value', 0, ...
    'Visible', 'off', ...
    'Callback', @scrollCallback);

handles.load = uicontrol( ...
    'Units', 'pixels', ...
    'Parent', handles.F, ...
    'Style', 'pushbutton', ...
    'String', 'Load', ...
    'Callback', @loadButton);

handles.save = uicontrol( ...
    'Units', 'pixels', ...
    'Parent', handles.F, ...
    'Style', 'pushbutton', ...
    'String', 'Save', ...
    'Callback', @saveButton);

handles.OK = uicontrol( ...
    'Units', 'pixels', ...
    'Parent', handles.F, ...
    'Style', 'pushbutton', ...
    'String', 'OK', ...
    'Callback', @OK);

handles.titles = gobjects(1, N);
handles.reset = gobjects(1, N);
handles.ET = gobjects(1, N);

for n = 1:N
    handles.titles(n) = uicontrol( ...
        'Units', 'pixels', ...
        'Parent', handles.panel, ...
        'Style', 'text', ...
        'HorizontalAlignment', 'right', ...
        'String', optNames(n));

    handles.reset(n) = uicontrol( ...
        'Units', 'pixels', ...
        'Parent', handles.panel, ...
        'Style', 'pushbutton', ...
        'String', '', ...
        'Callback', @(varargin) reset(n));

    reset(n);
end

% Initialize layout after all controls exist. Start at the top of the list.
resizeGUI([], [], true);

waitfor(handles.F);

    function parseET(src, n)
        fieldName = optNames{n};
        type = class(opts.(fieldName));
        set(handles.titles(n), 'ForegroundColor', 'k')
        try
            switch type
                case 'logical'
                    strs = get(src, 'String');
                    val = get(src, 'Value');
                    if iscell(strs)
                        optsOut.(fieldName) = eval(strs{val});
                    else
                        optsOut.(fieldName) = eval(strtrim(strs(val, :)));
                    end
                case 'double'
                    strVal = get(src, 'String');
                    optsOut.(fieldName) = eval(['[' char(strVal) ']']);
                case 'cell'
                    strs = get(src, 'String');
                    val = get(src, 'Value');
                    if iscell(strs)
                        selected = strs{val};
                    else
                        selected = strtrim(strs(val, :));
                    end
                    try
                        optsOut.(fieldName) = eval(selected);
                    catch
                        optsOut.(fieldName) = selected;
                    end
                case 'char'
                    optsOut.(fieldName) = char(get(src, 'String'));
                case 'string'
                    optsOut.(fieldName) = string(get(src, 'String'));
                otherwise
                    error('Unsupported option type: %s', type);
            end
        catch ME %#ok<NASGU>
            % Mark the label red when parsing fails.
            set(handles.titles(n), 'ForegroundColor', 'r')
        end
    end

    function saveButton(varargin)
        for k = 1:length(handles.ET)
            if isgraphics(handles.ET(k))
                parseET(handles.ET(k), k)
            end
        end

        [fn, path] = uiputfile('*.mat');
        if isequal(fn, 0) || isequal(path, 0)
            return
        end
        save(fullfile(path, fn), 'optsOut');
    end

    function loadButton(varargin)
        [fn, path] = uigetfile('*.mat');
        if isequal(fn, 0) || isequal(path, 0)
            return
        end
        S = load(fullfile(path, fn), 'optsOut');
        if ~isfield(S, 'optsOut')
            warning('Selected file does not contain optsOut.');
            return
        end
        optsOut = S.optsOut;
        opts = optsOut;
        for n = 1:length(optNames)
            reset(n);
        end
        updateControlPositions();
    end

    function OK(varargin)
        for k = 1:length(handles.ET)
            if isgraphics(handles.ET(k))
                parseET(handles.ET(k), k)
            end
        end
        if isgraphics(handles.F)
            delete(handles.F)
        end
    end

    function reset(n)
        fieldName = optNames{n};
        set(handles.titles(n), 'ForegroundColor', 'k')
        optsOut.(fieldName) = opts.(fieldName);

        if length(handles.ET) >= n && isgraphics(handles.ET(n))
            delete(handles.ET(n));
        end

        switch class(opts.(fieldName))
            case 'logical'
                handles.ET(n) = uicontrol( ...
                    'Units', 'pixels', ...
                    'Parent', handles.panel, ...
                    'Style', 'popupmenu', ...
                    'String', {'false', 'true'}, ...
                    'Value', double(opts.(fieldName)) + 1, ...
                    'Callback', @(src, evnt) parseET(src, n)); %#ok<INUSD>
            case 'double'
                handles.ET(n) = uicontrol( ...
                    'Units', 'pixels', ...
                    'Parent', handles.panel, ...
                    'Style', 'edit', ...
                    'String', num2str(opts.(fieldName)), ...
                    'Callback', @(src, evnt) parseET(src, n)); %#ok<INUSD>
            case 'char'
                handles.ET(n) = uicontrol( ...
                    'Units', 'pixels', ...
                    'Parent', handles.panel, ...
                    'Style', 'edit', ...
                    'String', opts.(fieldName), ...
                    'Callback', @(src, evnt) parseET(src, n)); %#ok<INUSD>
            case 'string'
                handles.ET(n) = uicontrol( ...
                    'Units', 'pixels', ...
                    'Parent', handles.panel, ...
                    'Style', 'edit', ...
                    'String', char(opts.(fieldName)), ...
                    'Callback', @(src, evnt) parseET(src, n)); %#ok<INUSD>
            case 'cell'
                handles.ET(n) = uicontrol( ...
                    'Units', 'pixels', ...
                    'Parent', handles.panel, ...
                    'Style', 'popupmenu', ...
                    'String', opts.(fieldName), ...
                    'Value', 1, ...
                    'Callback', @(src, evnt) parseET(src, n)); %#ok<INUSD>
            otherwise
                error('Unsupported option type: %s', class(opts.(fieldName)));
        end

        if isfield(tooltips, fieldName)
            set(handles.titles(n), 'TooltipString', tooltips.(fieldName));
            set(handles.ET(n), 'TooltipString', tooltips.(fieldName));
        end
        updateControlPositions();
    end

    function resizeGUI(~, ~, initializeAtTop)
        if nargin < 3
            initializeAtTop = false;
        end
        if ~isfield(handles, 'F') || ~isgraphics(handles.F)
            return
        end

        figPos = get(handles.F, 'Position');
        figWNow = max(figPos(3), figW);
        figHNow = max(figPos(4), 180);

        % Avoid allowing resize to an unusably narrow window.
        if figWNow ~= figPos(3) || figHNow ~= figPos(4)
            set(handles.F, 'Position', [figPos(1), figPos(2), figWNow, figHNow]);
        end

        panelH = max(40, figHNow - HOK - margin);
        panelW = max(100, figWNow - 2 * margin - sliderW - 4);
        set(handles.panel, 'Position', [margin, HOK, panelW, panelH]);

        set(handles.save, 'Position', [margin, 8, buttonW, buttonH]);
        set(handles.load, 'Position', [margin + buttonW + 8, 8, buttonW, buttonH]);
        set(handles.OK, 'Position', [figWNow - margin - buttonW, 8, buttonW, buttonH]);

        maxScrollOffset = max(0, contentH - panelH);
        if maxScrollOffset > 0
            set(handles.scroll, ...
                'Visible', 'on', ...
                'Position', [figWNow - margin - sliderW, HOK, sliderW, panelH], ...
                'Min', 0, ...
                'Max', maxScrollOffset);
            if maxScrollOffset > 1
                set(handles.scroll, 'SliderStep', [min(1, rowStep / maxScrollOffset), min(1, panelH / maxScrollOffset)]);
            else
                set(handles.scroll, 'SliderStep', [1, 1]);
            end
            if initializeAtTop
                scrollOffset = maxScrollOffset;
            else
                scrollOffset = min(scrollOffset, maxScrollOffset);
            end
            set(handles.scroll, 'Value', scrollOffset);
        else
            scrollOffset = 0;
            set(handles.scroll, 'Visible', 'off', 'Value', 0, 'Min', 0, 'Max', 1);
        end

        updateControlPositions();
        drawnow;
        movegui(handles.F, 'onscreen');
    end

    function scrollCallback(src, ~)
        scrollOffset = get(src, 'Value');
        updateControlPositions();
    end

    function updateControlPositions()
        if ~isfield(handles, 'titles')
            return
        end
        for n = 1:N
            y = rowY(n);
            if isgraphics(handles.titles(n))
                set(handles.titles(n), 'Position', [titlesX, y - 3, titlesW, H]);
            end
            if isgraphics(handles.reset(n))
                set(handles.reset(n), 'Position', [etX + etW + 6, y - 3, resetW, H]);
            end
            if length(handles.ET) >= n && isgraphics(handles.ET(n))
                switch class(opts.(optNames{n}))
                    case {'logical', 'cell'}
                        set(handles.ET(n), 'Position', [etX, y, etW, H]);
                    otherwise
                        set(handles.ET(n), 'Position', [etX, y - 3, etW, H]);
                end
            end
        end
    end

    function y = rowY(n)
        % Display options top-to-bottom in sorted order. The scrollbar value is
        % defined as an upward offset: maxScrollOffset shows the top of content,
        % zero shows the bottom.
        baseY = contentH - n * rowStep + 5;
        y = baseY - scrollOffset;
    end
end
