classdef drawROIs_REVIEW < handle
    %DRAWROIS_REVIEW Variant of drawROIs used by review_user_rois.
    %   Same image/axes/draw behavior as drawROIs, but the bottom button
    %   row replaces "Done" with two new buttons:
    %       - "Load ANNOTATIONS": pick any pipeline-style ANNOTATIONS.mat
    %         and populate this window's ROIs from the matching DMD entry.
    %       - "Resave All": back up any existing ANNOTATIONS.mat in the
    %         working directory, then write a new one capturing every
    %         currently-open review window.
    %   Both new buttons have their ButtonPushedFcn callbacks wired
    %   externally by review_user_rois (which is the only thing that knows
    %   about the full set of windows and the working directory).
    %
    %   Closing the windows without pressing "Resave All" leaves all files
    %   on disk untouched.

    properties
        hFigPanel; %image panel
        hButtonPanel; %button panel
        hF; %figure handle
        hGrid;
        hAx; %axes handle
        hButton1; hButton1c;
        hButton2;
        hButton3;
        hButtonLoadAnn; %Load ANNOTATIONS button (callback set by review_user_rois)
        hButtonResave;  %Resave All button   (callback set by review_user_rois)
        hButtonCh;
        hEditCh;

        hROIs;
        roiData;
        hIm;

        chDim = 3;
        showChannel = 1;

        IM;
        IMc;
        IMsk;
        IMact;
        IMavg;
        bleach;
        aData;
        dsFac;
        numChannels;
        fn;
        dr = [];
        drsave;
        fnsave;
        fnStem;
        DMDix; %which DMD index this window represents

        defaultChLabels = {'iGluSnFr'; 'jRGECO'};
        Clevel = 3;
    end

    methods
        function obj = drawROIs_REVIEW(IM, dr, fn, DMDix)
            obj.numChannels = size(IM,3);
            obj.dr = dr;
            obj.fn = fn;
            if nargin < 4 || isempty(DMDix)
                obj.DMDix = 1;
            else
                obj.DMDix = DMDix;
            end

            obj.hF = uifigure('name', ['ROI REVIEW - DMD ' num2str(obj.DMDix) ' - ' obj.fn], 'KeyPressFcn', @obj.kpf);
            obj.hGrid = uigridlayout(obj.hF);
            obj.hGrid.RowHeight = {50, '1x'};
            obj.hGrid.ColumnWidth = {'1x'};

            obj.hButtonPanel = uigridlayout(obj.hGrid);
            obj.hButtonPanel.RowHeight = {'1x'};
            obj.hButtonPanel.ColumnWidth = {120, 120, 120, 120, 170, 130};
            obj.hButton1     = uibutton(obj.hButtonPanel, 'Text', 'Draw Polygon',     'ButtonPushedFcn', @obj.drawPolygon);
            obj.hButton1c    = uibutton(obj.hButtonPanel, 'Text', 'Draw Circle',      'ButtonPushedFcn', @obj.drawCircle);
            obj.hButton2     = uibutton(obj.hButtonPanel, 'Text', 'Save ROIs',        'ButtonPushedFcn', @obj.saveROIs);
            obj.hButton3     = uibutton(obj.hButtonPanel, 'Text', 'Load ROIs',        'ButtonPushedFcn', @obj.loadROIs);
            obj.hButtonLoadAnn = uibutton(obj.hButtonPanel, 'Text', 'Load ANNOTATIONS', 'BackgroundColor', [0.85 0.92 1.00]);
            obj.hButtonResave  = uibutton(obj.hButtonPanel, 'Text', 'Resave All',       'BackgroundColor', [1.00 0.85 0.85]);

            obj.hAx = uiaxes(obj.hGrid);

            switch obj.numChannels
                case 1
                    IMdisp = IM./prctile(IM(:), 99.9);
                case 2
                    IM1 = IM(:,:,1)./prctile(IM(:,:,1), 99.9, 'all');
                    IM2 = IM(:,:,2)./prctile(IM(:,:,2), 99.9, 'all');
                    IMdisp = cat(3, IM1, IM2, IM2);
                otherwise
                    error('Unexpected number of channels (dim3) for drawROIs_REVIEW GUI');
            end

            obj.hIm = imshow(IMdisp,'parent', obj.hAx);
            set(obj.hIm, 'HitTest', 'off', 'pickableparts', 'none', 'clipping', 'on')
            set(obj.hAx, 'hittest', 'off', 'PickableParts', 'all')
        end

        function setCh(obj,chIx)
            obj.showChannel = chIx;
            updateDisplay(obj);
        end

        function updateDisplay(obj)
            set(obj.hIm, 'CData',cat(3, obj.IMact(:,:,obj.showChannel)./obj.Clevel, obj.IMavg(:,:,obj.showChannel), obj.IMavg(:,:,obj.showChannel)))
        end

        function saveROIs (obj, arg1, arg2)
            if isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_ROIs.mat']);
            if isempty(obj.fnsave) || ~any(obj.fnsave)
                return
            end
            roiData = {};
            for rix = 1:length(obj.hROIs)
                try
                    S = rmfield(get(obj.hROIs(rix)), {'Parent', 'Children', 'ContextMenu'});
                    roiData = [roiData {S}]; %#ok<AGROW>
                catch
                    %this object didn't exist
                end
            end
            obj.roiData = roiData;
            save([obj.drsave filesep obj.fnsave], 'roiData');

            disp('saving figure')
            exportgraphics(obj.hAx, [obj.drsave filesep obj.fnStem '_FIGURE.pdf'], 'ContentType', 'vector');
        end

        function loadROIsFromCell(obj, roiData)
            %LOADROISFROMCELL Populate the axes from a roiData cell array.
            %   Mirrors drawROIs.loadROIs but takes the cell array directly
            %   (no file I/O), and reattaches the right-click "Set Label"
            %   context menu so loaded ROIs can be relabeled.
            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            for rix = 1:length(roiData)
                % Build the ROI handle into a local variable first. Then
                % seed obj.hROIs with a direct assignment on the first ROI
                % (matching the drawPolygon / drawCircle pattern). This is
                % important: subscripted assignment into the empty double
                % `obj.hROIs = []` silently demotes the handle to its
                % numeric backing handle, which displays fine but later
                % breaks `createMask` (and any other method-style call) at
                % save time.
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('ellipse was incomplete, creating a dud ellipse to delete later')
                            newROI = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');
                        else
                            newROI = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.circle'
                        if isempty(roiData{rix}.Center)
                            disp(rix)
                            disp('circle was incomplete, creating a dud circle to delete later')
                            newROI = images.roi.Circle(obj.hAx,'Center', [5, 5], 'Label', 'dud');
                        else
                            newROI = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Label', roiData{rix}.Label);
                        end
                    case 'images.roi.polygon'
                        if isempty(roiData{rix}.Position)
                            disp(rix)
                            disp('Polygon was incomplete, creating a dud ellipse to delete later')
                            newROI = images.roi.Ellipse(obj.hAx, 'AspectRatio', 0.9, 'FixedAspectRatio', 0, 'Center', [5, 5], 'SemiAxes', [3, 3], 'RotationAngle', 215, 'Label', 'dud');
                        else
                            newROI = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Label', roiData{rix}.Label);
                        end
                    otherwise
                        error('bad ROI data type');
                end

                if rix == 1
                    obj.hROIs = newROI;
                else
                    obj.hROIs(rix) = newROI;
                end

                fnms = fieldnames(roiData{rix});
                for ix = 1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end

                cm = get(obj.hROIs(rix), 'ContextMenu');
                if ~isempty(cm)
                    uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(rix))));
                end
            end
        end

        function loadROIs (obj, arg1, arg2)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [fn, dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            if isequal(fn, 0)
                return
            end
            S = load([dr filesep fn], 'roiData');
            if ~isfield(S, 'roiData')
                warning('File does not contain ''roiData'': %s', [dr filesep fn]);
                return
            end
            obj.loadROIsFromCell(S.roiData);
        end

        function drawPolygon(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            if L==1
                obj.hROIs = drawpolygon(obj.hAx);
            else
                obj.hROIs(L) = drawpolygon(obj.hAx);
            end
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L))));
            obj.setLabel(obj.hROIs(L));
        end

        function drawCircle(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            if L==1
                obj.hROIs = drawellipse(obj.hAx);
            else
                obj.hROIs(L) = drawellipse(obj.hAx);
            end
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L))));
            obj.setLabel(obj.hROIs(L));
        end

        function setLabel(obj, hROI)
            label = inputdlg('Set a label');
            if ~isempty(label)
                set(hROI, 'Label', label{1});
            end
        end

        function kpf(obj, arg1, arg2)
            switch arg2.Key
                case 'a'
                    @obj.drawPolygon;
                case 'c'
                    @obj.drawCircle;
                case 'period'
                    obj.Clevel = obj.Clevel/1.5;
                    obj.updateDisplay
                case 'comma'
                    obj.Clevel = obj.Clevel.*1.5;
                    obj.updateDisplay;
            end
        end
    end
end
