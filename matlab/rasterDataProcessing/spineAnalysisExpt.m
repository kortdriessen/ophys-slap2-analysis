classdef spineAnalysisExpt < handle
    %like spineanalysis, but reads data directly from .dat files
    properties
        hFigPanel; %image panel
        hButtonPanel; %button panel
        hF; %figure handle
        hGrid;
        hAx; %axes handle
        hButton1; hButton1c;
        hButton2;
        hButton3;
        hButton4;
        hButtonCh;
        hEditCh;

        hROIs;
        hIm; %image object

        chDim = 3; %images are in XYCT order
        timeDim = 4; %images are in XYCT order
        showChannel = 1;

        %butterworth filter parameters for highpass
        b2;
        a2;

        IM;
        IMc; %correlation image; normalized to median 0 standard deviation 1
        IMsk; %skewness image; normalized to median 0 standard deviation 1
        IMact; %the image used to display activity
        IMavg; %average image, normalized to [0 1]
        bleach;
        aData; %alignment data
        dsFac; %downsampling factor of the image file, not used for Dat files
        numChannels;
        fn;
        dr = [];
        drsave;
        fnsave;
        fnStem;

        exptData;

        defaultChLabels = {'iGluSnFr'; 'jRGECO'};
        Clevel = 3; %adjustable threshold for display of correlation image, in standard deviations
    end

    methods
        function obj = spineAnalysisExpt(arg1)
            if ~nargin || isempty(arg1)
                [obj.fn, obj.dr] = uigetfile('*.mat');
            else
                error('incorrect input!')
            end

            %load motion correction data
            obj.exptSummary = load([obj.dr filesep obj.fn], 'exptSummary');
                
            obj.numChannels = 1;
            obj.dsFac = nan;
            obj.IMavg = IMgamma;
            obj.IMc = obj.exptSummary.
            obj.IMact = 

            %gui creation
            obj.hF = uifigure('name', ['Spine analysis - ' obj.fn], 'KeyPressFcn', @obj.kpf);
            obj.hGrid = uigridlayout(obj.hF);
            obj.hGrid.RowHeight = {50, '1x'};
            obj.hGrid.ColumnWidth = {'1x'};

            obj.hButtonPanel = uigridlayout(obj.hGrid);
            obj.hButtonPanel.RowHeight = {'1x'};
            obj.hButtonPanel.ColumnWidth = {100,100,100,100,100,100, 100,100,100};
            obj.hButton1 = uibutton(obj.hButtonPanel, 'Text', 'Draw Polygon', 'ButtonPushedFcn', @obj.drawPolygon);
            obj.hButton1c = uibutton(obj.hButtonPanel, 'Text', 'Draw Circle', 'ButtonPushedFcn', @obj.drawCircle);
            obj.hButton2 = uibutton(obj.hButtonPanel, 'Text', 'Save ROIs', 'ButtonPushedFcn', @obj.saveROIs);
            obj.hButton3 = uibutton(obj.hButtonPanel, 'Text', 'Load ROIs', 'ButtonPushedFcn', @obj.loadROIs);
            obj.hButton4 = uibutton(obj.hButtonPanel, 'Text', 'Analyze', 'ButtonPushedFcn', @obj.analyze);
            for chIx = 1:obj.numChannels
                obj.hButtonCh(chIx) = uibutton(obj.hButtonPanel, 'Text', ['Show Ch' num2str(chIx) ':'], 'ButtonPushedFcn', @(x1,x2)(obj.setCh(chIx)));
                obj.hEditCh(chIx) = uieditfield(obj.hButtonPanel, 'Value', obj.defaultChLabels{chIx});
            end

            obj.hAx = uiaxes(obj.hGrid);
            obj.hIm = imshow(cat(3, obj.IMact(:,:,obj.showChannel)./obj.Clevel, obj.IMavg(:,:,obj.showChannel), obj.IMavg(:,:,obj.showChannel)),'parent', obj.hAx);
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

        function sData = analyze(obj, arg1, arg2)
            ind =strfind(obj.fn, '_REGISTERED');
            datFn = [obj.dr obj.fn(1:ind-1) '.dat'];
            S2data = slap2.Slap2DataFile(datFn);
            meta = loadMetadata(datFn);
            linerateHz = 1/meta.linePeriod_s;
            analyzeHz = 100; %framerate at which to analyze
            dt = ceil(linerateHz/analyzeHz);
            numChannels = S2data.numChannels;
            numLines = S2data.totalNumLines;
            numFrames = floor(numLines/dt);
            frameIxs = dt*(1:numFrames);
            
            if ~isempty(obj.aData) %populate motion predictors
                motionC = interp1(obj.aData.DSframes,obj.aData.motionDSc, frameIxs, 'linear', 'extrap');
                motionR = interp1(obj.aData.DSframes,obj.aData.motionDSr, frameIxs, 'linear', 'extrap');
                motPreds = [motionC ; motionR ; motionC.^2 ; motionR.^2 ; motionC.^3 ; motionR.^3 ; motionC.*motionR ; ones(size(motionC))];
            end

            %generate masks for the ROIs
            tix = 0;
            sData.names ={}; %INITIALIZE
            for rix = 1:length(obj.hROIs)
                try
                    roiObj = findobj(obj.hROIs(rix));
                catch
                    continue
                end
                tix = tix+1; %tix will count the total # of valid ROIs
                name = roiObj.Label;
                assert(~isempty(name), 'ROI names cannot be empty!');
                assert(~any(strcmpi(sData.names, name)), 'A ROI name was duplicated!')
                sData.names{tix} = name;

                sData.mask{tix} = roiObj.createMask;
                [rr,cc] = find(sData.mask{tix});
                sData.bbox{tix} = [min(cc) min(rr) ; max(cc) max(rr)];
                Dmask{tix} = sData.mask{tix}(sData.bbox{tix}(1,2):sData.bbox{tix}(2,2), sData.bbox{tix}(1,1):sData.bbox{tix}(2,1)); %#ok<AGROW> 
            end
            nROIs = tix;

            %load the high time resolution data
            %extize
            for rix =nROIs:-1:1
                D{rix} = nan([sData.bbox{rix}(2,2)-sData.bbox{rix}(1,2)+1, sData.bbox{rix}(2,1)-sData.bbox{rix}(1,1)+1, obj.numChannels, numFrames]);
            end
            disp('Reading High-res data')
            for fix = 1:numFrames
                if mod(fix, 100)==1
                    disp(['Frame ' int2str(fix) ' of ' int2str(numFrames) '...'])
                end
                for cix = 1:obj.numChannels
                    fdata = S2data.getImage(cix, ceil(fix*dt), 2*dt, 1); %Channel, Frame, Timewindow, Zindex
                    for rix = 1:nROIs
                        Rq = obj.aData.cropRow+motionR(fix)-1+(sData.bbox{rix}(1,2):sData.bbox{rix}(2,2));
                        Cq = obj.aData.cropCol+motionC(fix)-1+(sData.bbox{rix}(1,1):sData.bbox{rix}(2,1));
                        [Cq,Rq] = meshgrid(Cq,Rq);
                        D{rix}(:,:,cix,fix) = interp2(fdata,Cq,Rq); 
                    end
                end
            end

            sData.trace0 = nan(numFrames, nROIs, obj.numChannels);
            sData.trace1 = nan(numFrames, nROIs, obj.numChannels);
            sData.trace2 = nan(numFrames, nROIs, obj.numChannels);
            for rix = 1:nROIs
                Dsz = size(D{rix});
                for cix = 1:obj.numChannels

                    DD = double(reshape(D{rix}(:,:,cix,:), [], Dsz(4)));
                    nans = isnan(DD);
                    numVal = sum(~nans, 2);
                    [~,p] = sort(numVal,'descend');
                    r = (1:length(p))'; r(p) = r;

                    %select pixels and frames to include
                    selpix = (numVal./max(numVal))>(r./max(r) + 0.05);
                    selframes = ~any(nans(selpix,:),1);
                    selpix = ~any(nans(:,selframes),2);
                    DD = DD(selpix,selframes); %discard nans
                    if isempty(DD)
                        continue
                    end

                    trace0= sum(DD(Dmask{rix}(selpix),:),1); %trace0 is the SUM over all selected pixels

                    %remove motion noise in the singular components of the movie
                    mDD = mean(DD,2);
                    [U,S,V] = svds(DD-mDD, floor(sum(selpix)/10));
                    Vcorr = V;
                    for vix = 1:size(V,2)
                        [~, Vcorr(:, vix)] = obj.decorrelateMotion(V(:,vix), motPreds(:,selframes));
                    end
                    DDcorr = DD-(U*S*Vcorr'); %motion-corrected data: the original data minus the reconstructed motion artefact
                    trace1= mean(DDcorr(Dmask{rix}(selpix), :),1); % the decorrelated trace computed from the mask alone

                    % %perform NMF
                    % nfacs = max(2, floor(sum(selpix)/12));
                    % [W,H] = nnmf(double(DDcorr),nfacs);
                    % 
                    % %find NMF components spatially correlated to trace0, temporally correlated, and with positive skewness
                    % HP = filtfilt(obj.b2,obj.a2,H')';
                    % t0HP = filtfilt(obj.b2,obj.a2,trace0);
                    % Cspace = corr(W, Dmask{rix}(selpix));
                    % Ctime = corr(HP', t0HP'); %high frequency correlation
                    % SK = skewness(HP, 0, 2);
                    % sel = Cspace>0.1 & Ctime>0.1 & SK>0.2;
                    % trace2 = mean(W(:,sel)*H(sel,:),1);

                    %compile output
                    sData.trace0(selframes, rix,cix) = trace0; %raw ROI, without decorrelating motion variables
                    sData.n0(rix,cix) = estimatenoise(trace0); %noise for trace0
                    sData.trace1(selframes, rix,cix) = trace1; %raw ROI, with motion variables decorrelated
                    sData.n1(rix,cix) = estimatenoise(trace1); %noise for trace1
                    % sData.trace2(selframes, rix,cix) = trace2; %NMF-based traces
                    % sData.n2(rix,cix) = estimatenoise(trace2); %noise for trace2
                    % sData.contamination(rix, cix) = 1-corr(trace0', trace2'); %a measure of how contaminated the raw trace is by other factors/noise
                end
            end

            sData.frametime = 1/analyzeHz;
            sData.bleach = obj.bleach; %bleaching curve for this field of view, in downsampled time

            %plot the data
            colors = [0.3 0.8 0.1; 0.8 0 0.3; 0.25 0.25 0.5]; %ch1 yellow
            T = (1:size(sData.trace0,1))/analyzeHz; %time vector
            figure,
            for pix = 1:nROIs
                for cix = 1:obj.numChannels
                    hAx(pix,cix) = subplot(nROIs*obj.numChannels, 1, (pix-1)*obj.numChannels+cix);
                    tNorm = sData.trace1(:,pix,cix)./sqrt(sData.n1(pix,cix)); %normalized trace
                    plot(T, tNorm, 'color', sqrt(colors(cix,:))); hold on;
                    %plot(T, nanfastsmooth(tNorm, 7, 3, 0.5), 'color', colors(cix,:), 'linewidth', 2);
                    ylabel([sData.names{pix} ' - Ch' int2str(cix)]);
                    if pix<nROIs
                        set(hAx(pix,cix), 'xticklabel', [], 'tickdir', 'out');
                    end
                end
            end
            xlabel('time (s)');
            linkaxes(hAx, 'x');

            %save figure
            disp('saving figure')
            exportgraphics(obj.hAx, [obj.drsave filesep obj.fnStem '_FIGURE.pdf'], 'ContentType', 'vector');

            %save output
            if isempty(obj.drsave) || ~any(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_TRACES.h5']);
            
            if obj.fnsave %if user didn't cancel
                obj.saveAsH5([obj.drsave filesep obj.fnsave(1:end-4) '.h5'] , sData);
            end


        end

        function saveAsH5(obj, fname, sData)
            if exist(fname, 'file')
                resp = questdlg(['The file ' fname ' already exists. Overwrite?'], 'Overwrite?', 'Yes', 'No', 'No');
                switch resp
                    case 'Yes'
                        delete(fname);
                    case 'No'
                        return
                end
            end
            nROIs = length(sData.names);

            szTrace = [size(sData.trace0,1) obj.numChannels]; %size of each trace

            %saves spine analyses data in H5 format for easier python analysis
            for roiIx = 1:nROIs
                roiName = sData.names{roiIx};

                %create data
                h5create(fname,['/fluo/raw/' roiName], [szTrace(2) szTrace(1)], 'Datatype', 'single', 'ChunkSize', min([szTrace(2) szTrace(1)], [1 2000]), 'Deflate', 5); % fluorescence from this ROI
                
                % write data
                h5write(fname,['/fluo/raw/' roiName],permute(sData.trace0(:,roiIx,:), [3 1 2]));

                %write attributes
                h5writeatt(fname,['/fluo/raw/' roiName],'fs', 1/obj.aData.frametime); %sampling frequency
                h5writeatt(fname,['/fluo/raw/' roiName],'noise', sData.n0(roiIx,:)); %estimated noise, per channel
                chans = {};
                for chIx = 1:obj.numChannels
                    chans = cat(1, chans, {get(obj.hEditCh(chIx), 'Value')});
                end
                h5writeatt(fname,['/fluo/raw/' roiName],'chans', chans); %channel names
                h5writeatt(fname,['/fluo/raw/' roiName],'nPix', sum(sData.mask{roiIx}(:))); %channel names
            end

            %save global information
            h5create(fname,"/fluo/motionC",size(motionC), 'Datatype', 'single','ChunkSize', size(motionC), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/motionC",motionC);
            h5create(fname,"/fluo/motionR",size(motionR), 'Datatype', 'single','ChunkSize', size(motionR), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/motionC",motionR);
            h5create(fname,"/fluo/alignmentError",size(obj.aData.aError), 'Datatype', 'single','ChunkSize', size(obj.aData.aError), 'Deflate', 5); % XY movement
            h5write(fname,"/fluo/alignmentError",obj.aData.aError);
            
            h5writeatt(fname,"/fluo/motionR",'fs', analyzeHz);
            h5writeatt(fname,"/fluo/motionC",'fs', analyzeHz);
            h5writeatt(fname,"/fluo/alignmentError",'fs', analyzeHz);

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
                    roiData = [roiData {S}];
                catch
                    %this object didn't exist
                end
            end
            save([obj.drsave filesep obj.fnsave], 'roiData');
        end

        function [corrected, traceMot] = decorrelateMotion(obj, trace, motPreds)

            %highpass filter trace
            trace = trace(:) - smooth(trace, ceil(8/obj.aData.frametime), 'lowess');
            beta = mvregress(motPreds',trace);
            traceMot = motPreds'*beta;
            corrected = trace' - traceMot';
        end

        function loadROIs (obj, arg1, arg2)
            if ~(any(obj.drsave)) || isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);
            load([dr filesep fn], 'roiData');

            try
                delete(obj.hROIs);
            catch
            end
            obj.hROIs = [];
            for rix = 1:length(roiData)
                switch roiData{rix}.Type
                    case 'images.roi.ellipse'
                        obj.hROIs(rix) = images.roi.Ellipse(obj.hAx, 'AspectRatio', roiData{rix}.AspectRatio, 'FixedAspectRatio', roiData{rix}.FixedAspectRatio, 'Center', roiData{rix}.Center, 'SemiAxes', roiData{rix}.SemiAxes, 'RotationAngle', roiData{rix}.RotationAngle, 'Label', roiData{rix}.Label);
                    case 'images.roi.circle'
                        obj.hROIs(rix) = images.roi.Circle(obj.hAx,'Center', roiData{rix}.Center, 'Label', roiData{rix}.Label);
                    case 'images.roi.polygon'
                        obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'Position', roiData{rix}.Position, 'Radius', roiData{rix}.Radius, 'Label', roiData{rix}.Label);
                    otherwise
                        error('bad ROI data type');
                end

                fnms = fieldnames(roiData{rix});
                for ix =1:length(fnms)
                    try
                        set(obj.hROIs(rix), fnms{ix}, roiData{rix}.(fnms{ix}));
                    catch
                    end
                end


                cm = get(obj.hROIs(rix), 'ContextMenu');
                uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(rix)))); %add a 'Set Label' context menu
            end
        end

        function drawPolygon(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            obj.hROIs(L) = drawpolygon(obj.hAx);
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function drawCircle(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            obj.hROIs(L) = drawellipse(obj.hAx);
            set(obj.hROIs(L), 'FaceAlpha', 0, 'labelAlpha', 0 , 'labelTextColor', 'y', 'linewidth', 0.5, 'MarkerSize', 1);
            cm = get(obj.hROIs(L), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(L)))); %add a 'Set Label' context menu
            obj.setLabel(obj.hROIs(L));
        end

        function setLabel(obj, hROI)
            label = inputdlg('Set a label');
            set(hROI, 'Label', label{1});
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