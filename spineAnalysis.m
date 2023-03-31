classdef spineAnalysis < handle
    properties
        hFigPanel; %image panel
        hButtonPanel; %button panel
        hF; %figure handle
        hGrid;
        hAx; %axes handle
        hButton1;
        hButton2;
        hButton3;
        hButton4;

        hROIs;
        hIm; %image object
       
        %butterworth filter parameters for highpass
        b2;
        a2;

        IM;
        IMc; %correlation image; normalized to median 0 standard deviation 1
        IMsk; %skewness image; normalized to median 0 standard deviation 1
        IMavg; %average image, normalized to [0 1]
        bleach;
        aData; %alignment data
        dsFac; %downsampling factor of the image file, used to properly set the framerate;
        fn;
        dr = [];
        drsave;
        fnsave;
        fnStem;

        Clevel = 3; %adjustable threshold for display of correlation image, in standard deviations
    end

    methods
        function obj = spineAnalysis(arg1)
            if ~nargin || isempty(arg1)
                [obj.fn, obj.dr] = uigetfile('*DOWNSAMPLED*.tif');
                A = ScanImageTiffReader([obj.dr obj.fn]);
                obj.IM = permute(A.data, [2 1 3]);
            elseif isstring(arg1)
                [obj.dr, obj.fn] = fileparts(arg1);
                A = ScanImageTiffReader([obj.dr obj.fn]);
                obj.IM = permute(A.data, [2 1 3]);
            elseif isnumeric(arg1) %treat input as an image
                obj.IM = arg1;
            else
                error('incorrect input!')
            end

            %load motion correction data
            fnStemEnd = strfind(obj.fn, '_REGISTERED') -1;
            obj.fnStem = obj.fn(1:fnStemEnd);
            try
                load([obj.dr filesep obj.fn(1:fnStemEnd) '_ALIGNMENTDATA.mat'], 'aData');
                obj.aData = aData;
            catch
                msgbox('Could not load alignment data; traces will not be decorrelated to motion.')
            end

            %infer downsampling from filename
            str = extractBetween(obj.fn, '_DOWNSAMPLED-', 'x.tif');
            if ~isempty(str)
                obj.dsFac = str2double(str{1});
            else
                obj.dsFac = 1;
            end


            %initialize butterworth filter
            [obj.b2,obj.a2] = butter(4, 0.05, 'high'); %highpass filter for generating correlation image

            %compute average image
            obj.IMavg = mean(obj.IM,3, 'omitnan');
            obj.IMavg = sqrt(max(0, obj.IMavg./prctile(obj.IMavg(:), 99.99)));

            %compute the bleaching curve
            selValid = ~any(isnan(obj.IM),3);
            tmp = reshape(obj.IM, [], size(obj.IM,3));
            obj.bleach = mean(tmp(selValid,:), 1, 'omitnan');

            %compute correlation image
            disp('computing correlation image');

            %select valid frames for cross-correlation analysis
            if ~isempty(obj.aData)
                err = aData.aError(1:obj.dsFac:end);
                errStd = sqrt(estimatenoise(err));
                valid = err < (ordfilt2(err, 40, ones(1,200), 'symmetric')+5*errStd);
            else
                valid = true(1, size(obj.IM,3));
            end

            
            noNan = double(obj.IM(:,:,valid));
            noNan = noNan - mean(noNan,3, 'omitnan');
            nanInds = isnan(noNan);
            noNan(nanInds)=0;
            HP = permute(filtfilt(obj.b2,obj.a2,permute(noNan, [3 1 2])), [2 3 1]); clear noNan;
            HP(nanInds) = nan;
            ss = sum(HP.^2,3, 'omitnan');
            vertC = sum(HP .* circshift(HP, [1 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0]));
            horzC = sum(HP .* circshift(HP, [0 1 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0]));
            C = mean(cat(3, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),3, 'omitnan');
            obj.IMc = (C-median(C(:), 'omitnan'))./std(C, 0,'all','omitnan');
            sk = skewness(HP,0,3);
            obj.IMsk = (sk-median(sk(:), 'omitnan'))./std(sk, 0,'all','omitnan');

            %gui creation
            obj.hF = uifigure('name', ['Spine analysis - ' obj.fn], 'KeyPressFcn', @obj.kpf);
            obj.hGrid = uigridlayout(obj.hF);
            obj.hGrid.RowHeight = {50, '1x'};
            obj.hGrid.ColumnWidth = {'1x'};
            
            obj.hButtonPanel = uigridlayout(obj.hGrid);
            obj.hButtonPanel.RowHeight = {'1x'};
            obj.hButtonPanel.ColumnWidth = {100,100,100,100};
            obj.hButton1 = uibutton(obj.hButtonPanel, 'Text', 'Add ROI', 'ButtonPushedFcn', @obj.drawPolygon);
            obj.hButton2 = uibutton(obj.hButtonPanel, 'Text', 'Save ROIs', 'ButtonPushedFcn', @obj.saveROIs);
            obj.hButton3 = uibutton(obj.hButtonPanel, 'Text', 'Load ROIs', 'ButtonPushedFcn', @obj.loadROIs);
            obj.hButton4 = uibutton(obj.hButtonPanel, 'Text', 'Analyze', 'ButtonPushedFcn', @obj.analyze);

            obj.hAx = uiaxes(obj.hGrid);
            obj.hIm = imshow(cat(3, obj.IMc./obj.Clevel, obj.IMavg, obj.IMavg),'parent', obj.hAx);
            set(obj.hIm, 'HitTest', 'off', 'pickableparts', 'none', 'clipping', 'on')
            set(obj.hAx, 'hittest', 'off', 'PickableParts', 'all')
            %imshow(cat(3, obj.IMc./obj.Clevel, obj.IMavg, obj.IMavg),'parent', obj.hAx);
        end

        function sData = analyze(obj, arg1, arg2)
            if ~isempty(obj.aData) %populate motion predictors
                motPreds = [obj.aData.motionC ; obj.aData.motionR ; obj.aData.motionC.^2 ; obj.aData.motionR.^2 ; obj.aData.motionC.^3 ; obj.aData.motionR.^3 ; obj.aData.motionC.*obj.aData.motionR ; ones(size(obj.aData.motionC))];
            end

            %generate masks for the ROIs
            tix = 0;
            for rix = 1:length(obj.hROIs)
                try
                    poly = findobj(obj.hROIs(rix));
                    Pos = poly.Position;
                catch
                    continue
                end
                if isempty(Pos)
                    continue;
                else
                    tix = tix+1; %tix will count the total # of valid ROIs
                    sData.poly{tix} = Pos;
                    sData.names{tix} = poly.Label;
                end
                sData.mask{tix} = poly.createMask;
                sData.bbox{tix} = [floor(min(Pos,[],1)) ; ceil(max(Pos,[],1))];
                Dmask{tix} = sData.mask{tix}(sData.bbox{tix}(1,2):sData.bbox{tix}(2,2), sData.bbox{tix}(1,1):sData.bbox{tix}(2,1));
            end
            nROIs = tix;

            %load the high time resolution data
            ind =strfind(obj.fn, '_DOWNSAMPLED');
            if ~isempty(ind)
                disp('Reading High-res data')
                fnRaw = [obj.fn(1:ind) 'RAW.tif'];
                nframes = obj.dsFac * size(obj.IM,3);
                t = Tiff([obj.dr filesep fnRaw], 'r');
                t.setDirectory(1);
            for fix = 1:nframes
                fdata = t.read;
                if fix==1
                    for rix =nROIs:-1:1
                        D{rix} = nan([sData.bbox{rix}(2,2)-sData.bbox{rix}(1,2)+1, sData.bbox{rix}(2,1)-sData.bbox{rix}(1,1)+1, nframes]);
                    end
                end
                for rix = 1:nROIs
                    D{rix}(:,:,fix) = fdata(sData.bbox{rix}(1,2):sData.bbox{rix}(2,2), sData.bbox{rix}(1,1):sData.bbox{rix}(2,1));
                end
                if fix<nframes
                    t.nextDirectory;
                end
            end
            end

            for rix = 1:nROIs
                Dsz = size(D{rix});
                DD = double(reshape(D{rix}, [], Dsz(3)));
                nans = isnan(DD);
                numVal = sum(~nans, 2);
                [~,p] = sort(numVal,'descend');
                r = (1:length(p))'; r(p) = r;

                %select pixels and frames to include
                selpix = (numVal./max(numVal))>(r./max(r) + 0.05);
                selframes = ~any(nans(selpix, :),1);
                selpix = ~any(nans(:,selframes),2);
                DD = DD(selpix,selframes); %discard nans
                trace0= mean(DD(Dmask{rix}(selpix), :),1);

                %remove motion noise in the singular components of the movie
                mDD = mean(DD,2);
                [U,S,V] = svds(DD-mDD, floor(sum(selpix)/10));
                Vcorr = V;
                for vix = 1:size(V,2)
                    [~, Vcorr(:, vix)] = obj.decorrelateMotion(V(:,vix), motPreds(:,selframes));
                end
                DDcorr = DD-(U*S*Vcorr'); %motion-corrected data: the original data minus the reconstructed motion artefact
                trace1= mean(DDcorr(Dmask{rix}(selpix), :),1); % the trace computed from the mask alone

                %perform NMF
                nfacs = max(2, floor(sum(selpix)/12));
                [W,H] = nnmf(double(DDcorr),nfacs);

                %find NMF components spatially correlated to trace0, temporally correlated, and with positive skewness
                HP = filtfilt(obj.b2,obj.a2,H')';
                t0HP = filtfilt(obj.b2,obj.a2,trace0);
                Cspace = corr(W, Dmask{rix}(selpix));
                Ctime = corr(HP', t0HP'); %high frequency correlation
                SK = skewness(HP, 0, 2);
                sel = Cspace>0.1 & Ctime>0.1 & SK>0.2; 
                trace2 = mean(W(:,sel)*H(sel,:),1);

%                 %denoise with motion variables %this is deprecated; we
%                 now do SVD-based motion correction before NMF
%                 if ~isempty(obj.aData)
%                     trace3 = obj.decorrelateMotion(trace2, motPreds(:,selframes));
%                     trace0 = obj.decorrelateMotion(trace0, motPreds(:,selframes));
%                     sData.motionSubtraced = true;
%                 else
%                     sData.motionSubtraced = false;
%                 end

                %compile output
                sData.trace0(selframes, rix) = trace0; %raw ROI, without decorrelating motion variables
                sData.n0(rix) = estimatenoise(trace0); %noise for trace0
                sData.trace1(selframes, rix) = trace1; %raw ROI, with motion variables decorrelated
                sData.n1(rix) = estimatenoise(trace1); %noise for trace1
                sData.trace2(selframes, rix) = trace2; %NMF-based traces
                sData.n2(rix) = estimatenoise(trace2); %noise for trace2
                
                sData.bleach = obj.bleach; %bleaching curve for this field of view
                
                sData.contamination(rix) = 1-corr(trace0', trace2'); %a measure of how contaminated the raw trace is by other factors/noise
            end

            %plot the data
            T = obj.aData.frametime * (1:size(sData.trace0,1)); %time vector
            figure, 
            for pix = 1:nROIs
                hAx(pix) = subplot(nROIs,1, pix);
                tNorm = sData.trace2(:,pix)./std(sData.trace2(:,pix)); %normalized trace
                plot(T, tNorm, 'color', [0.5 0.5 0.5]);
                hold on,
                plot(T, smooth(tNorm, 8), 'color', 'k');
                ylabel(sData.names{pix});
                
                if pix<nROIs
                    set(hAx(pix), 'xticklabel', [], 'tickdir', 'out');
                end
            end
            xlabel('time (s)');
            linkaxes(hAx, 'x');

            %save output
            if isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_TRACES.mat']);
        end

        function saveROIs (obj, arg1, arg2)
            if isempty(obj.drsave)
                obj.drsave = obj.dr;
            end
            [obj.fnsave, obj.drsave] = uiputfile([obj.drsave filesep obj.fnStem '_ROIs.mat']);
            if isempty(obj.fnsave)
                return
            end
            pos = {}; labels = {};
            for rix = 1:length(obj.hROIs)
                try
                    posTmp = get(obj.hROIs(rix), 'position');
                    if ~isempty(posTmp)
                        %convert the ROI to saveable data
                        pos = [pos {posTmp}];
                        labels = [labels {get(obj.hROIs(rix), 'Label')}];
                    end
                catch
                    %this object didn't exist
                end
            end
            save([obj.drsave filesep obj.fnsave], 'pos', 'labels');
        end

        function [corrected, traceMot] = decorrelateMotion(obj, trace, motPreds)

            %highpass filter trace
            trace = trace(:) - smooth(trace, ceil(8/obj.aData.frametime), 'lowess');
            beta = mvregress(motPreds',trace);
            traceMot = motPreds'*beta;
            corrected = trace' - traceMot'; 
        end

        function loadROIs (obj, arg1, arg2)
            [fn dr] = uigetfile([obj.drsave filesep '*_ROIS.mat']);

            load([dr filesep fn], 'pos', 'labels');
            delete(obj.hROIs);
            obj.hROIs = [];
            for rix = 1:length(pos)
                obj.hROIs(rix) = images.roi.Polygon(obj.hAx,'position', pos{rix}, 'label', labels{rix});
            end
        end

        function drawPolygon(obj, arg1, arg2)
            L = length(obj.hROIs)+1;
            obj.hROIs(L) = drawpolygon(obj.hAx);
            set(obj.hROIs(L), 'Label', int2str(L));
            cm = get(obj.hROIs(end), 'ContextMenu');
            uimenu(cm,'Text','Set Label','MenuSelectedFcn',@(arg1,arg2)(obj.setLabel(obj.hROIs(end)))); %add a 'Set Label' context menu
        end

        function setLabel(obj, hROI)
                label = inputdlg('Set a label');
                set(hROI, 'Label', label{1});
        end
        function kpf(obj, arg1, arg2)
            switch arg2.Key
                case 'a'
                    @obj.drawPolygon;
                case 'period'
                    obj.Clevel = obj.Clevel/1.5;
                    set(obj.hIm, 'Cdata', cat(3, obj.IMc./obj.Clevel, obj.IMavg, obj.IMavg));
                case 'comma'
                    obj.Clevel = obj.Clevel.*1.5;
                    set(obj.hIm, 'Cdata', cat(3, obj.IMc./obj.Clevel, obj.IMavg, obj.IMavg));

            end
        end
    end
end