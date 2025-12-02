function E = processTrial_vIM(dr, fnRaw, startLine, endLine, selPix, sources, discardFrames, alignData, meanIM, motOutput, roiData, params)
fn = fnRaw;
disp('Loading high-res data for file:')
disp([dr filesep fn])

switch params.microscope
    case 'SLAP2'
        if params.includeIntegrationROIs
            spTypeFlag = 0; %use all superpixel types
        else
            spTypeFlag = 1; %use only raster superpixels
        end

        %load the high time resolution data
        S2data = slap2.Slap2DataFile([dr filesep fn]);
        meta = loadMetadata([dr filesep fn]);
        numChannels = params.numChannels;

        linerateHz = 1/meta.linePeriod_s;
        dt = linerateHz/params.analyzeHz;
        frameLines = ceil(startLine:dt:endLine);
        nFrames= length(frameLines);
        selPx2D = any(selPix,3);
        sz = size(selPx2D);

        %upsample motion
        motionC = interp1(alignData.DSframes, alignData.motionDSc, frameLines, 'pchip', 'extrap') + motOutput(2);
        motionR = interp1(alignData.DSframes, alignData.motionDSr, frameLines, 'pchip', 'extrap') + motOutput(1);

        labeled = medfilt2(meanIM(:,:,params.activityChannel), [3 3]);
        labeled = ~isnan(meanIM(:,:,params.activityChannel)) & labeled>3*prctile(labeled(~isnan(labeled)), 25); %labeled pixels
        for cix = numChannels:-1:1
            tmp = double(meanIM(:,:,cix));
            mLabeled{cix} = tmp(labeled);
        end

        Ysz = [length(alignData.trimRows) length(alignData.trimCols)];
        selPx2D = selPx2D(1:size(alignData.viewC,1), 1:size(alignData.viewC,2));

        orderedChannels = [params.activityChannel:numChannels, 1:params.activityChannel-1];

        for fix = nFrames:-1:1
            for cix = 1:numel(orderedChannels)               
                if cix==1
                    [Y, Fr_tmp] = getImageWrapper(S2data, orderedChannels(cix), frameLines(fix), ceil(dt), 1, spTypeFlag);
                    Y = Y(alignData.trimRows, alignData.trimCols); 
                    Fr_tmp = Fr_tmp(alignData.trimRows, alignData.trimCols);
                    [Y, Fr_tmp2] = interpFrame(Y,alignData.viewC(1,:)+motionC(fix), alignData.viewR(:,1)+motionR(fix), Fr_tmp); %Y = interp2(1:Ysz(2), 1:Ysz(1), Y,alignData.viewC+motionC(fix), alignData.viewR+motionR(fix), 'linear', nan);
                    IMsel(:, fix) = Y(selPx2D);
                    Finvsel(:,fix) = Fr_tmp2(selPx2D);
                elseif cix==2
                    Y = getImageWrapper(S2data, orderedChannels(cix), frameLines(fix), ceil(dt), 1, spTypeFlag);
                    Y = Y(alignData.trimRows, alignData.trimCols);
                    Y = interpFrame(Y,alignData.viewC(1,:)+motionC(fix), alignData.viewR(:,1)+motionR(fix), Fr_tmp);
                    IM2sel(:, fix) = Y(selPx2D);
                end

                %compute global ROI activity
                mIM = meanIM(:,:,orderedChannels(cix));
                yLabeled = double(Y(labeled)); nans= isnan(yLabeled);
                FF = (sum(yLabeled(~nans))./sum(mLabeled{orderedChannels(cix)}(~nans))).*sum(mLabeled{orderedChannels(cix)});
                E.global.F(fix,cix) = FF;

                %compute user ROI activity
                for rix = 1:length(roiData)
                    mask = roiData{rix}.mask;
                    tmp1 = Y(mask); tmp2 = mIM(mask);   tmp2(isnan(tmp2)) = 0;
                    nans= isnan(tmp1);
                    Fpx{rix}(:,fix,cix) = tmp1;
                    FF = (sum(tmp1(~nans))./sum(tmp2(~nans))).*sum(tmp2);
                    E.ROIs.F_fullSpeed(rix, fix,cix) = FF;
                end
            end
        end

        %extract downsampled soma data
        if ~isempty(roiData) %if any of the rois is labeled 'soma'
            if params.roiHz==params.analyzeHz
                E.ROIs.F = E.ROIs.F_fullSpeed;
                somaLines =frameLines;
            else
                Fpx = {};
                dtSoma = linerateHz/params.roiHz;
                somaLines = ceil(startLine:dtSoma:endLine);
                nFramesSoma= length(somaLines);
                roiFtmp = nan(length(roiData),nFramesSoma,numel(orderedChannels));
                E.ROIs.F = nan(length(roiData),nFrames,numel(orderedChannels));
                for fix = nFramesSoma:-1:1
                    for cix = 1:numel(orderedChannels)
                        % Y = S2data.getImage(orderedChannels(cix), somaLines(fix), ceil(dtSoma), 1);
                        [Y, ~] = getImageWrapper(S2data, orderedChannels(cix), somaLines(fix), ceil(dtSoma), 1, spTypeFlag);
                        Y = Y(alignData.trimRows, alignData.trimCols);
                        Y = interp2(1:Ysz(2), 1:Ysz(1), Y,alignData.viewC+motionC(fix), alignData.viewR+motionR(fix), 'linear', nan);
                        mIM = meanIM(:,:,orderedChannels(cix));
                        for rix = 1:length(roiData)
                            mask = roiData{rix}.mask;
                            tmp1 = Y(mask); tmp2 = mIM(mask);   tmp2(isnan(tmp2)) = 0;
                            nans= isnan(tmp1);
                            Fpx{rix}(:,fix,cix) = tmp1;
                            FF = (sum(tmp1(~nans))./sum(tmp2(~nans))).*sum(tmp2);
                            roiFtmp(rix, fix,cix) = FF;
                        end
                    end
                end
                for cix = 1:numel(orderedChannels)
                    for rix = 1:length(roiData)
                        E.ROIs.F(rix,:,cix) = interp1(somaLines, roiFtmp(rix, :,cix), frameLines);
                    end
                end
            end
        end

        %perform SVD on user ROIs to denoise
        E.ROIs.Fsvd = nan(length(roiData), nFrames, numChannels);
        for rix = 1:length(roiData)
            for cix = 1:numel(orderedChannels)
                Dtmp = double(Fpx{rix}(:,:,cix));
                [UU,SS,VV,bg] = nansvd(Dtmp,3, 10, params.nanThresh);
                roiLikeness = (abs(mean(UU,1, 'omitnan'))./sqrt(mean(UU.^2,1, 'omitnan')))*SS;
                [~,selPC] = max(roiLikeness);
                Fsvd = mean(bg+(UU(:,selPC)*SS(selPC,selPC)*VV(:,selPC)'),1, 'omitnan');
                %upsample to desired framerate if necessary
                E.ROIs.Fsvd(rix,:,cix) = interp1(somaLines, Fsvd, frameLines);
            end
        end
        clear Fpx;

    case 'bergamo'
        activityChannel = params.activityChannel;
        numChannels = params.numChannels;
        fn = fnRaw;
        IM = networkScanImageTiffReader([dr filesep fn]);
        IM = double(IM);

        selPx2D = any(selPix,3);

        %rearrange IM into correct dimensions
        IM = reshape(IM, size(IM,1), size(IM,2), numChannels, []);

        if activityChannel>1
            IM = IM(:,:,[activityChannel:end, 1:activityChannel-1],:);
            disp('Reordering channels for analysis!')
        end
        IM1 = squeeze(IM(:,:,1,:));
        IMsel = interpArray(IM1, selPx2D, motOutput); %interpolate the movie at the shifted coordinates

        if numChannels==2
            IM2 =  squeeze(IM(:,:,2,:));
            clear IM;
            IM2sel = interpArray(IM2, selPx2D, motOutput); %interpolate the movie at the shifted coordinates
            clear IM2;
        else %1 channel
            clear IM;
        end

end

discard = interp1(1:numel(discardFrames), double(discardFrames(:)), linspace(1, numel(discardFrames), size(IMsel,2)))>0; %upsample the discard frames
IMsel(:,discard) = nan;     %throw away movement frames as above
E.global.F(discard,:) = nan;
E.ROIs.F_fullSpeed(:,discard,:) = nan;
E.ROIs.F(:, discard,:) = nan;
E.ROIs.Fsvd(:,discard,:) = nan;
if numChannels ==2
    IM2sel(:,discard) = nan;
end


%fill nans in IMsel
IMnans = isnan(IMsel);
fillMatrix = repmat(mean(IMsel,2, 'omitmissing'), 1, size(IMsel,2));
IMsel(IMnans) = fillMatrix(IMnans);
Finvsel(IMnans) = mean(Finvsel,'all','omitmissing').*1000;

[H,B,S,LS,F0,SNR] = extractTrial(IMsel,Finvsel, sources, selPx2D, params);
X = convn(S, params.k, 'same');
X(:,discard) = nan;
S(:,discard) = nan;

%populate E
E.dF.events(:,:,1) = S; %[source#, time, channel]
E.dF.denoised(:,:,1) = X; %[source#, time, channel]
E.dF.ls(:,:,1) = LS; %[source#, time, channel] LEAST SQUARES SOLVE
E.F0(:,:,1) = F0;
E.footprints = single(Wfull);
E.discardFrames = discard;

end

function [IM, freshness] = getImageWrapper(S2data, channel, frames, dt, zPlane, spTypeFlag)
if spTypeFlag
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane, spTypeFlag);
else
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane); %for backward compatibility
end
end


% function [IMsel, IM2sel, FRsel, globalF, roiF, roiPx] = readFrameAsync(S2data, orderedChannels, frameLines, dt, spTypeFlag, alignData, motionR, motionC, selPx2D, roiData, meanIM,labeled, mLabeled)
%             for cix = 1:numel(orderedChannels)               
%                 if cix==1
%                     [Y, Fr_tmp] = getImageWrapper(S2data, orderedChannels(cix), frameLines, dt, 1, spTypeFlag);
%                     Y = Y(alignData.trimRows, alignData.trimCols); 
%                     Fr_tmp = Fr_tmp(alignData.trimRows, alignData.trimCols);
%                     [Y, Fr_tmp2] = interpFrame(Y,alignData.viewC(1,:)+motionC, alignData.viewR(:,1)+motionR, Fr_tmp); %Y = interp2(1:Ysz(2), 1:Ysz(1), Y,alignData.viewC+motionC(fix), alignData.viewR+motionR(fix), 'linear', nan);
%                     IMsel = Y(selPx2D);
%                     FRsel = Fr_tmp2(selPx2D);
%                 elseif cix==2
%                     Y = getImageWrapper(S2data, orderedChannels(cix), frameLines, dt, 1, spTypeFlag);
%                     Y = Y(alignData.trimRows, alignData.trimCols);
%                     Y = interpFrame(Y,alignData.viewC(1,:)+motionC, alignData.viewR(:,1)+motionR, Fr_tmp);
%                     IM2sel = Y(selPx2D);
%                 end
% 
%                 %compute global ROI activity
%                 mIM = meanIM(:,:,orderedChannels(cix));
%                 yLabeled = double(Y(labeled)); nans= isnan(yLabeled);
%                 FF = (sum(yLabeled(~nans))./sum(mLabeled{orderedChannels(cix)}(~nans))).*sum(mLabeled{orderedChannels(cix)});
%                 globalF(cix) = FF;
% 
%                 %compute user ROI activity
%                 roiPx = {}; roiF = [];
%                 for rix = length(roiData):-1:1
%                     mask = roiData{rix}.mask;
%                     tmp1 = Y(mask); tmp2 = mIM(mask);   tmp2(isnan(tmp2)) = 0;
%                     nans= isnan(tmp1);
%                     roiPx{rix}(:,1,cix) = tmp1;
%                     FF = (sum(tmp1(~nans))./sum(tmp2(~nans))).*sum(tmp2);
%                     roiF(rix, 1,cix) = FF;
%                 end
%             end
% end