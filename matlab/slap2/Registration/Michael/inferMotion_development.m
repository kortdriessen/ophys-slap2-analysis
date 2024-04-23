function inferMotion(varargin)

if numel(varargin) < 1 || isempty(varargin{1})
    [fns_data, dr_data] = uigetfile({'*.meta;*.mat'}, 'Select Data File', 'multiselect', 'on');
    dataFile = [dr_data fns_data];
else
    dataFile = varargin{1};
end

if dataFile(end-3:end) == "meta"
    simulation = false;
elseif dataFile(end-3:end) == ".mat"
    simulation = true;
else
    error('Data must be .meta or .mat');
end

[filepath, ~, ~] = fileparts(dataFile);
filepath = convertStringsToChars(filepath);

if numel(varargin) < 2 || isempty(varargin{2})
    refStackFlag = false;
else
    refStackFlag = true;
    refStackDir = varargin{2};
    refStack = tiffreadVolume(refStackDir);
    if strcmpi(refStackDir(end-6:end-5),'CH')
        channel = str2num(refStackDir(end-4));
    else
        channel = 1;
    end
end 

if numel(varargin) < 3
    ds = 0;
else
    ds = varargin{3};
end

if numel(varargin) < 4
    robust = false;
    b = 3;
else
    robust = varargin{4};
    b = varargin{5};
end

%%
if ~refStackFlag
    [fns, dr] = uigetfile('*.tif', 'Select Reference Stack', 'multiselect', 'on');
    refStack = tiffreadVolume([dr fns]);
    
    if strcmpi(fns(end-6:end-5),'CH')
        channel = str2num(fns(end-4));
    else
        channel = 1;
    end
end

%%
% load('Z:\ophys\SLAP2\exp data\Beads\slap2_000000_2023-06-15_14-07-03\finalFilter.mat');
% refStackDil15 = convn(refStack,finalFilter,'same');

%%
if ~simulation
    hDataViewer = slap2.Slap2DataViewer(dataFile);
    hLowLevelDataFile = hDataViewer.hDataFile.hDataFile;
    fname = hLowLevelDataFile.filename;
    
    numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
    totalCycles = hLowLevelDataFile.numCycles;
    numCycles = floor(totalCycles / 2^ds);
    
    dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
    dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
    numFastZs = length(hLowLevelDataFile.fastZs);
end

%% get list of superpixels and extract data

if ~simulation
    allSuperPixelIDs = [];
    
    for lineSweepIdx = 1:numLinesPerCycle
        superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};
    
        if size(superPixIdxs,1) == 0; continue; end
    
        zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
        
        spIDs = superPixIdxs*100+zIdx;
        allSuperPixelIDs = [allSuperPixelIDs; spIDs];
    end
    
    [allSuperPixelIDs, ~, ic] = unique(allSuperPixelIDs);
    spSampleCt = accumarray(ic, 1);
    
    fprintf("%d superpixels detected\n", length(allSuperPixelIDs));
else
    load(dataFile);
    totalCycles = size(spData,2);
    numCycles = floor(totalCycles / 2^ds);
end

%% get roiMasks

fprintf("Calculating ROI Masks... ");
tic;

if ~simulation
    % using sparse matrix
    sparseMaskInds = [];
    allPixelReplacementMaps = hLowLevelDataFile.metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
    
    for i = 1:length(allSuperPixelIDs)
        tmpMask = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numFastZs);
    
        sp = allSuperPixelIDs(i);
        zIdx = rem(sp,100);
        superPixIdx = (sp - zIdx) / 100;
        pixelReplacementMap = allPixelReplacementMaps{zIdx};
    
        open = double(pixelReplacementMap(pixelReplacementMap(:,2) == superPixIdx,1))+1;
        openR = floor((open-1) /dmdPixelsPerRow)+1;
        openC = open - (openR-1) * dmdPixelsPerRow;
        openPixs = uint32(openR + (openC-1) * dmdPixelsPerColumn + double(zIdx-1) * dmdPixelsPerColumn * dmdPixelsPerRow);
    
        sparseMaskInds = [sparseMaskInds; openPixs, ones(size(openPixs))*i];
    end
    clear('tmpMask');
end

roiMasks = sparse(sparseMaskInds(:,1),sparseMaskInds(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs));
fprintf('done - took %f sec\n', toc);

%% make sure refStack is larger than imaged frame

bl = single(refStack)/100; %  reference image (X x Y x Z)
bl(bl < 0) = 0;
bl = padarray(bl,[floor((max(size(bl,1),dmdPixelsPerColumn)-size(bl,1))/2),...
                  floor((max(size(bl,2),dmdPixelsPerRow)-size(bl,2))/2),...
                  floor((max(size(bl,3),numFastZs)-size(bl,3))/2)],...
              mean(bl(:)),...
              'both');

if size(bl,1) < dmdPixelsPerColumn
    bl = padarray(bl,[1,0,0],mean(bl(:)),'pre');
end

if size(bl,2) < dmdPixelsPerRow
    bl = padarray(bl,[0,1,0],mean(bl(:)),'pre');
end

% bl15 = single(refStackDil15)/100; %  reference image (X x Y x Z)
% bl15 = padarray(bl15,[floor((max(size(bl15,1),dmdPixelsPerColumn)-size(bl15,1))/2),...
%                   floor((max(size(bl15,2),dmdPixelsPerRow)-size(bl15,2))/2),...
%                   floor((max(size(bl15,3),numFastZs)-size(bl15,3))/2)],...
%               mean(bl15(:)),...
%               'both');
% 
% if size(bl15,1) < dmdPixelsPerColumn
%     bl15 = padarray(bl15,[1,0,0],mean(bl15(:)),'pre');
% end
% 
% if size(bl15,2) < dmdPixelsPerRow
%     bl15 = padarray(bl15,[0,1,0],mean(bl15(:)),'pre');
% end

%% Calculate lookup table

xPre = 25; xPost = 25;
yPre = 25; yPost = 25;

zPre = 0;
zPost = 0;
blCropZ = bl(:,:,1+zPre:end-zPost);
% bl15CropZ = bl15(:,:,1+zPre:end-zPost);

xMotRange = xPre + xPost + 1;
yMotRange = yPre + yPost + 1;
zMotRange = size(blCropZ,3) - numFastZs + 1;

if ~simulation
    [lookupFile, lookupDir] = uigetfile('*.mat', 'Select Lookup Table');
end

if simulation || (sum(lookupFile == 0) && simulation) || (sum(lookupFile == 0) && ~exist([fname(1:end-5) '_LOOKUPTABLE.mat'],'file'))
    paddedBL = padarray(padarray(blCropZ,[yPre,xPre,0],mean(blCropZ(:)),'pre'),[yPost,xPost,0],'post');
    paddedSz = size(paddedBL);

%     paddedBL15 = padarray(padarray(bl15CropZ,[yPre,xPre,0],mean(bl15CropZ(:)),'pre'),[yPost,xPost,0],'post');
%     paddedSz15 = size(paddedBL15);

    fprintf("Calculating lookup table... ")
    
    tic;

%     longSPs = [];
%     shortSPs = [];

    likelihood_means = single(zeros(yMotRange, xMotRange, zMotRange,length(allSuperPixelIDs))); % X x Y x Z x no. of superpixels
    for n = 1:length(allSuperPixelIDs)
        openPixs = sparseMaskInds(sparseMaskInds(:,2) == n, 1);

%         shortSPs = [shortSPs n];

%         if numel(openPixs) > 9
%             longSPs = [longSPs n];
%         else
%             shortSPs = [shortSPs n];
%         end

%         nActivePixs = numel(openPixs) + 4;
%         extraPixs = floor((nActivePixs - refStackDil)/2);

        refPix = ceil(length(openPixs)/2);

        [r, c, d] = ind2sub([dmdPixelsPerColumn dmdPixelsPerRow numFastZs], openPixs);
        for y = 1:yMotRange
            for x = 1:xMotRange
                for z = 1:zMotRange
%                     likelihood_means(y,x,z,n) = sum(paddedBL((d+z-2)*paddedSz(1)*paddedSz(2)+(c+x-2)*paddedSz(1)+(r+y-1)));
%                     if numel(openPixs) > 9
%                         likelihood_means(y,x,z,n) = bl15(r(refPix)+y-yPre-1, c(refPix)+x-xPre-1, d(refPix)+z-1);
%                     else
                        likelihood_means(y,x,z,n) = bl(r(refPix)+y-yPre-1, c(refPix)+x-xPre-1, d(refPix)+z-1);
%                     end
                end
            end
        end
    end
    clear('tmpMask');
    
    fprintf('done - took %f sec\n', toc);
    
    if exist('fname')
        save([fname(1:end-5) '_LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
    else
        save([filepath '\LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
    end
elseif sum(lookupFile == 0)
    fprintf("Loading lookup table... ")
    load([fname(1:end-5) '_LOOKUPTABLE.mat']);
    fprintf('done\n');
else
    fprintf("Loading lookup table... ")
    load([lookupDir lookupFile]);
    fprintf('done\n');
end

%%  Calculate log likelihoods and infer motion from MAP
likelihood_means(likelihood_means == 0) = 1e-20;

log_means = log(likelihood_means);
sqrt_means = permute(sqrt(likelihood_means),[4 1 2 3]);

if simulation; numCycles = floor(size(spData,2) / 2^ds); end;

motion = zeros(numCycles,3);
brightness = zeros(numCycles,1);

Mx = -1;
My = -1;
Mz = -1;
searchRadius = 4;

dataMatrix = zeros(length(allSuperPixelIDs),numCycles);
expectedMatrix = zeros(length(allSuperPixelIDs),numCycles);

fprintf("Inferring motion... ")
tic

% residualMovie = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numCycles,'single');

for cycleIdx = 1:numCycles
    scaled_likelihood_means = [];
    scaled_log_means = [];
    if simulation
        data = zeros(length(allSuperPixelIDs),1);
        for c = 1:2^ds
            data = data + spData(:,(cycleIdx-1)*2^ds+c) ./ 100;
        end
    else
        data = zeros(length(allSuperPixelIDs),1);
        for c = 1:2^ds
            allLineData = hLowLevelDataFile.getLineData(1:numLinesPerCycle, ((cycleIdx-1)*2^ds+c)*ones(numLinesPerCycle,1));
            for lineSweepIdx = 1:numLinesPerCycle
        
                superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};
        
                if size(superPixIdxs,1) == 0; continue; end
        
                lineData = allLineData{lineSweepIdx};
                zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
        
                spID = superPixIdxs*100 + uint32(zIdx);
                [~,spIdx] = ismember(spID,allSuperPixelIDs);
        
                data(spIdx(spIdx>0)) = data(spIdx(spIdx>0)) + single(lineData(spIdx>0,channel));
            end
        end
        data = data ./ spSampleCt ./ 100;
    end

    if cycleIdx == 1 
        motionLimits = zeros(xMotRange,yMotRange,zMotRange);
        ySearch = 1:yMotRange;
        xSearch = 1:xMotRange;
        zSearch = 1:zMotRange;
    else
        ySearch = max(1,motion(cycleIdx-1,1) - searchRadius):min(xMotRange,motion(cycleIdx-1,1) + searchRadius);
        xSearch = max(1,motion(cycleIdx-1,2) - searchRadius):min(yMotRange,motion(cycleIdx-1,2) + searchRadius);
        zSearch = max(1,motion(cycleIdx-1,3) - searchRadius):min(zMotRange,motion(cycleIdx-1,3) + searchRadius);
        [Xs, Ys, Zs] = meshgrid(xSearch - motion(cycleIdx-1,2),ySearch - motion(cycleIdx-1,1),zSearch - motion(cycleIdx-1,3));
        motionLimits = zeros(size(Xs)); % ~((Xs) .^ 2 + (Ys) .^ 2 + (Zs) .^ 2 <= searchRadius .^ 2);
    end

    sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,:);
    sub_sqrt_means = sqrt_means(:,ySearch,xSearch,zSearch);
    sub_log_means = log_means(ySearch,xSearch,zSearch,:);

    nonzeroData = data(:);
    nonzero_sq_means = sub_sqrt_means(:,:,:,:);

    scalingFactor = ones(size(sub_likelihood_means,1:3));
%     scalingFactorLong = ones(size(sub_likelihood_means,1:3));
%     scalingFactorShort = ones(size(sub_likelihood_means,1:3));
    for y = 1:length(ySearch)
        for x = 1:length(xSearch)
            for z = 1:length(zSearch)
                if motionLimits(y,x,z)  == 0
%                     scalingFactor(y,x,z) = (nonzero_sq_means(:,y,x,z) \ sqrt(nonzeroData)) .^ 2;
                    scalingFactor(y,x,z) = sum(nonzeroData) / sum(sub_likelihood_means(y,x,z,:));

%                     scalingFactorLong(y,x,z) = (nonzero_sq_means(longSPs,y,x,z) \ sqrt(nonzeroData(longSPs))) .^ 2;
%                     scalingFactorShort(y,x,z) = (nonzero_sq_means(shortSPs,y,x,z) \ sqrt(nonzeroData(shortSPs))) .^ 2;
                    
%                     scalingFactorLong(y,x,z) = sum(nonzeroData(longSPs)) / sum(sub_likelihood_means(y,x,z,longSPs));
%                     scalingFactorShort(y,x,z) = sum(nonzeroData(shortSPs)) / sum(sub_likelihood_means(y,x,z,shortSPs));
% 
%                     if cycleIdx > 1
%                         scalingFactorLong(y,x,z) = max(min(scalingFactorLong(y,x,z),brightnessLong(cycleIdx-1)*1.001),brightnessLong(cycleIdx-1)*0.995);
%                         scalingFactorShort(y,x,z) = max(min(scalingFactorShort(y,x,z),brightnessShort(cycleIdx-1)*1.001),brightnessShort(cycleIdx-1)*0.995);
%                     end
                end
            end
        end
    end
    scaled_likelihood_means = zeros(size(sub_likelihood_means));
    scaled_likelihood_means = scalingFactor .* sub_likelihood_means;
%     scaled_likelihood_means(:,:,:,longSPs) = scalingFactorLong .* sub_likelihood_means(:,:,:,longSPs);
%     scaled_likelihood_means(:,:,:,shortSPs) = scalingFactorShort .* sub_likelihood_means(:,:,:,shortSPs);

    scaled_log_means = zeros(size(sub_log_means));
    scaled_log_means = log(scalingFactor) + sub_log_means;
%     scaled_log_means(:,:,:,longSPs) = log(scalingFactorLong) + sub_log_means(:,:,:,longSPs);
%     scaled_log_means(:,:,:,shortSPs) = log(scalingFactorShort) + sub_log_means(:,:,:,shortSPs);

    if robust
        log_likelihood = sum(max(-b, min(b,(reshape(data,[1,1,1,length(data)]) - scaled_likelihood_means) ./ sqrt(scaled_likelihood_means))).^2,4);
    else
        log_likelihood = sum(reshape(data,[1,1,1,length(data)]) .* scaled_log_means - scaled_likelihood_means,4);
    end

    log_likelihood(log_likelihood >= Inf) = -1e10;

    mat = log_likelihood - motionLimits .* 1e10;

    if robust
        [M, I] = min(abs(mat(:)));
    else
        [M, I] = max(mat(:));
    end

    [My, Mx, Mz, scale] = ind2sub(size(mat),I);
    motion(cycleIdx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)];
    brightness(cycleIdx) = scalingFactor(My, Mx, Mz);
%     brightnessShort(cycleIdx) = scalingFactorShort(My, Mx, Mz);
%     brightnessLong(cycleIdx) = scalingFactorLong(My,Mx,Mz);

    dataMatrix(:,cycleIdx) = data;
    expectedMatrix(:,cycleIdx) = scaled_likelihood_means(My,Mx,Mz,:);

%     residualMovie(:,:,cycleIdx) = single(reshape(roiMasks * double(data - expectedMatrix(:,cycleIdx)) * 100,[dmdPixelsPerColumn dmdPixelsPerRow]));
end

disp(['done - took ' num2str(toc/numCycles) ' sec per frame'])

%% Upsample if downsampled
% totalCycles = 183200;
frames = 1:totalCycles;
dsFrames = frames(1:2^ds:(2^ds*numCycles));

dsMotion = motion;
motion = interp1(dsFrames,dsMotion,frames);

dsBrightness = brightness;
brightness = interp1(dsFrames,dsBrightness,frames,'linear','extrap');

% dsBrightnessShort = brightnessShort;
% brightnessShort = interp1(dsFrames,dsBrightnessShort,frames,'linear','extrap');

% dsBrightnessLong = brightnessLong;
% brightnessLong = interp1(dsFrames,dsBrightnessLong,frames,'linear','extrap');

%% Plot results

figure(401); clf;
ax1 = subplot(3,1,1);
plot(motion(:,1)); hold on;
if simulation; plot((motFinal(:,1)-mean(motFinal(:,1)))+mean(motion(:,1)),'-.'); end; hold off;
title('y motion')
legend('inferred','true','Location','best')

ax2 = subplot(3,1,2);
plot(motion(:,2)); hold on;
if simulation; plot((motFinal(:,2)-mean(motFinal(:,2)))+mean(motion(:,2)),'-.'); end; hold off;
title('x motion')
legend('inferred','true','Location','best')

ax3 = subplot(3,1,3);
plot(motion(:,3)); hold on;
if simulation; plot(motFinal(:,3)-mean(motFinal(:,3))+mean(motion(:,3)),'-.'); end; hold off;
title('z motion')
legend('inferred','true','Location','best')

linkaxes([ax1,ax2])
linkaxes([ax1,ax2,ax3],'x')

saveas(gcf,[filepath '\inferredMotion.fig']);

%% Save out data

inferMotionOut.motion = motion - [xPre+1 yPre+1 0];
inferMotionOut.brightness = brightness';
% inferMotionOut.brightnessShort = brightnessShort';
% inferMotionOut.brightnessLong = brightnessLong';
inferMotionOut.dataMatrix = dataMatrix;
inferMotionOut.expectedMatrix = expectedMatrix;
inferMotionOut.sparseMaskInds = sparseMaskInds;

if ~simulation
    if robust
        save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_ROBUST_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
    else
        save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
    end
else
    if robust
        save([filepath '\INFER_MOTION_OUT_ROBUST.mat'],'inferMotionOut','-v7.3');
    else
        save([filepath '\INFER_MOTION_OUT.mat'],'inferMotionOut','-v7.3');
    end
end