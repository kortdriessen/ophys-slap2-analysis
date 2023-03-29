simulation = false;
robust = false;
b = 3;
maxSuperPixels = 2000;

refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_666161_2023-03-28_16-07-46\zstack_20230328_163403_DMD1-REFERENCE_CH1.tif");

% refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_664321_2023-03-22_10-02-10\refStack_20230322_124704_DMD1-REFERENCE_CH2.tif");

% refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_664321_2023-03-22_10-02-10\refStack_20230322_124704_DMD1-REFERENCE_CH2.tif");
% refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_664321_2023-03-22_10-02-10\refStack_20230322_124704_DMD2-REFERENCE_CH2.tif");
% refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_664321_2023-03-15_11-56-43\refStack1_20230315_130400_DMD1-REFERENCE_CH2.tif");
% refStack = tiffreadVolume("Z:\ophys\SLAP2\exp data\Mice\slap2_664321_2023-03-15_11-56-43\refStack2_20230315_130853_DMD1-REFERENCE_CH2.tif");
size(refStack)

hDataViewer = slap2.Slap2DataViewer();

%%
hLowLevelDataFile = hDataViewer.hDataFile.hDataFile;

numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
numCycles = hLowLevelDataFile.numCycles;

dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;

%% get list of superpixels and extract data

allSuperPixelIDs = [];

for lineSweepIdx = 1:numLinesPerCycle
    superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};

    if size(superPixIdxs,1) == 0; continue; end

    zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
    
    spIDs = superPixIdxs*100+zIdx;
    allSuperPixelIDs = [allSuperPixelIDs; spIDs];

%     if mod(lineSweepIdx,10000) == 0; disp(length(allSuperPixelIDs)); allSuperPixelIDs = unique(allSuperPixelIDs); end;
end

[allSuperPixelIDs, ~, ic] = unique(allSuperPixelIDs);
spSampleCt = accumarray(ic, 1);

fprintf("%d superpixels detected\n", length(allSuperPixelIDs));


%% get roiMasks

spToUse = 1:length(allSuperPixelIDs);

if length(allSuperPixelIDs) > maxSuperPixels
    roiMasks = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,length(hLowLevelDataFile.fastZs),1);
    
    fprintf("Too many superpixels - limiting to %d\n", maxSuperPixels);
    s = RandStream('mlfg6331_64'); 
    cellLocation = [381 755];
    captureRadius = 25;
    ct = 1;
    for i = 1:length(allSuperPixelIDs)
        if ct > (2*captureRadius+1)^2; spToUse(spToUse >= i) = []; break; end
        sp = allSuperPixelIDs(i);
        zIdx = rem(sp,100);
        superPixIdx = (sp - zIdx) / 100;
    
        pixelReplacementMap = hLowLevelDataFile.metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
        pixelReplacementMap = pixelReplacementMap{zIdx};
    
        tmp = zeros(dmdPixelsPerRow, dmdPixelsPerColumn*2);
        tmp(superPixIdx+1) = 1;
        tmp(pixelReplacementMap(:,1)+1) = tmp(pixelReplacementMap(:,2)+1);

        if sum(tmp(cellLocation(2)-captureRadius:cellLocation(2)+captureRadius, cellLocation(1)-captureRadius:cellLocation(1)+captureRadius), 'all') < 1;
            spToUse(spToUse == i) = [];
            continue;
        end

        roiMasks(:,:,zIdx,ct) = tmp(:,1:dmdPixelsPerColumn)';
        ct = ct + 1;
    end
else
    roiMasks = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,length(hLowLevelDataFile.fastZs),length(spToUse));
    
    for i = 1:length(allSuperPixelIDs)
        sp = allSuperPixelIDs(i);
        zIdx = rem(sp,100);
        superPixIdx = (sp - zIdx) / 100;
    
        pixelReplacementMap = hLowLevelDataFile.metaData.AcquisitionContainer.ParsePlan.pixelReplacementMaps;
        pixelReplacementMap = pixelReplacementMap{zIdx};
    
        tmp = zeros(dmdPixelsPerRow, dmdPixelsPerColumn*2);
        tmp(superPixIdx+1) = 1;
        tmp(pixelReplacementMap(:,1)+1) = tmp(pixelReplacementMap(:,2)+1);
    
        roiMasks(:,:,zIdx,i) = tmp(:,1:dmdPixelsPerColumn)';
    end
end

% sparse matrix

%% make sure refStack is larger than imaged frame

bl = single(refStack)/100; %  reference image (X x Y x Z)
bl = padarray(bl,[floor((max(size(bl,1),size(roiMasks,1))-size(bl,1))/2),...
                  floor((max(size(bl,2),size(roiMasks,2))-size(bl,2))/2),...
                  floor((max(size(bl,2),size(roiMasks,2))-size(bl,2))/2)],...
              mean(bl(:)),...
              'both');

if size(bl,1) < size(roiMasks,1)
    bl = padarray(bl,[1,0,0],mean(bl(:)),'pre');
end

if size(bl,2) < size(roiMasks,2)
    bl = padarray(bl,[0,1,0],mean(bl(:)),'pre');
end

% bl = padarray(padarray(bl(1:end-20,201:end,:),[0 200 0],mean(bl(:)),'post'),[20 0 0], mean(bl(:)),'pre');

%% Calculate lookup table

xPre = 40; xPost = 40;
yPre = 40; yPost = 40;

xMotRange = xPre + xPost + 1;
yMotRange = yPre + yPost + 1;
zMotRange = size(bl,3) - size(roiMasks,3) + 1;

fname = hLowLevelDataFile.filename;
if ~exist([fname(1:end-5) '_LOOKUPTABLE.mat'],'file') % length(allSuperPixelIDs) > maxSuperPixels ||
    fprintf("Calculating lookup table... ")
    
    tic;

    likelihood_means = single(zeros(yMotRange, xMotRange, zMotRange,length(spToUse))); % X x Y x Z x no. of superpixels
    parfor n = 1:length(spToUse)
        likelihood_means(:,:,:,n) = single(convn(padarray(padarray(bl,[yPre,xPre,0],mean(bl(:)),'pre'),[yPost,xPost,0],'post'),roiMasks(end:-1:1,end:-1:1,end:-1:1,n),'valid'));
    end
    
    fprintf('done - took %f sec\n', toc);
    
    save([fname(1:end-5) '_LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
else
    fprintf("Loading lookup table... ")
    load([fname(1:end-5) '_LOOKUPTABLE.mat']);
    fprintf('done\n');
end

%%  Calculate log likelihoods and infer motion from MAP

log_means = log(likelihood_means);
sqrt_means = permute(sqrt(likelihood_means),[4 1 2 3]);

if simulation; numCycles = size(spData,2); end;

motion = zeros(numCycles,3);
brightness = zeros(numCycles,1);

Mx = -1;
My = -1;
Mz = -1;
searchRadius = 2;

dataMatrix = zeros(length(spToUse),numCycles);
expectedMatrix = zeros(length(spToUse),numCycles);

fprintf("Inferring motion... ")
tic

for cycleIdx = 1:numCycles
    if simulation; data = spData(:,cycleIdx) ./ 100;
    else
        data = zeros(length(spToUse),1);
        allLineData = hLowLevelDataFile.getLineData(1:numLinesPerCycle, cycleIdx);
        for lineSweepIdx = 1:numLinesPerCycle
    
            superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};
    
            if size(superPixIdxs,1) == 0; continue; end
    
            lineData = allLineData{lineSweepIdx};
            zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
    
            spID = superPixIdxs*100 + uint32(zIdx);
            [~,spIdx] = ismember(spID,allSuperPixelIDs(spToUse));
    
            data(spIdx(spIdx>0)) = data(spIdx(spIdx>0)) + single(lineData(spIdx>0));
        end
        data = data ./ spSampleCt(spToUse) ./ 100;
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
        [Ys, Xs, Zs] = meshgrid(ySearch - motion(cycleIdx-1,1),xSearch - motion(cycleIdx-1,2),zSearch - motion(cycleIdx-1,3));
        motionLimits = ~((Xs) .^ 2 + (Ys) .^ 2 + (Zs) .^ 2 <= searchRadius .^ 2);
    end

    sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,:);
    sub_sqrt_means = sqrt_means(:,ySearch,xSearch,zSearch);
    sub_log_means = log_means(ySearch,xSearch,zSearch,:);

    
%     scalingFactor = sqrt(reshape(likelihood_means,[xMotRange*yMotRange*zMotRange length(data)]))' \ sqrt(data);

    scalingFactor = ones(size(sub_likelihood_means,1:3));
    for y = 1:length(ySearch)
        for x = 1:length(xSearch)
            for z = 1:length(zSearch)
                if motionLimits(y,x,z)  == 0
                    scalingFactor(y,x,z) = (sub_sqrt_means(:,y,x,z) \ sqrt(data)) .^ 2;
                end
            end
        end
    end
    scaled_likelihood_means = scalingFactor .* sub_likelihood_means;
    scaled_log_means = log(scalingFactor) + sub_log_means;
    
%     tmp_lh_means = likelihood_means;
%     tmp_lh_means(tmp_lh_means < 0.25) = 0;
% 
%     scaled_likelihood_means = likelihood_means ./ sum(tmp_lh_means,4) .* sum(data);
%     scaled_likelihood_means(likelihood_means < 0.25) = likelihood_means(likelihood_means < 0.25);
%     scaled_log_means = log(scaled_likelihood_means);

    if robust
        log_likelihood = sum(abs(max(-b, min(b,(reshape(data,[1,1,1,length(data)]) - likelihood_means) ./ sqrt(likelihood_means)))),4);
    else
        log_likelihood = sum(bsxfun(@times,reshape(data,[1,1,1,length(data)]),scaled_log_means) - scaled_likelihood_means,4);
    end

    log_likelihood(log_likelihood >= Inf) = -1e8;

    mat = log_likelihood - motionLimits .* 1e8;

    if robust
        [M, I] = min(abs(mat(:)));
    else
        [M, I] = max(mat(:));
    end

    [My, Mx, Mz, scale] = ind2sub(size(mat),I);
    motion(cycleIdx,:) = [ySearch(My); xSearch(Mx); zSearch(Mz)];
    brightness(cycleIdx) = scalingFactor(My, Mx, Mz);
%     if cycleIdx > 1 && abs(motion(cycleIdx,1) - motion(cycleIdx - 1,1)) > r
%         disp('debug');
%     end
%     disp(motion(cycleIdx,:));

    dataMatrix(:,cycleIdx) = data;
    expectedMatrix(:,cycleIdx) = scaled_likelihood_means(My,Mx,Mz,:);
end

disp(['done - took ' num2str(toc/numCycles) ' sec per frame'])

%% Plot results

figure(401); clf;
subplot(4,1,1);
plot(motion(:,1)); hold on;
if simulation; plot(mot_hist(:,1)-mot_hist(1,1),'-.'); end; hold off;
title('y motion')
legend('inferred','true','Location','southeast')

subplot(4,1,2);
plot(motion(:,2)); hold on;
if simulation; plot(mot_hist(:,2)-mot_hist(1,2),'-.'); end; hold off;
title('x motion')
legend('inferred','true')

subplot(4,1,3);
plot(motion(:,3)); hold on;
if simulation; plot(mot_hist(:,3)-mot_hist(1,3),'-.'); end; hold off;
title('z motion')
legend('inferred','true')