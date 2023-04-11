simulation = false;
robust = false;
b = 3;
ds = 0;

%%
[fns, dr] = uigetfile('*.tif', 'Select Reference Stack', 'multiselect', 'on');
refStack = tiffreadVolume([dr fns]);

size(refStack)

%%
if ~simulation
    hDataViewer = slap2.Slap2DataViewer();
    hLowLevelDataFile = hDataViewer.hDataFile.hDataFile;
    fname = hLowLevelDataFile.filename;
    
    numLinesPerCycle = hLowLevelDataFile.header.linesPerCycle;
    totalCycles = hLowLevelDataFile.numCycles;
    numCycles = totalCycles / 2^ds;
    
    dmdPixelsPerColumn = hLowLevelDataFile.metaData.dmdPixelsPerColumn;
    dmdPixelsPerRow = hLowLevelDataFile.metaData.dmdPixelsPerRow;
    numFastZs = length(hLowLevelDataFile.fastZs);
end

%% get list of superpixels and extract data

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


%% get roiMasks

% using sparse matrix
fprintf("Calculating ROI Masks... ");
tic;
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
roiMasks = sparse(sparseMaskInds(:,1),sparseMaskInds(:,2),1,dmdPixelsPerColumn * dmdPixelsPerRow * numFastZs,length(allSuperPixelIDs));
fprintf('done - took %f sec\n', toc);

%% make sure refStack is larger than imaged frame

bl = single(refStack)/100; %  reference image (X x Y x Z)
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

%% Calculate lookup table

xPre = 25; xPost = 25;
yPre = 25; yPost = 25;

zPre = 0;
zPost = 0;
blCropZ = bl(:,:,1+zPre:end-zPost);

xMotRange = xPre + xPost + 1;
yMotRange = yPre + yPost + 1;
zMotRange = size(blCropZ,3) - numFastZs + 1;

[lookupFile, lookupDir] = uigetfile('*.mat', 'Select Lookup Table');

if sum(lookupFile == 0) && ~exist([fname(1:end-5) '_LOOKUPTABLE.mat'],'file')
    paddedBL = padarray(padarray(blCropZ,[yPre,xPre,0],mean(blCropZ(:)),'pre'),[yPost,xPost,0],'post');
    paddedSz = size(paddedBL);

    fprintf("Calculating lookup table... ")
    
    tic;

    likelihood_means = single(zeros(yMotRange, xMotRange, zMotRange,length(allSuperPixelIDs))); % X x Y x Z x no. of superpixels
    parfor n = 1:length(allSuperPixelIDs)
        openPixs = sparseMaskInds(sparseMaskInds(:,2) == n, 1);
        [r, c, d] = ind2sub([dmdPixelsPerColumn dmdPixelsPerRow numFastZs], openPixs);
        for y = 1:yMotRange
            for x = 1:xMotRange
                for z = 1:zMotRange
                    likelihood_means(y,x,z,n) = sum(paddedBL((d+z-2)*paddedSz(1)*paddedSz(2)+(c+x-2)*paddedSz(1)+(r+y-1)));
                end
            end
        end
    end
    clear('tmpMask');
    
    fprintf('done - took %f sec\n', toc);
    
    save([fname(1:end-5) '_LOOKUPTABLE.mat'],'likelihood_means','-v7.3');
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

if simulation; numCycles = size(spData,2); end;

motion = zeros(numCycles,3);
brightness = zeros(numCycles,1);

Mx = -1;
My = -1;
Mz = -1;
searchRadius = 2;

dataMatrix = zeros(length(allSuperPixelIDs),numCycles);
expectedMatrix = zeros(length(allSuperPixelIDs),numCycles);

fprintf("Inferring motion... ")
tic

for cycleIdx = 1:numCycles
    if simulation; data = spData(:,cycleIdx) ./ 100;
    else
        data = zeros(length(allSuperPixelIDs),1);
        for c = 1:2^ds
            allLineData = hLowLevelDataFile.getLineData(1:numLinesPerCycle, (cycleIdx-1)*2^ds+c);
            for lineSweepIdx = 1:numLinesPerCycle
        
                superPixIdxs = hLowLevelDataFile.lineSuperPixelIDs{lineSweepIdx};
        
                if size(superPixIdxs,1) == 0; continue; end
        
                lineData = allLineData{lineSweepIdx};
                zIdx = hLowLevelDataFile.lineFastZIdxs(lineSweepIdx);
        
                spID = superPixIdxs*100 + uint32(zIdx);
                [~,spIdx] = ismember(spID,allSuperPixelIDs);
        
                data(spIdx(spIdx>0)) = data(spIdx(spIdx>0)) + single(lineData(spIdx>0));
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
        motionLimits = ~((Xs) .^ 2 + (Ys) .^ 2 + (Zs) .^ 2 <= searchRadius .^ 2);
    end

    sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,:);
    sub_sqrt_means = sqrt_means(:,ySearch,xSearch,zSearch);
    sub_log_means = log_means(ySearch,xSearch,zSearch,:);

    nonzeroData = data(data > 0);
    nonzero_sq_means = sub_sqrt_means(data > 0,:,:,:);

    scalingFactor = ones(size(sub_likelihood_means,1:3));
    for y = 1:length(ySearch)
        for x = 1:length(xSearch)
            for z = 1:length(zSearch)
                if motionLimits(y,x,z)  == 0
                    scalingFactor(y,x,z) = (nonzero_sq_means(:,y,x,z) \ sqrt(nonzeroData)) .^ 2;
                end
            end
        end
    end
    scaled_likelihood_means = scalingFactor .* sub_likelihood_means;
    scaled_log_means = log(scalingFactor) + sub_log_means;

    if robust
        log_likelihood = sum((max(-b, min(b,(reshape(data,[1,1,1,length(data)]) - scaled_likelihood_means) ./ sqrt(scaled_likelihood_means)))).^2,4);
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

    dataMatrix(:,cycleIdx) = data;
    expectedMatrix(:,cycleIdx) = scaled_likelihood_means(My,Mx,Mz,:);
end

disp(['done - took ' num2str(toc/numCycles) ' sec per frame'])

%% Plot results

figure(401); clf;
subplot(3,1,1);
plot(motion(:,1)); hold on;
if simulation; plot(-(motFinal(:,1)-mean(motFinal(:,1)))+mean(motion(:,1)),'-.'); end; hold off;
title('y motion')
legend('inferred','true','Location','best')

subplot(3,1,2);
plot(motion(:,2)); hold on;
if simulation; plot(-(motFinal(:,2)-mean(motFinal(:,2)))+mean(motion(:,2)),'-.'); end; hold off;
title('x motion')
legend('inferred','true','Location','best')

subplot(3,1,3);
plot(motion(:,3)); hold on;
if simulation; plot(motFinal(:,3)-mean(motFinal(:,3))+mean(motion(:,3)),'-.'); end; hold off;
title('z motion')
legend('inferred','true','Location','best')

%% Save out data

inferMotionOut.motion = motion - [xPre+1 yPre+1 0];
inferMotionOut.brightness = brightness;
inferMotionOut.dataMatrix = dataMatrix;
inferMotionOut.expectedMatrix = expectedMatrix;
inferMotionOut.sparseMaskInds = sparseMaskInds;

% save([fname(1:end-5) sprintf('_INFER_MOTION_OUTPUT_DS_%dx.mat',2^ds)],'inferMotionOut','-v7.3');
save('INFER_MOTION_OUT.mat','inferMotionOut','-v7.3');
