trialDirs = {'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV1\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV2\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV3\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV4\'};
brightnesses = [0.0625 0.125 0.25 0.5 1 2 4 8];

K = 10;
Nt = 600;

frameRate = 50; % fps
spikeTau = 20; % ms
lambda = 1.5; % in spikes per s


for trial = 1:length(trialDirs)
    fprintf('Trial %d\n',trial);
    dr = trialDirs{trial};
    fns = dir([dr '*_CONFIG2-REFERENCE_CH2.tif']);
    refStack11 = tiffreadVolume([dr fns(1).name]);
    
    fns = dir([dr '*_CONFIG3-REFERENCE_CH2.tif']);
    refStack15 = tiffreadVolume([dr fns(1).name]);

    load('Z:\ophys\SLAP2\exp data\Beads\slap2_000000_2023-06-15_14-07-03\finalFilter.mat');
%     refStack15 = convn(refStack11,finalFilter,'same');

    [dmdPixelsPerColumn, dmdPixelsPerRow, ~] = size(refStack11);
    numFastZs = 1;

    if ~exist([dr '\simulation\simulatedActivity.mat'])
        % select and create ROIs
    
        A11 = [];
        A15 = [];
        phi = [];
    
    
        fig = figure(121); clf;
        imshow3D(refStack11);
        hAx = fig.CurrentAxes;
    
        for k = 1:K
            ellipse = drawellipse(hAx);
            pause();
    
            mask = imgaussfilt(double(createMask(ellipse)));
    
            footprint11 = mask .* refStack11(:,:,9);
            mask15 = convn(mask,finalFilter,'same');
            mask15 = mask15 ./ max(mask15(:));
            footprint15 = mask15 .* refStack15(:,:,9);
    
            A11(:,end+1) = footprint11(:);
            A15(:,end+1) = footprint15(:);
    
            spikes = zeros(Nt*10,1);
        
            t = 1;
            while t <= Nt*10
                isi = round(random('Exponential',1/(lambda/frameRate/10)));
                if t + isi <= Nt*10
                    spikes(t+isi) = 5;
                end
                t = t+isi;
            end
            
            transients = conv(spikes,exp(-(0:500)/(spikeTau * frameRate * 10 / 1000)));
            phi(:,end+1) = transients(1:10:Nt*10);
        end
        save([dr '\simulation\simulatedActivity.mat'],'A11','A15','phi');
    else
        load([dr '\simulation\simulatedActivity.mat']);
    end
    
    load([dr 'simulation\simulationROI.mat']);
    [rasterROI,blockROI,rasterTargetHz,blockTargetHz] = rois.translateRoiShapes([dmdPixelsPerColumn dmdPixelsPerRow]);
    
    roiMasks = zeros([size(blockROI, 1:2) 0]);
    spData = [];
    allSuperPixelIDs = [];
    spCounter = 0;
    longSPs = [];
    shortSPs = [];
    
    refRs = [];
    refCs = [];
    refPixs = [];
    
    for roiIdx = 1:max(blockROI(:))
        for colIdx = find(sum(blockROI == roiIdx))
            if sum(blockROI(:,colIdx) == roiIdx) == 0; 
                continue; 
            end;
    
            spCounter = spCounter + 1;
            allSuperPixelIDs(end+1) = spCounter;
            roiMasks(:,colIdx,spCounter) = uint8(blockROI(:,colIdx) == roiIdx);
    
            openPixs = find(roiMasks(:,:,spCounter));
    
            if numel(openPixs) > 9
                longSPs = [longSPs spCounter];
            else
                shortSPs = [shortSPs spCounter];
            end
    
            [refR, refC] = ind2sub(size(blockROI),median(openPixs));
            refRs = [refRs refR];
            refCs = [refCs refC];
    
            refPixs = [refPixs median(openPixs)];
        end
    end
    
    numSuperPixels = size(roiMasks,3);
    allSuperPixelIDs = 1:spCounter;
    
    [r,c] = find(reshape(roiMasks,[dmdPixelsPerColumn*dmdPixelsPerRow*numFastZs length(allSuperPixelIDs)]));
    sparseMaskInds = [r, c];

    tmp = load([dr 'simulation\simulatedMotion.mat']);
    motFinal = tmp.motFinal; clear tmp;

    motFinal(:,3) = 9 * ones(size(motFinal(:,3))); % override Z motion

    Nt = size(motFinal,1);
    
    for b = 1:length(brightnesses)
        bright = brightnesses(b);
        
        pbleach = bright .* exp(-(1:Nt)/1000);
        
        spDataClean = nan(numSuperPixels,Nt);

        activityMovie11 = reshape(A11 * phi',[dmdPixelsPerColumn dmdPixelsPerRow Nt]);
        activityMovie15 = reshape(A15 * phi',[dmdPixelsPerColumn dmdPixelsPerRow Nt]);
        
        for i = 1:Nt
            refFrame11 = activityMovie11(:,:,i) + refStack11(:,:,motFinal(i,3));
            refFrame15 = activityMovie15(:,:,i) + refStack15(:,:,motFinal(i,3));

            spDataClean(longSPs,i) = refFrame15(refPixs(longSPs)+motFinal(i,1)+motFinal(i,2)*dmdPixelsPerColumn) .* pbleach(i);
            spDataClean(shortSPs,i) = refFrame11(refPixs(shortSPs)+motFinal(i,1)+motFinal(i,2)*dmdPixelsPerColumn) .* pbleach(i);
        end
        
        spData = max(0, normrnd(poissrnd(spDataClean ./ 100) * 100,spDataClean/50));

        if ~exist([dr sprintf('simulation\\bright_%d',log2(bright))],'dir')
            mkdir([dr sprintf('simulation\\bright_%d',log2(bright))]);
        end
        
        fprintf('Saving brightness %d...',log2(bright));
        save([dr sprintf('simulation\\bright_%d\\simulated_SLAP2_data_activity.mat',log2(bright))],'sparseMaskInds','dmdPixelsPerRow','dmdPixelsPerColumn','numFastZs','spData','allSuperPixelIDs','motFinal')
        disp('done')
    end
end

