trialDirs = {'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV1\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV2\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV3\',
             'Z:\ophys\SLAP2\exp data\Mice\slap2_664318_2023-08-01_11-07-05\FOV4\'};
brightnesses = [0.0625 0.125 0.25 0.5 1 2 4];
activity = true;

K = 10;
Nt = 3000;

frameRate = 100; % fps
spikeTau = 20; % ms
lambdas = [1 4 16 64 256]; % in spikes per s


for trial = 1%:length(trialDirs)
    fprintf('Trial %d\n',trial);
    dr = trialDirs{trial};

    % read reference stacks
    fns = dir([dr '*_CONFIG3-REFERENCE_CH2.tif']);
    refStack15 = tiffreadVolume([dr fns(1).name]);
    
    fns = dir([dr '*_CONFIG1-REFERENCE_CH2.tif']);
    refStack5 = tiffreadVolume([dr fns(1).name]);

    % load filter to go from dil 5 -> dil 15
    %   or create filter if doesn't exist
    if exist([dr 'Dil5To15Filter.mat'])
        load([dr 'Dil5To15Filter.mat']);
    else
        filterlucy = deconvlucy(refStack15,refStack5,300);

        % select central part of the filter
        filterfinal = filterlucy(392:410,637:645,11);
        % flip in all directions and average to make symmetric
        filterfinal = (filterfinal + filterfinal(end:-1:1,:) + filterfinal(:,end:-1:1) + filterfinal(end:-1:1,end:-1:1))/4;
        % scale filter intensity
        filterfinal = filterfinal / sum(filterfinal(:)) * sum(filterlucy(:));
        
        % display filter and wait to be checked
        figure; imagesc(filterfinal); daspect([1 1 1]); colormap gray
        pause();

        save([dr 'Dil5To15Filter.mat'], 'filterfinal');
    end

    load("Z:\ophys\SLAP2\exp data\Beads\slap2_000000_2023-06-15_14-07-03\dil15_psf_ds.mat");

    [dmdPixelsPerColumn, dmdPixelsPerRow, numZs] = size(refStack15);
    numFastZs = 1;

    % load selected superpixels
    load([dr 'simulation\simulationROI.mat']);
    [rasterROI,blockROI,rasterTargetHz,blockTargetHz] = rois.translateRoiShapes([dmdPixelsPerColumn dmdPixelsPerRow]);
    
    roiMasks = zeros([size(blockROI, 1:2) 0]);
    spData = [];
    allSuperPixelIDs = [];
    spCounter = 0;
    
    refRs = [];
    refCs = [];
    refPixs = [];
    
    % convert to SLAP2 format for roiMasks
    for roiIdx = 1:max(blockROI(:))
        currROIs = (blockROI == roiIdx);
        ROICols = find(sum(currROIs));
        tmpRoiMasks = zeros([800 1280 length(ROICols)]);
        spCounter = 0;
        for colIdx = ROICols
            spCounter = spCounter + 1;
            allSuperPixelIDs(end+1) = spCounter;
            tmpRoiMasks(:,colIdx,spCounter) = currROIs(:,colIdx);
    
            openPixs = find(tmpRoiMasks(:,:,spCounter));
    
            [refR, refC] = ind2sub(size(blockROI),median(openPixs));
            refRs = [refRs refR];
            refCs = [refCs refC];
    
            refPixs = [refPixs median(openPixs)];
        end
        roiMasks = cat(3, roiMasks, tmpRoiMasks);
    end

    numSuperPixels = size(roiMasks,3);
        
    % save sparseMaskInds as pixel position x superpixel ID
    [r,c] = find(reshape(roiMasks,[dmdPixelsPerColumn*dmdPixelsPerRow*numFastZs length(allSuperPixelIDs)]));
    sparseMaskInds = [r, c];

    % load motion vector
    tmp = load([dr 'simulation\simulatedMotion.mat']);
    motFinal = tmp.motFinal; clear tmp;
    
    for lambda = lambdas(1)
        if activity % make activity
            if ~exist([dr '\simulation\simulatedActivity_A.mat'])
                A15 = [];
                masks = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,numZs,K);
    
                % select and create ROIs for activity
                % display dil 5 reference stack
                fig = figure(121); clf;
                imshow3D(refStack5);
                hAx = fig.CurrentAxes;
            
                for k = 1:K % loop through K sources to simulate
                    point = drawpoint(hAx);
                    pause();

                    sourcePosition = round(point.Position);
                    currPlane = round(fig.Children(end-1).Value);

                    tmpMask = zeros(dmdPixelsPerColumn, dmdPixelsPerRow, numZs);
                    tmpMask(sourcePosition(2):sourcePosition(2)+1,sourcePosition(1):sourcePosition(1)+1,currPlane) = 1;
            
                    masks(:,:,:,k) = tmpMask;

                    clear('tmpMask');
            
                    % soften mask edges and make dil 5 footprint
                    % mask = imgaussfilt(masks(:,:,k));
                    % footprint5 = mask .* refStack5(:,:,9);
                    % 
                    % % convolve to get the dil 15 footprint
                    % footprint15 = convn(footprint5,filterfinal,'same');

                    footprint15 = convn(masks(:,:,:,k),psf15,'same');
                    footprint15 = footprint15 ./ footprint15(sourcePosition(2),sourcePosition(1),currPlane) .* refStack15(sourcePosition(2),sourcePosition(1),currPlane);
    
                    % save footprint
                    A15(:,end+1) = footprint15(:);
                end
                save([dr '\simulation\simulatedActivity_A.mat'],'A15','masks');
            else
                load([dr '\simulation\simulatedActivity_A.mat']);
            end

            if ~exist([dr sprintf('simulation\\simulatedActivity_phi_%d.mat',log(lambda)/log(4))])
                phi = [];

                 % generate a simulated spike train at 10x frequency
                for k = 1:K
                    spikes = zeros(Nt*10,1);
                    t = 1;
                    while t <= Nt*10
                        isi = round(random('Exponential',1/(lambda/frameRate/10)));
                        if t + isi <= Nt*10
                            spikes(t+isi) = 5;
                        end
                        t = t+isi;
                    end
                    
                    % convolve with the iglusnfr transient shape
                    transients = conv(spikes,exp(-(0:500)/(spikeTau * frameRate * 10 / 1000)));
    
                    % subsample and save glusnfr time trace
                    phi(:,end+1) = transients(1:10:Nt*10);
                end
                save([dr sprintf('simulation\\simulatedActivity_phi_%d.mat',log(lambda)/log(4))],'phi');
            else
                load([dr sprintf('simulation\\simulatedActivity_phi_%d.mat',log(lambda)/log(4))]);
            end
        end
    
    %     if activity
    %         motFinal(:,3) = 9 * ones(size(motFinal(:,3))); % override Z motion
    % %         motFinal(:,1:2) = zeros(size(motFinal(:,1:2))); % override XY motion
    %     end
        
        % simulate at different brightnesses
        for b = 5; %1:length(brightnesses)
            bright = brightnesses(b);
            
            % simulate photobleaching as a single exponential
            pbleach = bright .* exp(-(1:Nt)/1000);
            
            spDataClean = nan(numSuperPixels,Nt);
            movie15 = zeros(dmdPixelsPerColumn,dmdPixelsPerRow,Nt);
    
            if activity
                activityMovie15 = reshape(A15 * phi',[dmdPixelsPerColumn dmdPixelsPerRow numZs Nt]);
            else
                activityMovie15 = zeros([dmdPixelsPerColumn dmdPixelsPerRow Nt]);
            end
    
            for i = 1:Nt
                % add signal and background
                refFrame15 = activityMovie15(:,:,motFinal(i,3),i) + refStack15(:,:,motFinal(i,3));
                
                % full pixel resolution movie
                movie15(:,:,i) = refFrame15 .* pbleach(i);
    
                % downsampled superpixel recording
                spDataClean(:,i) = refFrame15(refPixs+motFinal(i,1)+motFinal(i,2)*dmdPixelsPerColumn) .* pbleach(i);
            end
            
            % add noise to the superpixel recording
            spData = max(0, normrnd(poissrnd(spDataClean ./ 100) * 100,spDataClean/50));
    
            if ~exist([dr sprintf('simulation\\bright_%d',log2(bright))],'dir')
                mkdir([dr sprintf('simulation\\bright_%d',log2(bright))]);
            end
    
            if ~exist([dr sprintf('simulation\\bright_%d\\activity_%d',log2(bright),log(lambda)/log(4))],'dir')
                mkdir([dr sprintf('simulation\\bright_%d\\activity_%d',log2(bright),log(lambda)/log(4))]);
            end
    
            fprintf('Saving brightness %d...',log2(bright));
            % save out simulated variables
            if activity
                save([dr sprintf('simulation\\bright_%d\\activity_%d\\simulated_SLAP2_data_activity.mat',log2(bright),log(lambda)/log(4))],'sparseMaskInds','dmdPixelsPerRow','dmdPixelsPerColumn','numFastZs','spData','allSuperPixelIDs','motFinal');
            else
                save([dr sprintf('simulation\\bright_%d\\simulated_SLAP2_data.mat',log2(bright))],'sparseMaskInds','dmdPixelsPerRow','dmdPixelsPerColumn','numFastZs','spData','allSuperPixelIDs','motFinal');
            end
            disp('done')
        end
    end
end

