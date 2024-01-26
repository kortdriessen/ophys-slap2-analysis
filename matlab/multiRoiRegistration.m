function multiRoiRegistration(alignHz)
maxshift = 80;
clipShift = 5;%the maximum allowable shift per frame
alpha = 0.005; %exponential time constant for template
if ~nargin || isempty(alignHz)
    alignHz = 33; %we will align data at this timescale, Hz
end
    
[fns, dr] = uigetfile('*.dat', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

for f_ix = 1:length(fns)
    fn = fns{f_ix};
    disp(['Aligning: ' [dr filesep fn]])
    
    S2data = slap2.Slap2DataFile([dr filesep fn]);
    meta = loadMetadata([dr filesep fn]);
    linerateHz = 1/meta.linePeriod_s;
    dt = ceil(linerateHz/alignHz);
    numChannels = S2data.numChannels;
    numLines = S2data.totalNumLines;
    
    %sanity checks
    assert(length(S2data.hDataFile.fastZs)==1); %single plane acquisitions only

    %imageData = S2data.getImage(cIx, tIx, dt, zIx);
    %S2data.hDataFile.numCycles

    % A = ScanImageTiffReader([dr filesep fn]);
    % Ad = single(A.data);
    % Ad = permute(reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []), [2 1 3 4]);
    % sz = size(Ad);

    % %temporary; in future need to set invalid data as -1 rather than 0 in the tiff files...
    % nanMask = imopen(all(all(Ad<=0,3),4), ones(3));
    % Ad(repmat(nanMask, [1 1 numChannels size(Ad,4)])) = nan;

    % %downsample to align
    % if nargin<1
    %     %compute downsampling factor
    %     % the movie is downsampled using averaging in time by a factor of 2^ds_time
    %     target_Hz = 20; %downsample to a minimum speed of 25 Hz
    %     ds_time = max(0, floor(log2(1/frametime/target_Hz))); 
    % end
    % dsFac = 2^ds_time;


    %%%%%Make an initial template
    %crosscorrelate each initial frame to each other
    disp('generating template')
    nInitFrames = min(40, floor(numLines/dt));
    for fix = nInitFrames:-1:1
       Y(:,:,1,fix) = S2data.getImage(1, ceil(fix*dt), dt, 1);
       Y(:,:,2,fix) = S2data.getImage(2, ceil(fix*dt), dt, 1);
    end

    Y = squeeze(sum(Y,3));
    
    %make data smaller for alignment 
    trimRows = [find(~all(isnan(Y(:,:,1)),2), 1, 'first'):find(~all(isnan(Y(:,:,1)),2), 1, 'last')]; 
    trimCols = [find(~all(isnan(Y(:,:,1)),1), 1, 'first'):find(~all(isnan(Y(:,:,1)),1), 1, 'last')]; 
    Y = Y(trimRows, trimCols,:);
    sz = size(Y);

    R = ones(nInitFrames);
    motion = zeros(2,nInitFrames,nInitFrames);
    for f1 = 1:nInitFrames
        for f2 = (f1+1):nInitFrames
            [motion(:,f1,f2), R(f1,f2)] = xcorr2_nans(Y(:,:,f2), Y(:,:,f1), [0 ; 0], 3);
            motion(:,f2,f1) = -motion(:,f1,f2);
            R(f2,f1) = R(f1,f2);
        end
    end
    [bestR, maxind] = max(median(R));
    frameInds = find(R(:,maxind)>=bestR);
    
    [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation
    template = nan(2*maxshift+sz(1), 2*maxshift+sz(2));
    for fix = 1:length(frameInds)
        template(:,:,fix) = interp2(1:sz(2), 1:sz(1), Y(:,:,frameInds(fix)),viewC-motion(2,frameInds(fix), maxind), viewR-motion(1,frameInds(fix), maxind), 'linear', nan);
    end
    template = mean(template, 3); 
    T0 = template;
    clear Y Yhp;

    initR = 0; initC = 0;
    nDSframes= floor(numLines/dt); %number of downsampled frames
    DSframes = (1:nDSframes)*dt;

    motionDSr = nan(1,nDSframes); 
    motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration
    aRankCorrDS = nan(1,nDSframes); %rank correlation, a better measure of alignment quality

    %output TIF
    fnwrite = [dr filesep fn(1:end-4) '_REGISTERED_DOWNSAMPLED-' int2str(alignHz) 'Hz.tif'];
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);

    disp('Registering:');
    for DSframe = 1:nDSframes
        M1 = S2data.getImage(1, ceil(DSframe*dt), dt, 1); %moving image Ch1
        M1 = M1(trimRows, trimCols);
        if numChannels==2
            M2 =  S2data.getImage(2, ceil(DSframe*dt), dt, 1); %moving image Ch2
            M2 = M2(trimRows, trimCols);
            M = M1+M2;
        else
            M = M1;
        end
          
        if ~mod(DSframe, 1000)
            disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
        end

        Ttmp = mean(cat(3, T0,template),3, 'omitnan');
        T = Ttmp(maxshift-initR + (1:sz(1)), maxshift-initC+(1:sz(2)));
        
        [motOutput, corrCoeff] = xcorr2_nans(M, T, [0 ; 0], clipShift);
        motionDSr(DSframe) = initR+motOutput(1);
        motionDSc(DSframe) = initC+motOutput(2);
        aErrorDS(DSframe) = 1-corrCoeff^2;

        A1 = interp2(1:sz(2), 1:sz(1), M1,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
        fTIF.WriteIMG(single(A1));
        if numChannels==2
            A2 = interp2(1:sz(2), 1:sz(1), M2,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
            fTIF.WriteIMG(single(A2));
            A =A1+A2;
        else
            A = A1;
        end

        sel = ~isnan(A);
        selCorr = sel & ~isnan(template);
        aRankCorrDS(DSframe) = corr(A(selCorr), template(selCorr), 'type', 'Spearman');

        nantmp = sel & isnan(template);
        template(nantmp) = A(nantmp);
        template(sel) = (1-alpha)*template(sel) + alpha*(A(sel));
        
        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    end
    fTIF.close;

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = 1/alignHz;
    aData.DSframes = DSframes;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.aError = aErrorDS;
    aData.aRankCorrDS = aRankCorrDS;
    aData.cropRow = trimRows(1)-maxshift; %offset to add to ROIs to index into original recording
    aData.cropCol = trimCols(1)-maxshift; %offset to add to ROIs to index into original recording
    save([dr filesep fn(1:end-4) '_ALIGNMENTDATA.mat'], 'aData');
end

disp('done multiRoiRegistration.')
end

function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
end

function [Y, readsuccess, done] = readFrames(tiffObj, nFramesToRead, nChannels)
nReads = nFramesToRead*nChannels;
done = false;
Y = single(tiffObj.read);
Y(:,:,2:nReads) = nan;
try
    for r = 2:nReads
        tiffObj.nextDirectory;
        Y(:,:,r) = tiffObj.read;
    end
catch ME
    readsuccess = false;
    return
end
Y = reshape(Y, size(Y,1), size(Y,2), nChannels, nFramesToRead);
readsuccess = true;
try
    tiffObj.nextDirectory;
catch
    disp('reached End of File');
    done = true;
end
end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end