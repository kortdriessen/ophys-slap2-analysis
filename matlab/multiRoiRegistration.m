function multiRoiRegistration(alignHz, overwriteExisting)

if ~nargin || isempty(alignHz)
    alignHz = 50; %we will align data at this timescale, Hz
end
if nargin<2
    overwriteExisting = false; %by default, don't repeat alignment for files that are already processed
end
    
[fns, dr] = uigetfile('*.dat', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

aData.alignHz = alignHz;
aData.maxshift = 50;
aData.clipShift = 5;%the maximum allowable shift per frame
aData.alpha = 0.005; %exponential time constant for template

parfor f_ix = 1:length(fns)
    alignAsync(dr, fns{f_ix}, aData, overwriteExisting);
end

disp('done multiRoiRegistration.')
end

function alignAsync(dr, fn, aData, overwriteExisting)
disp(['Aligning: ' [dr filesep fn]])
    fnwrite = [dr filesep fn(1:end-4) '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
    fnAdata = [dr filesep fn(1:end-4) '_ALIGNMENTDATA.mat'];

    if ~overwriteExisting && exist(fnAdata, 'file') && exist(fnwrite, 'file')
        disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
        return
    end

    S2data = slap2.Slap2DataFile([dr filesep fn]);
    meta = loadMetadata([dr filesep fn]);
    linerateHz = 1/meta.linePeriod_s;
    dt = linerateHz/aData.alignHz;
    numChannels = S2data.numChannels;
    numLines = S2data.totalNumLines;
    
    %sanity checks
    assert(length(S2data.hDataFile.fastZs)==1); %single plane acquisitions only

    %%%%%Make an initial template
    %crosscorrelate each initial frame to each other
    disp('generating template')
    nInitFrames = min(40, floor(numLines/dt));
    for fix = nInitFrames:-1:1
        for cix = numChannels:-1:1
            Y(:,:,cix,fix) = S2data.getImage(cix, ceil(fix*dt), ceil(dt), 1);
        end
    end

    Y = squeeze(sum(Y,3));
    
    %make data smaller for alignment 
    trimRows = find(~all(isnan(Y(:,:,1)),2), 1, 'first'):find(~all(isnan(Y(:,:,1)),2), 1, 'last'); 
    trimCols = find(~all(isnan(Y(:,:,1)),1), 1, 'first'):find(~all(isnan(Y(:,:,1)),1), 1, 'last'); 
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
    
    [viewR, viewC] = ndgrid((1:(sz(1)+2*aData.maxshift))-aData.maxshift, (1:(sz(2)+2*aData.maxshift))-aData.maxshift); %view matrices for interpolation
    template = nan(2*aData.maxshift+sz(1), 2*aData.maxshift+sz(2));
    for fix = 1:length(frameInds)
        template(:,:,fix) = interp2(1:sz(2), 1:sz(1), Y(:,:,frameInds(fix)),viewC-motion(2,frameInds(fix), maxind), viewR-motion(1,frameInds(fix), maxind), 'linear', nan);
    end
    template = mean(template, 3); 
    T0 = template;
    clear Y Yhp;

    initR = 0; initC = 0;
    nDSframes= floor(numLines/dt); %number of downsampled frames
    DSframes = ceil((1:nDSframes)*dt);

    motionDSr = nan(1,nDSframes); 
    motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration
    aRankCorrDS = nan(1,nDSframes); %rank correlation, a better measure of alignment quality

    %output TIF
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);

    registrationFailed = false;
    disp('Registering:');
    try
    for DSframe = 1:nDSframes
        M1 = S2data.getImage(1, ceil(DSframe*dt), ceil(dt), 1); %moving image Ch1
        M1 = M1(trimRows, trimCols);
        if numChannels==2
            M2 =  S2data.getImage(2, ceil(DSframe*dt), ceil(dt), 1); %moving image Ch2
            M2 = M2(trimRows, trimCols);
            M = M1+M2;
        else
            M = M1;
        end
          
        if ~mod(DSframe, 1000)
            disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
        end

        Ttmp = mean(cat(3, T0,template),3, 'omitnan');
        T = Ttmp(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)));
        
        [motOutput, corrCoeff] = xcorr2_nans(M, T, [0 ; 0], aData.clipShift);
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
        template(sel) = (1-aData.alpha)*template(sel) + aData.alpha*(A(sel));
        
        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    end
    catch ME
        disp(ME);
        registrationFailed = true;
    end

    fTIF.close;

    if std(motionDSc)>5 || std(motionDSr)>5
        registrationFailed = true;
    end
    if registrationFailed
        disp(['REGISTRATION ERROR OCCURRED FOR FILE: ' fn newline 'YOU MAY NEED TO QC THIS FILE!' newline 'CONTINUING...'])
        return
    end

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = 1/aData.alignHz;
    aData.DSframes = DSframes;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.aError = aErrorDS;
    aData.aRankCorrDS = aRankCorrDS;
    aData.cropRow = trimRows(1)-aData.maxshift; %offset to add to ROIs to index into original recording
    aData.cropCol = trimCols(1)-aData.maxshift; %offset to add to ROIs to index into original recording

    %CONVERTING DATAFILE IMAGES INTO THE SAVED TIFF IMAGE SPACE:
    aData.trimRows = trimRows; %used to remap images from the datafile into the space of the saved tiffs
    aData.trimCols = trimCols;%used to remap images from the datafile into the space of the saved tiffs
    aData.viewC = viewC;%used to remap images from the datafile into the space of the saved tiffs
    aData.viewR = viewR;%used to remap images from the datafile into the space of the saved tiffs
    %EXAMPLE CODE
    %  Y = S2data.getImage(channel, lineTime, deltaTime, zPos);
    %  Ytrimmed = Y(aData.trimRows, aData.trimCols);
    %  sz = size(Ytrimmed);
    %  Yshifted = interp2(1:sz(2), 1:sz(1), Ytrimmed,aData.viewC+motionC, aData.viewR+motionR, 'linear', nan);

    save(fnAdata, 'aData');
end

function meta = loadMetadata(datFilename)
    ix = strfind(datFilename, 'DMD'+digitsPattern(1));
    metaFilename = [datFilename(1:ix+3) '.meta'];
    meta = load(metaFilename, '-mat');
end