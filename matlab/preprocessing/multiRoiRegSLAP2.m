function multiRoiRegSLAP2(fullPathToTrialTable, paramsIn)

if ~nargin
    [fn, dr] = uigetfile('*trialTable*.mat');
else
    [dr, fn, ext] = fileparts(fullPathToTrialTable); fn = [fn ext];
end

%PARAMETER SETTING
if nargin>1
    params = setParams('multiRoiRegSLAP2', paramsIn);
else
    params = setParams('multiRoiRegSLAP2');
end

params.startTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

%load the trial Table, which sets correspondences between the two DMDs
load([dr filesep fn], 'trialTable');

%set up parallelization
p = gcp('nocreate'); % If no pool, do not create new one.
if isempty(p)
    poolsize = 1;
else
    poolsize = p.NumWorkers;
end
nWorkers = min([params.nWorkers, numel(trialTable.filename), feature('numcores')]);
if poolsize<nWorkers || ~strcmpi(class(p), 'parallel.ProcessPool')
    delete(gcp('nocreate'));
    if params.nWorkers<15
        warning('You are using few parallel workers!');
    end
    parpool('processes',nWorkers); %limit the number of workers to avoid running out of RAM %4-30-24, lowering processes again to prevent another error (18 --> 15)
end
nDMDs = size(trialTable.filename,1);
nTrials = size(trialTable.trueTrialIx,2);

for DMD_ix = 1:nDMDs
    fnRegDS = cell(1,nTrials);
    fnAdata = cell(1,nTrials);
    firstLine = nan(1,nTrials);
    if nWorkers>1
        parfor f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix)]= alignAsync(dr, trialTable, params, f_ix, DMD_ix);
        end
    else
        for f_ix = 1:nTrials
            [fnRegDS{f_ix}, fnAdata{f_ix}, firstLine(f_ix)]= alignAsync(dr, trialTable, params, f_ix, DMD_ix);
        end
    end
    if params.isReVolt && ~isfield(trialTable, 'firstLineOriginal')
        trialTable.firstLineOriginal = trialTable.firstLine;
        tOffset = firstLine - trialTable.firstLineOriginal(DMD_ix,:);
        trialTable.firstLine = trialTable.firstLine + tOffset;
    end
    trialTable.fnRegDS(DMD_ix,:) = fnRegDS;
    trialTable.fnAdata(DMD_ix,:) = fnAdata;
end
%during alignment of some data we discard initial frames

params.endTime = char(datetime('now','TimeZone','local','Format','yyyy-MM-dd''T''HH:mm:ss.SSSZZZZZ'));

trialTable.alignParams = params;
save([dr filesep fn], "trialTable")

disp('done multiRoiRegistration.')
end

function [fnwrite, fnAdata, firstLine] = alignAsync(dr, trialTable, params, f_ix, DMD_ix)
if params.includeIntegrationROIs
    spTypeFlag = 0; %use all superpixel types
else
    spTypeFlag = 1; %use only raster superpixels
end

fn = trialTable.filename{DMD_ix,f_ix};
fnW = ['E' int2str(trialTable.epoch(f_ix)) 'T' int2str(f_ix) 'DMD' int2str(DMD_ix)];
firstLine = trialTable.firstLine(DMD_ix,f_ix);
lastLine = trialTable.lastLine(DMD_ix, f_ix);
aData = params;

disp(['Aligning: ' [dr filesep fn]])
fnwrite = [dr filesep fnW '_REGISTERED_DOWNSAMPLED-' int2str(aData.alignHz) 'Hz.tif'];
fnAdata = [dr filesep fnW '_ALIGNMENTDATA.mat'];

if ~params.overwriteExisting && exist(fnAdata, 'file') && exist(fnwrite, 'file')
    disp([fn ' is already aligned; skipping' newline 'To force realign, pass TRUE as second argument']);
    return
end

retry=0;
while retry<10
    try
        S2data = slap2.Slap2DataFile([dr filesep fn]);
        retry = 10;
    catch
        retry = retry+1
        continue
    end
end
meta = loadMetadata([dr filesep fn]);
linerateHz = 1/meta.linePeriod_s;
dt = linerateHz/aData.alignHz;
numChannels = S2data.numChannels;
minSamps = 10; %minimimum number of samples to include in template

if params.isReVolt
    numChannels = 1;
    redChannel =2;
    if isfield(trialTable,'firstLineOriginal') &&  ~isnan(trialTable.firstLineOriginal(DMD_ix,f_ix))
        %first line already adjusted
    else
        % load frames until you see the light turn on on channel 2
        f0 = firstLine+1000; minI = []; maxI = [];
        fEnd = round(0.8*firstLine + 0.2*lastLine);
        nSamps = 15;
        span = fEnd-f0;
        while (fEnd-f0)>(0.4*nSamps*dt)
            frames = round(linspace(f0,fEnd,nSamps));
            ii = nan(1,nSamps);
            for fix = 1:nSamps
                ii(fix) = mean(getImageWrapper(S2data, redChannel, frames(fix), ceil(dt), 1, spTypeFlag), 'all', 'omitnan'); %moving image Ch1
            end
            if isempty(minI)
                minI = min(ii, [], 'omitnan');
                maxI = max(ii, [], 'omitnan');
            end
            ixEnd = min(numel(ii), find(ii>(0.2*minI + 0.8*maxI), 1,'first')+1);
            ix0 = max(1,find(ii(1:ixEnd)<(0.65*minI + 0.35*maxI), 1, 'last')-2);
            if (frames(ixEnd)-frames(ix0)) >=span
                warning('trouble zooming in on time of light turn on.')
                break %stop zooming in
            else
                span = (frames(ixEnd)-frames(ix0));
            end
            if isempty(ix0) || isempty(ixEnd)
                warning(['isReVolt flag was set but laser on time could not be detected for file:' fn '. skipping...'])
                return
            end
            fEnd = frames(ixEnd); f0 = frames(ix0);
        end
        iMid = find(ii>(min(ii)+max(ii))/2, 1, 'first');
        firstLine = round(interp1(ii(iMid + [-1 0]), frames(iMid+[-1 0]), (min(ii)+max(ii))/2));
    end
end

%sanity checks
if isprop(S2data, 'hDataFile')
    assert(length(S2data.hDataFile.fastZs)==1); %single plane acquisitions only
    metaZ = S2data.hDataFile.fastZs;
else
    assert(length(S2data.hMultiDataFiles.fastZs)==1)
    metaZ = S2data.hMultiDataFiles.fastZs;
end

%%%%%Make an initial template
%crosscorrelate each initial frame to each other
disp('generating template')
initFrames =   round(linspace(firstLine, lastLine, 42)); initFrames = initFrames(2:end-1);
nInitFrames = length(initFrames);

if nInitFrames==0
    disp(['File ' fn ' was very short! Skipping alignment']);
    return
end
for fix = nInitFrames:-1:1
    for cix = numChannels:-1:1
        [Y(:,:,cix,fix), freshness(:,:,fix)] = getImageWrapper(S2data, cix, initFrames(fix), ceil(dt), 1, spTypeFlag);
    end
end

if params.isReVolt
    Y = squeeze(Y(:,:,1,:));
else
    Y = squeeze(sum(Y,3)); %use both channels;
end

%make data smaller for alignment
trimRows = find(~all(isnan(Y(:,:,1)),2), 1, 'first'):find(~all(isnan(Y(:,:,1)),2), 1, 'last');
trimCols = find(~all(isnan(Y(:,:,1)),1), 1, 'first'):find(~all(isnan(Y(:,:,1)),1), 1, 'last');
Y = Y(trimRows, trimCols,:); freshness = freshness(trimRows, trimCols,:);
sz = size(Y);

R = ones(nInitFrames);
motion = zeros(2,nInitFrames,nInitFrames);
for f1 = 1:nInitFrames
    for f2 = (f1+1):nInitFrames
        [motion(:,f1,f2), R(f1,f2)] = xcorr2_nans_weighted(Y(:,:,f2), freshness(:,:,f2), Y(:,:,f1), [0 ; 0], 3); 
        motion(:,f2,f1) = -motion(:,f1,f2);
        R(f2,f1) = R(f1,f2);
    end
end
[bestR, maxind] = max(median(R));
frameInds = find(R(:,maxind)>=bestR);

assert(aData.maxshift==round(aData.maxshift), 'params.maxshift must be an integer');
[viewR, viewC] = ndgrid((1:(sz(1)+2*aData.maxshift))-aData.maxshift, (1:(sz(2)+2*aData.maxshift))-aData.maxshift); %view matrices for interpolation
tFrames = nan(2*aData.maxshift+sz(1), 2*aData.maxshift+sz(2), numel(frameInds));
for fix = 1:numel(frameInds)
    tFrames(:,:,fix) = interpFrame(Y(:,:,frameInds(fix)), viewC(1,:)-motion(2,frameInds(fix), maxind), viewR(:,1)-motion(1,frameInds(fix), maxind), freshness(:,:,frameInds(fix)));
end
tSum = sum(tFrames,3, 'omitnan');
tN = sum(~isnan(tFrames),3);
template = sqrt(tSum./tN);
template(tN<minSamps) = nan;

if params.refStackTemplate
    if params.isReVolt
        error('refStack alignment not implemented for reVolt imaging')
    end
    fullTemplate = nan(size(trialTable.refStack{DMD_ix}.IM,[2 1]));
    fullTemplate((min(trimRows)-aData.maxshift):(max(trimRows)+aData.maxshift),(min(trimCols)-aData.maxshift):(max(trimCols)+aData.maxshift)) = template;

    if numel(trialTable.refStack{DMD_ix}.channels) == 2
        refStack = (trialTable.refStack{DMD_ix}.IM(:,:,1:2:end) + trialTable.refStack{DMD_ix}.IM(:,:,2:2:end))/2;
    else
        refStack = trialTable.refStack{DMD_ix}.IM;
    end
    templateShifts = xcorr2_nans(fullTemplate,refStack(:,:,floor(end/2)+1)',[0;0],aData.maxshift);
    T0 = imtranslate(permute(refStack,[2 1 3]),[templateShifts(2:-1:1),0]);
    T0 = T0((min(trimRows)-aData.maxshift):(max(trimRows)+aData.maxshift),(min(trimCols)-aData.maxshift):(max(trimCols)+aData.maxshift),:);
    disp('template generated from reference stack')
else
    T0 = template;
end

clear Y Yhp;

initR = 0; initC = 0;
DSframes = ceil(firstLine:dt:lastLine);
nDSframes= length(DSframes); %number of downsampled frames

%for spatial downsampling, used to calculate alignment quality and for quicker visualization
dsTimes = 2;
dsSz = floor(size(template)./(2^dsTimes));
A_ds = nan([dsSz nDSframes]);

motionDSr = nan(1,nDSframes);
motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
if params.refStackTemplate
    motionDSz = nan(1,nDSframes); %matrices to store the inferred motion
end
aErrorDS = nan(1,nDSframes); %alignment error output by dftregistration

%output TIF
pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);

V1 = nan(size(viewC,1),size(viewC,2),nDSframes,'single'); %variance factor; multiply the image value by this to get variance
registrationFailed = false;
%T = T0(aData.maxshift + (1:sz(1)), aData.maxshift+(1:sz(2)));
disp('Registering:');
try
    for DSframeIx = 1:nDSframes
        [M1, freshness] = getImageWrapper(S2data, 1, DSframes(DSframeIx), ceil(dt), 1, spTypeFlag); %moving image Ch1
        M1 = M1(trimRows, trimCols);
        freshness = freshness(trimRows, trimCols); %image freshness, the effective number of samples that contribute to the measurement
        if numChannels==2
            M2 =  getImageWrapper(S2data, 2, DSframes(DSframeIx), ceil(dt), 1, spTypeFlag); %moving image Ch2
            M2 = M2(trimRows, trimCols);
            M = sqrt(M1+M2);
        else
            M = sqrt(M1);
        end

        if ~mod(DSframeIx, 1000)
            disp([int2str(DSframeIx) ' of ' int2str(nDSframes)]);
        end

        if params.refStackTemplate
            T = T0(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)),:);
            [motOutput, corrCoeff] = xcorr2_nans3d(M, T, [0 ; 0], aData.clipShift);
            motionDSz(DSframeIx) = motOutput(3);
        else
            Ttmp = mean(cat(3, T0,template),3, 'omitnan');
            T = Ttmp(aData.maxshift-initR + (1:sz(1)), aData.maxshift-initC+(1:sz(2)));
            %[motOutput, corrCoeff] = xcorr2_nans(M, T, [0 ; 0], aData.clipShift);
            [motOutput, corrCoeff] = xcorr2_nans_weighted(M, freshness, T, [0 ; 0], aData.clipShift);
        end

        motionDSr(DSframeIx) = initR+motOutput(1);
        motionDSc(DSframeIx) = initC+motOutput(2);
        aErrorDS(DSframeIx) = 1-corrCoeff^2;

        %compute aligned image and variance factor
        [A1,V1(:,:,DSframeIx)] = interpFrame(M1, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);
        
        % %check for regrets
        % if numChannels==2
        %     A2 = interpFrame(M2, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);
        %     M = sqrt(A1+A2);
        % else
        %     M = sqrt(A1);
        % end
        % [motOutput2, corrCoeff2] = xcorr2_nans_weighted(M, 1./V1(:,:,DSframeIx), template, [0  ; 0], 3);
        % %[motOutput2, corrCoeff2] = xcorr2_nans(A1, template, [0  ; 0], 3);
        % if any(abs(motOutput2)>0.5)
        %     motOutput
        %     motOutput2
        %     DSframeIx  
        %     keyboard
        %     motionDSr(DSframeIx) = motionDSr(DSframeIx)+motOutput2(1);
        %     motionDSc(DSframeIx) = motionDSc(DSframeIx)+motOutput2(2);
        %     [A1,V1(:,:,DSframeIx)] = interpFrame(M1, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);
        % end

        fTIF.WriteIMG(single(A1));
        if numChannels==2
            A2 = interpFrame(M2, viewC(1,:)+motionDSc(DSframeIx), viewR(:,1)+motionDSr(DSframeIx), freshness);%interp2(1:sz(2), 1:sz(1), M2,viewC+motionDSc(DSframeIx), viewR+motionDSr(DSframeIx), 'linear', nan);
            fTIF.WriteIMG(single(A2));
            A =A1+A2;
        else
            A = A1;
        end

        %downsample in space
        dsTmp = A;
        for dsIx = 1:dsTimes
            dsTmp = dsTmp(1:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(1:2:2*floor(end/2), 2:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 1:2:2*floor(end/2)) + dsTmp(2:2:2*floor(end/2), 2:2:2*floor(end/2));
        end
        A_ds(:,:,DSframeIx) = dsTmp;

        sel = ~isnan(A);
        tSum(sel) = tSum(sel)+A(sel);
        % tSum = tSum*(1-aData.alpha);
        % tN = tN*(1-aData.alpha);
        tN(sel) = tN(sel)+1;
        template(sel) = sqrt(tSum(sel)./tN(sel));
        template(tN<minSamps) = nan;

        initR = round(motionDSr(DSframeIx));
        initC = round(motionDSc(DSframeIx));
    end
catch ME
    disp(ME);
    registrationFailed = true;
end

fTIF.close;

%Compute alignment error
offset = 10;
recNegErr = nan(1,nDSframes);
nChunks = ceil(nDSframes./(aData.alignHz*10)); %align 10 seconds at a time
chunkEdges = round(linspace(1, nDSframes+1, nChunks));
for chunkIx = 1:length(chunkEdges)-1
    t_ixs = chunkEdges(chunkIx):chunkEdges(chunkIx+1)-1;
    template = median(A_ds(:,:,t_ixs), 3, 'omitmissing');
    nanFrac = mean(isnan(A_ds(:,:,t_ixs)), 3);
    template(nanFrac>0.2) = nan;
    template_gamma = sqrt(max(0,template)+offset);
    template = repmat(template, 1,1,length(t_ixs));
    template(isnan(A_ds(:,:,t_ixs))) = nan;
    recNegErr(1,t_ixs) = sqrt(squeeze(mean((max(0, (template-A_ds(:,:,t_ixs))./template_gamma).^2), [1 2], 'omitnan')./mean(max(0,(template./template_gamma).^2, 'includenan'), [1 2], 'omitnan')));
end

if std(motionDSc)>1.5 || std(motionDSr)>1.5
    registrationFailed = true;
    warning(['Too much motion in file: ' fn]);
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
aData.varFacDS = V1; %variance factor; multiply pixel intensity by this to get a number proportional to the pixel's variance
if params.refStackTemplate
    aData.motionDSz = motionDSz;
end
aData.aError = aErrorDS;
aData.Z = metaZ;
%aData.aRankCorrDS = aRankCorrDS;
aData.recNegErr = recNegErr;
aData.cropRow = trimRows(1)-aData.maxshift; %offset to add to ROIs to index into original recording
aData.cropCol = trimCols(1)-aData.maxshift; %offset to add to ROIs to index into original recording

disp('Getting online motion correction offsets')
if isprop(S2data, 'hDataFile')
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hDataFile, DSframes);
else
    [aData.onlineXshift, aData.onlineYshift, aData.onlineZshift] = getOnlineMotion(S2data.hMultiDataFiles, DSframes);
end
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

save(fnAdata, 'aData', '-v7.3');
end

function meta = loadMetadata(datFilename)
ix = strfind(datFilename, 'DMD'+digitsPattern(1));
metaFilename = [datFilename(1:ix+3) '.meta'];
meta = load(metaFilename, '-mat');
end

function [IM, freshness] = getImageWrapper(S2data, channel, frames, dt, zPlane, spTypeFlag)
if spTypeFlag
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane, spTypeFlag);
else
    [IM,~,freshness] = S2data.getImage(channel, frames, dt, zPlane); %for backward compatibility
end
end