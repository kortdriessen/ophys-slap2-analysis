function multiRoiRegistration3D(alignHz)
maxShift = 50;
zClipShift = 3;

if ~nargin || isempty(alignHz)
    alignHz = 33; %we will align data at this timescale, Hz
end

[refFn, refDr] = uigetfile('*REFERENCE*.tif', 'Select a Reference Image');
A = ScanImageTiffReader([refDr filesep refFn]);
if contains(refFn, '_CH1')
    alignChan = 1;
elseif contains(refFn, '_CH2')
    alignChan = 2;
else
    error('Nonstandard reference image name!')
end

ref= permute(A.data, [2 1 3]);

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

    nDSframes= floor(numLines/dt); %number of downsampled frames
    DSframes = (1:nDSframes)*dt;

    motionDSr = nan(1,nDSframes);
    motionDSc = nan(1,nDSframes); 
    motionDSz = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes);

    %compute things that only need to be computed once 
    Y =  S2data.getImage(alignChan, 1, dt, 1);
    trimRows = max(1,find(~all(isnan(Y(:,:,1)),2), 1, 'first') - maxShift):min(size(Y, 1), find(~all(isnan(Y(:,:,1)),2), 1, 'last')+maxShift); 
    trimCols = max(1,find(~all(isnan(Y(:,:,1)),1), 1, 'first')-maxShift):min(size(Y, 2),find(~all(isnan(Y(:,:,1)),1), 1, 'last')+maxShift); 
    Y = Y(trimRows, trimCols);
    sz = size(Y);
    [viewR, viewC] = ndgrid(1:sz(1), 1:sz(2)); %view matrices for interpolation
    MM = nan([sz 2]);

    %output TIF
    fnwrite = [dr filesep fn(1:end-4) '_REGISTERED_DOWNSAMPLED-' int2str(alignHz) 'Hz.tif'];
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM; 250nm
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);

    disp('Registering:');
    for DSframe = 1:nDSframes
        if DSframe == 1
            zRange = max(1, ceil(size(ref,3)/2)-4):(min(size(ref,3), ceil(size(ref,3)/2)+4));
            initR = 0; initC = 0; clipShift = 30;
        else
            clipShift = max(5, 30-5*Dsframe);
            zRange = max(1,round(motionDSz(DSframe-1))-zClipShift):min(size(ref,3),round(motionDSz(DSframe-1))+zClipShift);
            initR = max(-maxShift, min(maxShift, round(motionDSr(DSframe-1))));
            initC = max(-maxShift, min(maxShift, round(motionDSc(DSframe-1))));
            %T3d = ref((1:sz(1))+maxShift+floor(padsz(1)/2)+mod(padsz(1),2)+initR, (1:sz(2))+maxShift+floor(padsz(2)/2)+mod(padsz(2),2)+initC, :);
        end
        
        for ch = [1 2]
            tmp = S2data.getImage(ch, ceil(DSframe*dt), dt, 1); %moving image
            MM(:,:,ch) = tmp(trimRows, trimCols);
        end
        M = MM(:,:,alignChan);
       
        
        bestC = -1; corrCoeff = nan(1,max(zRange)); 
        for z = zRange
            if ~mod(DSframe, 1000)
                disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
            end

            T = ref(trimRows-initR,trimCols-initC,z);

            [motOutput, corrCoeff(z)] = xcorr2_nans(M, T, [0 ; 0], clipShift);
            
            if corrCoeff(z)>bestC
                bestC = corrCoeff(z);
                bestZ = z;
                motionDSr(DSframe) = initR+motOutput(1);
                motionDSc(DSframe) = initC+motOutput(2);
                aErrorDS(DSframe) = 1-corrCoeff(z)^2;
            end
        end
        
        %upsample z motion
        if any(bestZ==zRange([1 end]))
            motionDSz = bestZ;
        else %superresolution
             ratio = min(1e6,(corrCoeff(bestZ) - corrCoeff(bestZ-1))/(corrCoeff(bestZ) - corrCoeff(bestZ+1)));
             dZ = (1-ratio)/(1+ratio)/2;
             motionDSz(DSframe) = bestZ-dZ;
        end

        for ch = 1:2
            A = interp2(1:sz(2), 1:sz(1), MM(:,:,ch),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
            fTIF.WriteIMG(single(A));
        end
    end
    fTIF.close;

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = 1/alignHz;
    aData.DSframes = DSframes;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.motionDSz = motionDSz;
    aData.aError = aErrorDS;
    aData.cropRow = trimRows(1); %offset to add to ROIs to index into original recording
    aData.cropCol = trimCols(1); %offset to add to ROIs to index into original recording
    save([dr filesep fn(1:end-4) '_ALIGNMENTDATA.mat'], 'aData');
end

disp('done multiRoiRegistration.')
end

function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
end

function meta = loadMetadata(datFilename)
ix = strfind(datFilename, 'DMD'+digitsPattern(1));
metaFilename = [datFilename(1:ix+3) '.meta'];
meta = load(metaFilename, '-mat');
end