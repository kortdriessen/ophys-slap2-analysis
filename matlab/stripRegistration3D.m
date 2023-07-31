function stripRegistration3D(ds_time) %% same as stripRegistration but aligning to a reference volume to allow for Z registration
tic;
maxshift = 80;
clipShift = 10;%the maximum allowable shift per frame
zClipShift = 2; %max allowable shift per frame in Z (either up or down)

[stackfns, stackdr] = uigetfile('*.tif', 'Select Reference Volume', 'multiselect', 'off');

refStack = tiffreadVolume([stackdr stackfns]);

refStackHP = refStack - imgaussfilt(refStack,4);

[fns, dr] = uigetfile('*.tif', 'Select Lookup Table', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

for f_ix = 1:length(fns)
    fn = fns{f_ix};
    disp(['Aligning: ' [dr filesep fn]])
    %tiff object to read metadata
    w = warning;
    warning('off', 'all')
    tiffObj = Tiff([dr filesep fn], 'r');
    tiffObj.setDirectory(1);
    for frame = 1:10 %compute the framerate from the metadata by reading a few frames
        meta = jsondecode(getTag(tiffObj,'ImageDescription'));
        timeStamp(frame) = meta.timestamp;
        tiffObj.nextDirectory;
    end
    numChannels = meta.numChannels;
    frametime = median(diff(timeStamp(1:numChannels:end)));
    clear tiffObj;
    warning(w);

    A = ScanImageTiffReader([dr filesep fn]);
    Ad = single(A.data);
    Ad = permute(reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []), [2 1 3 4]);
    sz = size(Ad);

    %tiffStack object to read data from disk memory mapped
    
    %disp('Reading TIFF headers...')
    %TIFFStack([dr filesep fn], [], [numChannels]); 
    %downsample to align
    if nargin<1
        ds_time = 3; % the movie is downsampled using averaging in time by a factor of 2^ds_time
    end
    dsFac = 2^ds_time;

    initR = 0; initC = 0;
    nDSframes= floor(sz(4)./(dsFac)); %number of downsampled frames
    motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); motionDSz = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes); % aCorrDS = nan(1,nDSframes);
%     [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation


    fullField = padarray(refStackHP,maxshift * [1 1 0],'both');
   
    disp('Registering:');
    for DSframe = 1:nDSframes
        readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));

        M = downsampleTime(Ad(:,:,:, readFrames), ds_time); 
        M = squeeze(sum(M,3)); %merge colors
        M = M-imgaussfilt(M, 4); %highpass

        if ~mod(DSframe, 100)
            disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
        end

        if DSframe == 1
            padsz = size(refStack, 1:2) - sz(1:2);
            M = padarray(M, floor(padsz/2),'both'); %pad M to match refStack
            M = padarray(M, mod(padsz,2),'pre');

            zRange = 1:size(refStackHP,3);

            T3d = refStackHP;
        else
            zRange = max(1,motionDSz(DSframe-1)-zClipShift):min(size(refStack,3),motionDSz(DSframe-1)+zClipShift);

            T3d = fullField((1:sz(1))+maxshift+floor(padsz(1)/2)+mod(padsz(1),2)+initR, (1:sz(2))+maxshift+floor(padsz(2)/2)+mod(padsz(2),2)+initC, :);
        end

        for z = zRange;
            T = T3d(:,:,z);

            if DSframe == 1
%                 [output, Greg] = dftregistration(fft2(single(T)),fft2(M),4);
                output = dftregistration(fft2(single(T)),fft2(M),4);
            else
%                 [output, Greg] = dftregistration_clipped(fft2(single(T)),fft2(M),4, clipShift);
                output = dftregistration_clipped(fft2(single(T)),fft2(M),4, clipShift);
            end

%             corrMetric = corr2(real(ifft2(Greg)),T); % using correlation gives same result as using error

            if isnan(aErrorDS(DSframe)) || aErrorDS(DSframe) > output(1) % aCorrDS(DSframe) < corrMetric
%                 aCorrDS(DSframe) = corrMetric;
                aErrorDS(DSframe) = output(1);
                motionDSr(DSframe) = initR + output(3);
                motionDSc(DSframe) = initC + output(4);
                motionDSz(DSframe) = z;
            end
        end
        
        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    end

    %upsample the shifts and compute a tighter field of view
    tDS = (1:nDSframes).*(dsFac) - (2^(ds_time-1)) + 0.5; 
    motionC = interp1(tDS, motionDSc, 1:((2^ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
    motionR = interp1(tDS, motionDSr, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
    motionZ = interp1(tDS, motionDSz, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
    aError = interp1(tDS, aErrorDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
%     aCorr = interp1(tDS, aCorrDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
    maxshiftC = max(abs(motionC-motionC(1))); maxshiftR = max(abs(motionR-motionR(1)));
    [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshiftR))-maxshiftR, (1:(sz(2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation
    
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

    toc;
    disp('Saving results');

    %save a downsampled aligned recording
    fnwrite = [dr filesep fn(1:end-4) '_3DREGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    for DSframe = 1:nDSframes
        readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));
        YY = downsampleTime(Ad(:,:,:, readFrames), ds_time); 
        for ch = 1:numChannels
            B =  interp2(1:sz(2), 1:sz(1), YY(:,:,ch),viewC-motionDSc(DSframe)+motionDSc(1), viewR-motionDSr(DSframe)+motionDSr(1), 'linear', nan)';
            fTIF.WriteIMG(single(B));
        end
    end
    fTIF.close;
    
    %save an original-time-resolution recording
    fnwrite = [dr filesep fn(1:end-4) '_3DREGISTERED_RAW.tif'];
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    for frame = 1:length(motionC)
        for ch = 1:numChannels
            B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC-motionC(frame)+motionC(1), viewR-motionR(frame)+motionR(1), 'linear', nan)';
            fTIF.WriteIMG(single(B));
        end
    end
    fTIF.close;

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = frametime;
    aData.motionC =  motionC;
    aData.motionR = motionR;
    aData.motionZ = motionZ;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.motionDSz = motionDSz;
    aData.aError = aError;
%     aData.aCorr = aCorr;
    save([dr filesep fn(1:end-4) '_3DALIGNMENTDATA.mat'], 'aData');
end

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