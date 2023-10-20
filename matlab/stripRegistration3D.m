function stripRegistration3D(ds_time) %% same as stripRegistration but aligning to a reference volume to allow for Z registration
maxshift = 80;
maxInitOffset = 80; %
clipShift = 15;%the maximum allowable shift per frame
zClipShift = 2; %max allowable shift per frame in Z (either up or down)

[stackfns, stackdr] = uigetfile('*.tif', 'Select Reference Volume', 'multiselect', 'off');

if strcmpi(stackfns(end-6:end-5),'CH')
    channel = str2num(stackfns(end-4));
else
    channel = 1;
end

refStack = tiffreadVolume([stackdr stackfns]);

refStackHP = refStack; % - imgaussfilt(refStack,4);

[fns, dr] = uigetfile('*.tif', 'Select Movie', 'multiselect', 'on');
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
    %     aErrorDS = nan(1,nDSframes);
    aCorrDS = nan(1,nDSframes);
    %     [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation


    fullField = padarray(refStackHP,maxshift * [1 1 0],'both');

    disp('Registering:');
    for DSframe = 1:nDSframes
        readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));

        M = downsampleTime(Ad(:,:,:, readFrames), ds_time);
        %         M = squeeze(sum(M,3)); %merge colors
        M = squeeze(M(:,:,channel,:));
        %M = M-imgaussfilt(M, 4); %highpass

        if ~mod(DSframe, 100)
            disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
        end

        if DSframe == 1
            padsz = size(refStack, 1:2) - sz(1:2);

            zRange = 3:size(refStackHP,3)-2;

            T3d = refStackHP;
            keep = false(2*size(T3d,1)-1, 2*size(T3d,2)-1);
            keep(ceil(end/2)+ [-maxInitOffset:maxInitOffset], ceil(end/2)+ [-maxInitOffset:maxInitOffset]) = true;
        else
            zRange = max(1,round(motionDSz(DSframe-1))-zClipShift):min(size(refStack,3),round(motionDSz(DSframe-1))+zClipShift);
            initR = max(-maxshift, min(maxshift, round(motionDSr(DSframe-1)))); 
            initC = max(-maxshift, min(maxshift, round(motionDSc(DSframe-1))));
            T3d = fullField((1:sz(1))+maxshift+floor(padsz(1)/2)+mod(padsz(1),2)+initR, (1:sz(2))+maxshift+floor(padsz(2)/2)+mod(padsz(2),2)+initC, :);
        end

        
        corrMetric = nan(1,max(zRange));
        for z = zRange
            T = T3d(:,:,z);

            if DSframe == 1
                c = normxcorr2(M,T);
                c(~keep) = nan; % restrict initial offset

                [rPeak, cPeak] = find(c==max(c(:)));

                initR = rPeak-(size(M,1)+floor(padsz(1)/2)+mod(padsz(1),2));
                initC = cPeak-(size(M,2)+floor(padsz(2)/2)+mod(padsz(2),2));

                T = fullField((1:sz(1))+maxshift+floor(padsz(1)/2)+mod(padsz(1),2)+initR, (1:sz(2))+maxshift+floor(padsz(2)/2)+mod(padsz(2),2)+initC,z);
            end

            [output, ~] = dftregistration_clipped(fft2(single(T)),fft2(M),4, clipShift);
            %shiftedFrame = real(ifft2(Greg));

            corrMetric(z) = 1-output(1);
            % subShiftedFrame = shiftedFrame(abs(shiftedFrame) > 1e-3);
            % subT = T(abs(shiftedFrame) > 1e-3);
            % 
            % corrMetric = corr2(subShiftedFrame,subT);

            if isnan(aCorrDS(DSframe)) || aCorrDS(DSframe) < corrMetric(z)
                aCorrDS(DSframe) = corrMetric(z);
                motionDSr(DSframe) = initR + output(3);
                motionDSc(DSframe) = initC + output(4);
                motionDSz(DSframe) = z;
            end
        end

        %superres Z estimation
        zz = motionDSz(DSframe);
        if zz>zRange(1) && zz<zRange(end)
            ratio = min(1e6,(corrMetric(zz) - corrMetric(zz-1))/(corrMetric(zz) - corrMetric(zz+1)));
            dZ = (1-ratio)/(1+ratio)/2;
            motionDSz(DSframe) = zz-dZ;
        end
    end

    %upsample the shifts and compute a tighter field of view
    tDS = (1:nDSframes).*(dsFac) - (2^(ds_time-1)) + 0.5;
    motionC = interp1(tDS, motionDSc, 1:((2^ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
    motionR = interp1(tDS, motionDSr, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
    motionZ = interp1(tDS, motionDSz, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
    %     aError = interp1(tDS, aErrorDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
    aCorr = interp1(tDS, aCorrDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
    maxshiftC = max(abs(motionC-motionC(1))); maxshiftR = max(abs(motionR-motionR(1)));
    [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshiftR))-maxshiftR, (1:(sz(2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation

    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

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

    % %save an original-time-resolution recording
    % fnwrite = [dr filesep fn(1:end-4) '_3DREGISTERED_RAW.tif'];
    % fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    % for frame = 1:length(motionC)
    %     for ch = 1:numChannels
    %         B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC-motionC(frame)+motionC(1), viewR-motionR(frame)+motionR(1), 'linear', nan)';
    %         fTIF.WriteIMG(single(B));
    %     end
    % end
    % fTIF.close;

    %save alignment metadata
    aData.numChannels = numChannels;
    aData.frametime = frametime;
    aData.motionC =  motionC;
    aData.motionR = motionR;
    aData.motionZ = motionZ;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.motionDSz = motionDSz;
    %     aData.aError = aError;
    aData.aCorr = aCorr;
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