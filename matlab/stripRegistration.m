function stripRegistration
maxshift = 25;
alpha = 0.01; %exponential time constant for template

[fns, dr] = uigetfile('*.tif', 'multiselect', 'on');
if ~iscell(fns)
    fns = {fns};
end

for f_ix = 1:length(fns)
    fn = fns{f_ix};
    A = ScanImageTiffReader([dr filesep fn]);
    Y_original = permute(single(A.data), [2 1 3]); % convert to single precision
    Y= Y_original;                 

    %compute the framerate from the metadata
    desc = A.descriptions;
    for frame = 10:-1:1
        meta = jsondecode(desc{frame});
        ts(frame) = meta.timestamp;
    end
    frametime = median(diff(ts));

    %downsample to align
    ds_time = 3; % the movie is downsampled using averaging in time by a factor of 2^ds_time
    ds_space = 0;% the movie is downsampled using averaging in space by a factor of 2^ds_space
    for ix = 1:ds_time
        Y = Y(:,:,1:2:(2*floor(end/2)))+ Y(:,:,2:2:end);
    end
    for ix = 1:ds_space
        Y = Y(:,1:2:(2*floor(end/2)),:)+ Y(:,2:2:end,:);
        Y = Y(1:2:(2*floor(end/2)),:,:)+ Y(2:2:end,:,:);
    end
    DSfacs = [2.^ds_space 2.^ds_space 2.^ds_time];
    szY = size(Y);

    %highpass in space
    Yhp = Y-imgaussfilt(Y, 4);

    %make an initial template with normcorr
    initFrames = 200;
    options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',initFrames,'max_shift',maxshift,'us_fac',4,'init_batch',initFrames, 'correct_bidir', false);
    F = normcorre(Yhp(:,:,1:initFrames),options_rigid);
    F = mean(F,3); %fixed image
    template = nan(2*maxshift+szY(1), 2*maxshift+szY(2));
    template(maxshift+(1:szY(1)), maxshift+(1:szY(2)))=F;
    T0 = template; T00 = zeros(size(template));

    nDSframes= szY(3); %number of downsampled frames
    %aligned = nan(2*maxshift+szY(1), 2*maxshift+szY(2), nDSframes);
    %alignedY = nan(2*maxshift+szY(1), 2*maxshift+szY(2), nDSframes);
    initR = 0; initC = 0;
    motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
    aErrorDS = nan(1,nDSframes);
    [viewR, viewC] = ndgrid((1:(szY(1)+2*maxshift))-maxshift, (1:(szY(2)+2*maxshift))-maxshift); %view matrices for interpolation
    for DSframe = 1:nDSframes 
        M = Yhp(:,:,DSframe);
        Ttmp = mean(cat(3, T0, T00,template),3, 'omitnan');
        T = Ttmp(maxshift-initR + (1:szY(1)), maxshift-initC+(1:szY(2)));
        
        output = dftregistration(fft2(M),fft2(single(T)),4);
        motionDSr(DSframe) = initR+output(3);
        motionDSc(DSframe) = initC+output(4);
        aErrorDS(DSframe) = output(1);

        A = interp2(1:szY(2), 1:szY(1), M,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
        %B = interp2(1:szY(2), 1:szY(1), Y(:,:,DSframe),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
        %aligned(:,:,DSframe) = A;
        %alignedY(:,:,DSframe) = B;
        
        sel = ~isnan(A);
        nantmp = sel & isnan(template);
        template(nantmp) = A(nantmp);
        template(sel) = (1-alpha)*template(sel) + alpha*(A(sel));
        
        initR = round(motionDSr(DSframe));
        initC = round(motionDSc(DSframe));
    end

    %upsample the shifts and compute a tighter field of view
    tDS = (1:nDSframes).*(2.^ds_time) - (2^(ds_time-1)) + 0.5; 
    motionC = interp1(tDS, motionDSc, 1:((2^ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
    motionR = interp1(tDS, motionDSr, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
    aError = interp1(tDS, aErrorDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
    maxshiftC = max(abs(motionC)); maxshiftR = max(abs(motionR));
    [viewR, viewC] = ndgrid((1:(size(Y_original,1)+2*maxshiftR))-maxshiftR, (1:(size(Y_original,2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation
    
    pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

    %save a downsampled aligned recording
    fnwrite = [dr filesep fn(1:end-4) '_REGISTERED_DOWNSAMPLED-' int2str(DSfacs(3)) 'x.tif'];
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    for DSframe = 1:length(motionDSc)
        B =  interp2(1:szY(2), 1:szY(1), Y(:,:,DSframe),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan); %single(interp2(1:size(Y_original,2), 1:size(Y_original,1), Y_original(:,:,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan));
        fTIF.WriteIMG(B');
    end
    fTIF.close;
    
    %save an original-time-resolution recording
    fnwrite = [dr filesep fn(1:end-4) '_REGISTERED_RAW.tif'];
    fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
    for frame = 1:length(motionC)
        B = interp2(1:size(Y_original,2), 1:size(Y_original,1), Y_original(:,:,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan);
        fTIF.WriteIMG(B');
    end
    fTIF.close;

    %save alignment metadata
    aData.frametime = frametime;
    aData.motionC =  motionC;
    aData.motionR = motionR;
    aData.motionDSc = motionDSc;
    aData.motionDSr = motionDSr;
    aData.aError = aError;
    save([dr filesep fn(1:end-4) '_ALIGNMENTDATA.mat'], 'aData');
end

end