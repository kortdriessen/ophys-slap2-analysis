function simulateActivity
%simulate a movie with some number of synapses to be retrieved using summarizeBergamo_Peaks
numChannels = 1;
rng(1);

%simulate data
cd('C:\Users\kaspar.podgorski\OneDrive - Allen Institute\Documents\GitHub\ophys-slap2-analysis\matlab\rasterDataProcessing\Bergamo\simulations\data');

fids = dir('*.tif');
fns = {fids.name};

darkrate = 0.02;
maxshift = 30;
IMsz = [132 60];
ds_time = 3;
dsFac = 2^ds_time;
clipShift = 5;
alpha = 0.0005; %exponential time constant for updating template
frametime = 0.0023;
brightness = 1;
bleachTau = 10000;
T = 10000;
motionAmp = 100;
tau = 8*1.33;
activityThresh = 0.12;
sigma = 1.33;
nsites = 80;
kernel = exp(-(0:ceil(8*tau))/tau);
minspike = 0.2;
maxspike = 5;
spikeAmp = 1;
sw = ceil(3*sigma);
skernel = zeros(2*sw+1);
skernel(sw+1,sw+1) = 1;
skernel = imgaussfilt(skernel, [sigma sigma]);
skernel = skernel./max(skernel(:));

for fix = 1:length(fns)
    fn = fns{fix};

    A = ScanImageTiffReader(fn);
    IMavg = mean(A.data,3,'omitmissing');
    %IMavg = IMavg(floor((size(IMavg,1)-IMsz(1))/2) +(1:IMsz(1)), floor((size(IMavg,2)-IMsz(2))/2) +(1:IMsz(2)));
    BG = prctile(IMavg(~isnan(IMavg)),30);
    IMavg = max(IMavg-BG, 0);
    IMavg = IMavg./prctile(IMavg(:), 99);
    selR = floor((size(IMavg,1)-IMsz(1))/2) +(1:IMsz(1));
    selC = floor((size(IMavg,2)-IMsz(2))/2) +(1:IMsz(2));

    %simulate synapses
    tmp = IMavg>min(prctile(IMavg, 97), 4*(mean(IMavg(:))));
    tmp([1:selR(1) selR(end):end], :) = false;
    tmp(:, [1:selC(1) selC(end):end]) = false;
    tmp = find(tmp);

    releaseSites = randsample(tmp, nsites);
    [rr,cc] = ind2sub(size(IMavg), releaseSites);
    dr = rand(size(rr))-0.5;
    dc = rand(size(cc))-0.5;

    GT.R = rr+dr-selR(1)+1; %in the coordinates of the saved image
    GT.C = cc+dc-selC(1)+1; %in the coordinates of the saved image
    for trialIx = 1:5
        aData = struct(); %extize;
        fnstem = ['SIMULATION_Trial' int2str(trialIx) fn(1:end-30)];

        B = brightness*exp(-(1:T)./bleachTau);
        activity = zeros(nsites, T);
        spikes = rand(size(activity))<smoothdata(rand(size(activity))<activityThresh,2, 'movmean', 40).^2;
        activity(spikes) =min(maxspike, max(minspike,spikeAmp*randn(1, sum(spikes(:)))));
        activity = convn(activity, kernel, 'same');

        movie = repmat(IMavg, [1 1 T]);
        for siteN = 1:nsites
            S = IMavg(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw));
            temp = S.*imtranslate(skernel, [dc(siteN) dr(siteN)]).* reshape(activity(siteN,:), 1,1,[]);
            movie(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw), :) = movie(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw), :) + temp;
        end
        movie = single(movie);

        envelope = sin(cumsum(randn(1,T)/20)).^2;
        motionPC1 =  smooth(envelope .* sin(smoothdata(randn(1,T).^3,2, 'movmean', 40)/10) .* motionAmp, 5);
        motionPC2 = smooth(envelope .* sin(smoothdata(randn(1,T).^3,2, 'movmean', 40)/10) .* motionAmp,5);
        GT.motionR = 0.8*motionPC1 + 0.4*motionPC2;
        GT.motionC = 0.2*motionPC1 - 0.2*motionPC2;

        for frameIx = T:-1:1
            tmp = imtranslate(movie(:,:,frameIx), [GT.motionC(frameIx), GT.motionR(frameIx)]);
            excessNoise = max(0.5, min(2, 1+randn(length(selR), length(selC))/2));
            Ad(:,:,1,frameIx) = poissrnd(tmp(selR,selC).*B(frameIx) + darkrate).*excessNoise;
        end
        Ad = reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []);
        sz = size(Ad);

        GT.activity = activity;

        T0 = padarray(IMavg(selR,selC), maxshift * [1 1], 0);
        template = T0;

        initR = 0; initC = 0;
        nDSframes= floor(sz(4)./(dsFac)); %number of downsampled frames
        motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
        aErrorDS = nan(1,nDSframes);
        [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshift))-maxshift, (1:(sz(2)+2*maxshift))-maxshift); %view matrices for interpolation

        disp('Registering:');
        for DSframe = 1:nDSframes
            readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));

            M = downsampleTime(Ad(:,:,:, readFrames), ds_time);
            M = squeeze(sum(M,3)); %merge colors

            if ~mod(DSframe, 1000)
                disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
            end

            Ttmp = mean(cat(3, T0,template),3, 'omitnan');
            T1 = Ttmp(maxshift-initR + (1:sz(1)), maxshift-initC+(1:sz(2)));

            %output = dftregistration(fft2(M),fft2(single(T)),4);
            output = dftregistration_clipped(fft2(M),fft2(single(T1)),4, clipShift);
            motionDSr(DSframe) = initR+output(3);
            motionDSc(DSframe) = initC+output(4);
            aErrorDS(DSframe) = output(1);

            if abs(motionDSr(DSframe))<maxshift && abs(motionDSc(DSframe))<maxshift
                A = interp2(1:sz(2), 1:sz(1), M,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
                sel = ~isnan(A);
                nantmp = sel & isnan(template);
                template(nantmp) = A(nantmp);
                template(sel) = (1-alpha)*template(sel) + alpha*(A(sel));
                initR = round(motionDSr(DSframe));
                initC = round(motionDSc(DSframe));
            end
        end

        %upsample the shifts and compute a tighter field of view
        tDS = (1:nDSframes).*(dsFac) - (2^(ds_time-1)) + 0.5;
        motionC = interp1(tDS, motionDSc, 1:((2^ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
        motionR = interp1(tDS, motionDSr, 1:((2^ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
        aError = interp1(tDS, aErrorDS, 1:((2^ds_time)*nDSframes), 'nearest', 'extrap');
        maxshiftC = max(abs(motionC)); maxshiftR = max(abs(motionR));
        [viewR, viewC] = ndgrid((1:(sz(1)+2*maxshiftR))-maxshiftR, (1:(sz(2)+2*maxshiftC))-maxshiftC); %view matrices for interpolation

        pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

        %save a downsampled aligned recording
        mkdir('./SIMULATIONS/');
        fnwrite = ['./SIMULATIONS/' fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];
        fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
        for DSframe = 1:nDSframes
            readFrames = (DSframe-1)*(dsFac) + (1:(dsFac));
            YY = downsampleTime(Ad(:,:,:, readFrames), ds_time);
            for ch = 1:numChannels
                B =  interp2(1:sz(2), 1:sz(1), YY(:,:,ch),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan)';
                fTIF.WriteIMG(single(B));
            end
        end
        fTIF.close;

        %save an original-time-resolution recording
        fnwrite = ['./SIMULATIONS/' fnstem '_REGISTERED_RAW.tif'];
        fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
        for frame = 1:length(motionC)
            for ch = 1:numChannels
                B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan)';
                fTIF.WriteIMG(single(B));
            end
        end
        fTIF.close;

        %save alignment metadata
        aData.numChannels = numChannels;
        aData.frametime = frametime;
        aData.motionC =  motionC;
        aData.motionR = motionR;
        aData.motionDSc = motionDSc;
        aData.motionDSr = motionDSr;
        aData.aError = aError;
        save(['./SIMULATIONS/' fnstem '_ALIGNMENTDATA.mat'], 'aData');
        save(['./SIMULATIONS/' fnstem '_SIMPARAMS.mat'], 'GT');
    end
end
end

function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
end