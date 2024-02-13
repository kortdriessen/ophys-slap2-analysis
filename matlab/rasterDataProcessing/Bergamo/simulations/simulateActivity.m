function simulateActivity
%simulate a movie with some number of synapses to be retrieved using summarizeBergamo_Peaks
numChannels = 1;
rng(1);

%simulate data
cd('C:\Users\kaspar.podgorski\OneDrive - Allen Institute\Documents\GitHub\ophys-slap2-analysis\matlab\rasterDataProcessing\Bergamo\simulations\data');

fids = dir('*.tif');
fns = {fids.name};

SimDescription = 'Standard2';

params.SimDescription = SimDescription;
params.darkrate = 0.02;
params.maxshift = 30;
params.IMsz = [132 60];
params.ds_time = 3;
params.dsFac = 2^params.ds_time;
params.clipShift = 5;
params.alpha = 0.0005; %exponential time constant for updating template
params.frametime = 0.0023;
params.brightness = 2;
params.bleachTau = 10000;
params.T = 10000;
params.motionAmp = 50;
params.tau = 8*1.33;
params.activityThresh = 0.12;
params.sigma = 1.33;
params.photonScale = 1000; % a multiplier from photons to digitizer signal, to test that analysis is invariant to scaling
params.nsites = 50;
params.minspike = 0.3;
params.maxspike = 4;
params.spikeAmp = 2;

kernel = exp(-(0:ceil(8*params.tau))/params.tau);

sw = ceil(3*params.sigma);
skernel = zeros(2*sw+1);
skernel(sw+1,sw+1) = 1;
skernel = imgaussfilt(skernel, [params.sigma params.sigma]);
skernel = skernel./max(skernel(:));

for fix = 1:length(fns)
    fn = fns{fix};

    A = ScanImageTiffReader(fn);
    IMavg = mean(A.data,3,'omitnan');
    %IMavg = IMavg(floor((size(IMavg,1)-params.IMsz(1))/2) +(1:params.IMsz(1)), floor((size(IMavg,2)-params.IMsz(2))/2) +(1:params.IMsz(2)));
    BG = prctile(IMavg(~isnan(IMavg)),30);
    IMavg = max(IMavg-BG, 0);
    IMavg = IMavg./prctile(IMavg(:), 99);
    selR = floor((size(IMavg,1)-params.IMsz(1))/2) +(1:params.IMsz(1));
    selC = floor((size(IMavg,2)-params.IMsz(2))/2) +(1:params.IMsz(2));

    %%%%simulate synapses
    %valid pixels for release sites
    tmp = medfilt2(IMavg, [3 3]);
    tmp = tmp>min(prctile(tmp(:), 97), 4*(mean(tmp(:))));
    tmp([1:selR(1) selR(end):end], :) = false;
    tmp(:, [1:selC(1) selC(end):end]) = false;
    tmp = find(tmp);

    releaseSites = randsample(tmp, params.nsites);
    [rr,cc] = ind2sub(size(IMavg), releaseSites);
    dr = rand(size(rr))-0.5; %subpixel location
    dc = rand(size(cc))-0.5; %subpixel location

    GT.R = rr+dr-selR(1)+1; %in the coordinates of the saved image
    GT.C = cc+dc-selC(1)+1; %in the coordinates of the saved image
    for trialIx = 1:5
        aData = struct(); %extize;
        fnstem = ['SIMULATION_' fn(1:11) SimDescription '_Trial' int2str(trialIx)];

        B = params.brightness*exp(-(1:params.T)./params.bleachTau);
        activity = zeros(params.nsites, params.T);
        spikes = rand(size(activity))<smoothdata(rand(size(activity))<params.activityThresh,2, 'movmean', 40).^2;
        activity(spikes) =min(params.maxspike, max(params.minspike,params.spikeAmp*randn(1, sum(spikes(:)))));
        activity = convn(activity, kernel, 'same');

        movie = repmat(IMavg, [1 1 params.T]);
        idealFilts = zeros([size(IMavg) params.nsites]);
        for siteN = params.nsites:-1:1
            S = IMavg(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw));
            sFilt = S.*imtranslate(skernel, [dc(siteN) dr(siteN)]);
            temp = sFilt.* reshape(activity(siteN,:), 1,1,[]);
            movie(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw), :) = movie(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw), :) + temp;
            idealFilts(rr(siteN)+(-sw:sw),cc(siteN)+ (-sw:sw), siteN) = sFilt;
        end
        movie = single(movie);

        envelope = sin(cumsum(randn(1,params.T)/20)).^2;
        motionPC1 =  smooth(envelope .* sin(smoothdata(randn(1,params.T).^3,2, 'movmean', 40)/10) .* params.motionAmp, 5);
        motionPC2 = smooth(envelope .* sin(smoothdata(randn(1,params.T).^3,2, 'movmean', 40)/10) .* params.motionAmp,5);
        GT.motionR = 0.8*motionPC1 + 0.4*motionPC2;
        GT.motionC = 0.2*motionPC1 - 0.2*motionPC2;

        for frameIx = params.T:-1:1
            tmp = imtranslate(movie(:,:,frameIx), [GT.motionC(frameIx), GT.motionR(frameIx)]);
            excessNoise = max(0.5, min(2, 1+randn(length(selR), length(selC))/2));
            Ad(:,:,1,frameIx) = poissrnd(tmp(selR,selC).*B(frameIx) + params.darkrate).*excessNoise.*params.photonScale;
        end
        Ad = reshape(Ad, size(Ad,1), size(Ad,2), numChannels, []);
        sz = size(Ad);

        T0 = padarray(IMavg(selR,selC), params.maxshift * [1 1], 0);
        template = T0;

        GT.activity = activity;     

        initR = 0; initC = 0;
        nDSframes= floor(sz(4)./(params.dsFac)); %number of downsampled frames
        motionDSr = nan(1,nDSframes); motionDSc = nan(1,nDSframes); %matrices to store the inferred motion
        aErrorDS = nan(1,nDSframes);
        aRankCorr = nan(1,nDSframes);
        [viewR, viewC] = ndgrid((1:(sz(1)+2*params.maxshift))-params.maxshift, (1:(sz(2)+2*params.maxshift))-params.maxshift); %view matrices for interpolation

        disp('Registering:');
        for DSframe = 1:nDSframes
            readFrames = (DSframe-1)*(params.dsFac) + (1:(params.dsFac));

            M = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
            M = squeeze(sum(M,3)); %merge colors

            if ~mod(DSframe, 1000)
                disp([int2str(DSframe) ' of ' int2str(nDSframes)]);
            end

            Ttmp = mean(cat(3, T0,template),3, 'omitnan');
            T1 = Ttmp(params.maxshift-initR + (1:sz(1)), params.maxshift-initC+(1:sz(2)));

            %output = dftregistration(fft2(M),fft2(single(T)),4);
            output = dftregistration_clipped(fft2(M),fft2(single(T1)),4, params.clipShift);
            motionDSr(DSframe) = initR+output(3);
            motionDSc(DSframe) = initC+output(4);
            aErrorDS(DSframe) = output(1);

            if abs(motionDSr(DSframe))<params.maxshift && abs(motionDSc(DSframe))<params.maxshift
                A = interp2(1:sz(2), 1:sz(1), M,viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan);
                sel = ~isnan(A);

                % A2 = A; A2(isnan(A2)) = mean(A2(:), 'omitnan');
                % T2 = template; T2(isnan(T2)) = mean(T2(:), 'omitnan');
                % A2 = imgaussfilt(A2,1);
                % T2 = imgaussfilt(T2,1);

                selCorr = ~(isnan(A) | isnan(template));
                aRankCorr(DSframe) = corr(A(selCorr), template(selCorr), 'type', 'Spearman');
                
                nantmp = sel & isnan(template);
                template(nantmp) = A(nantmp);
                template(sel) = (1-params.alpha)*template(sel) + params.alpha*(A(sel));
                initR = round(motionDSr(DSframe));
                initC = round(motionDSc(DSframe));
            end
        end

        %upsample the shifts and compute a tighter field of view
        tDS = (1:nDSframes).*(params.dsFac) - (2^(params.ds_time-1)) + 0.5;
        motionC = interp1(tDS, motionDSc, 1:((2^params.ds_time)*nDSframes), 'pchip','extrap'); %upsample the movement to the original timebase
        motionR = interp1(tDS, motionDSr, 1:((2^params.ds_time)*nDSframes), 'pchip', 'extrap'); %upsample the movement to the original timebase
        aError = interp1(tDS, aErrorDS, 1:((2^params.ds_time)*nDSframes), 'nearest', 'extrap');
        params.maxshiftC = ceil(max(abs(motionC))); params.maxshiftR = ceil(max(abs(motionR)));
        [viewR, viewC] = ndgrid((1:(sz(1)+2*params.maxshiftR))-params.maxshiftR, (1:(sz(2)+2*params.maxshiftC))-params.maxshiftC); %view matrices for interpolation

        pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

        tt = padarray(idealFilts(selR,selC,:), [params.maxshiftR params.maxshiftC 0], 0);
        GT.ROIs = tt;
        IF = reshape(tt, [],params.nsites); %ideal filters

        %save a downsampled aligned recording
        mkdir('./SIMULATIONS/');
        mkdir(['./SIMULATIONS/'  SimDescription filesep]);
        
        fnwrite = ['./SIMULATIONS/' SimDescription filesep fnstem '_REGISTERED_DOWNSAMPLED-' int2str(params.dsFac) 'x.tif'];
        fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
        for DSframe = 1:nDSframes
            readFrames = (DSframe-1)*(params.dsFac) + (1:(params.dsFac));
            YY = downsampleTime(Ad(:,:,:, readFrames), params.ds_time);
            for ch = 1:numChannels
                B =  interp2(1:sz(2), 1:sz(1), YY(:,:,ch),viewC+motionDSc(DSframe), viewR+motionDSr(DSframe), 'linear', nan)';
                fTIF.WriteIMG(single(B));
            end
        end
        fTIF.close;

        %save an original-time-resolution recording
        GT.ROI_activity = nan(params.nsites,params.T);
        fnwrite = ['./SIMULATIONS/' SimDescription filesep fnstem '_REGISTERED_RAW.tif'];
        fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
        for frame = 1:length(motionC)
            for ch = 1:numChannels
                B = interp2(1:sz(2), 1:sz(1), Ad(:,:,ch,frame),viewC+motionC(frame), viewR+motionR(frame), 'linear', nan)';
                fTIF.WriteIMG(single(B));

                tmp = reshape(B', 1,[]);
                for siteIx = 1:params.nsites
                    support = IF(:,siteIx)>0;
                    GT.ROI_activity(siteIx,frame) = tmp(1,support)*IF(support,siteIx); 
                end
            end
        end
        fTIF.close;

        gt{trialIx} = GT;
        %save alignment metadata
        aData.numChannels = numChannels;
        aData.frametime = params.frametime;
        aData.motionC =  motionC;
        aData.motionR = motionR;
        aData.motionDSc = motionDSc;
        aData.motionDSr = motionDSr;
        aData.aError = aError;
        aData.aRankCorr = aRankCorr;
        save(['./SIMULATIONS/' SimDescription filesep fnstem '_ALIGNMENTDATA.mat'], 'aData');
    end
    save(['./SIMULATIONS/' SimDescription filesep fnstem(1:end-7) '_SIMPARAMS.mat'], 'gt');
end

end

function Y = downsampleTime(Y, ds_time)
for ix = 1:ds_time
    Y = Y(:,:,:,1:2:(2*floor(end/2)))+ Y(:,:,:,2:2:end);
end
end