%working on solver for SLAP2 glutamate recordings
%TO DO
%Split out solver from simulation
%implement parallelization across subproblems
%
%Baseline estimation-
%   add regressors?
%
% NaN handling!!
clear GT params;
rng(0);

doParallel = false;
if doParallel
    delete(gcp('nocreate'));
    parpool('Processes'); %rather than Threads, needed for plotting
end

%Parameters
sz = [200,200];
num_sources = 100;
num_pixels = prod(sz);
num_time_points = 3000;
gen_Hz = 400;
tau_s = 0.03;
tau = gen_Hz*tau_s; %used for generating GT
k = [zeros(1, ceil(6*tau)) exp(-(0:(ceil(6*tau)))/tau)];
k = k./sum(k);

%inference parameters
params.microscope = 'simulation';
params.tau_full = tau;
params.dXY = 3;
params.selRadius = 2*params.dXY;
params.nmfIter = 2;
params.sigma_px = 1.33;
params.baselineWindow_Glu_s = 4;
params.denoiseWindow_s = 0.1;
params.analyzeHz = gen_Hz;
params.denoiseWindow_samps = params.denoiseWindow_s.*params.analyzeHz;
params.baselineWindow_samps = params.baselineWindow_Glu_s .* params.analyzeHz;

params = setParamsExtractTrial(params);

%simulate background
Bim = imgaussfilt(rand(sz), 2); Bim = max(0, Bim-prctile(Bim(:),20))+(0.1.*Bim);
GT.B = Bim .* reshape(1+exp(-(1:num_time_points)/num_time_points), [1,1,num_time_points]);
doMotion = true;
if doMotion
    Bmotion = smoothdata(imgaussfilt(randn([2*sz num_time_points]) , 4*params.selRadius+1), 3,"movmean",5); %fast, nonlocal noise
    Bmotion = Bmotion(ceil(sz(1)/2)+(1:sz(1)), ceil(sz(2)/2)+(1:sz(2)), :);
    Bmotion = 1+ 0.05.*Bmotion./std(Bmotion(:));
    GT.B = GT.B.*Bmotion;
    shiftR = rand(1,num_time_points)/4;
    shiftC = rand(1,num_time_points)/4;
    for frameN = 1:num_time_points
        GT.B(:,:,frameN) = imtranslate(GT.B(:,:,frameN), [shiftR(frameN) shiftC(frameN)]);
    end
end

%source locations
tmp = rand(sz).*Bim; sourceR = nan(1,num_sources); sourceC = nan(1, num_sources);
tmp([1 end],:) = 0;
tmp(:, [1 end]) = 0;
tmp = tmp.*islocalmax2(tmp);
[r,c,v] = find(tmp);
[sorted,sortorder] = sort(v, 'descend');
sources.R = r(sortorder(1:num_sources));
sources.C = c(sortorder(1:num_sources));

%selected pixels
selPix = false(sz);
selPix(sub2ind(sz, sources.R,sources.C)) = true;
selPix = imdilate(selPix, strel('disk', params.selRadius));

for ix = num_sources:-1:1
    tmp = zeros(sz);
    tmp(sources.R(ix), sources.C(ix)) = 1;
    GT.Hs(:,ix) = tmp(selPix);
    tmp = imgaussfilt(double(tmp), 1.5).*Bim;
    GT.H(:,ix) = tmp(selPix);
end

%simulate spikes
spikeprob = smoothdata(rand(num_sources,num_time_points),2, 'movmean', params.denoiseWindow_samps);
spikeprob = max(0,(spikeprob./mean(spikeprob) - 1));
spikeprob =0.1.*spikeprob./max(spikeprob(:));
spikes = rand(size(spikeprob)).*(rand(size(spikeprob))<spikeprob);
GT.S = 100*spikes;
%GT.S = 0.*GT.S; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%normalize H to maintain scale
normFac = sum(GT.H, 1);
GT.H = GT.H./normFac;
GT.Hs = GT.Hs./normFac;

GT.B = reshape(GT.B, [],num_time_points); GT.B = GT.B(selPix(:),:);
GT.X = convn(GT.S, k, 'same');
GT.Y = GT.B + GT.H*GT.X;

%simulate observations
F = 0.5.*(1+10*rand(size(GT.Y,1),1)) .* (1+rand(size(GT.Y))); %F=Freshness; a weighting factor proportional to the number of independent observations averaged into a measurement
Y_obs = poissrnd(F.*GT.Y)./F; %F=Freshness

%add in photon multiplier
mult = 70;
Y_obs = Y_obs.*mult;
GT.B = GT.B.*mult;
GT.X = GT.X.*mult;
GT.S = GT.S.*mult;
%SOLVE!
if doParallel
    resultFuture = extractTrial(Y_obs,1./F, sources, selPix, params,GT);
    [H,B,S,LS,F0,SNR] = fetchOutputs(resultFuture); 
    X = convn(S,params.k, 'same');
else
    [H,B,S,LS,F0,SNR] = extractTrial(Y_obs,1./F, sources, selPix, params,GT);
end
    
