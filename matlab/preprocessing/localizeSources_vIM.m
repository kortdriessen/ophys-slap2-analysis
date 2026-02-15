function [summaryEroded, P] = localizeSources_vIM(IM, vIM, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
nTimePoints = size(IM,3);
tau = params.tau_s.*params.alignHz; %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels
baselineWindow = ceil(params.baselineWindow_Glu_s.*params.alignHz);
denoiseWindow = ceil(params.denoiseWindow_s.*params.alignHz);
nans = isnan(IM);

sz = size(IM);
valid = mean(nans,3)<params.nanThresh; %a pixel must be imaged at least (1-nanThresh) of the time to be included
if ~any(valid)
    warning('Recording had no valid pixels; likely too much motion')
    P = [];
    summaryEroded = nan(sz(1:2));
    return
end
if nargin<4
    doPlot = false;
end

firstValidFrames = find(any(~isnan(IM), [1 2]),400, 'first');
varIM = var(IM(:,:,firstValidFrames),0,3,"omitmissing");

%initialize filtered image
IMf = IM; clear IM;
IMf(repmat(~valid, 1, 1, nTimePoints)) = nan;
nans = isnan(IMf);

%fill in missing values
IMf= reshape(IMf, sz(1)*sz(2), []);
nanFrac = mean(isnan(IMf),2);
incomplete = nanFrac>0 & nanFrac<1;
IMs = nan(size(IMf));
IMs(incomplete,:) = smoothdata(IMf(incomplete,:), 2, 'movmean', baselineWindow, 'omitnan');
IMf(reshape(nans, size(IMf))) = IMs(reshape(nans, size(IMf))); clear IMs
IMf = reshape(IMf, sz(1),sz(2), []);

% IMb = reshape(computeF0(reshape(smoothdata(IMf, 3, 'movmean', denoiseWindow, 'omitnan'),[],size(IMf,3))',denoiseWindow,baselineWindow),size(IMf));
% IMb = smoothdata(smoothdata(IMf, 3, 'movmean', denoiseWindow, 'omitnan'), 3, 'movmedian', baselineWindow, 'omitnan');  

if params.microscope == "SLAP2"
    %Highpass filter in time; This must occur before DoG to avoid edge artifacts
    IMb = smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan');
    IMf = IMf - IMb;   %- smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan');

    if isempty(vIM) %for raster imaging
        vIM = ones(sz(1:2));
    end
    varIM(nanFrac>0.4) = nan;
    
    Vb = prctile(varIM, 1, 'all'); %estimate the variance of a 'dim' pixel, due to electronic and dark noise

    stdIM = sqrt((100*IMb+Vb).*vIM); %compute standard deviation
else
    IMfden = smoothdata(IMf, 3, 'movmean', denoiseWindow, 'omitnan');
    %Highpass filter in time; This must occur before DoG to avoid edge artifacts
    IMb = smoothdata(IMfden, 3, 'movmedian', baselineWindow, 'omitnan');
    IMf = IMf - IMb;   %- smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan');

    % MAD-based robust standard deviation estimate
    stdIM = movmad(IMfden - IMb,baselineWindow,3,'omitmissing') ./ 0.6741891400433162.*denoiseWindow;
end
%divide by uncertainty to get a Z-score
IMf = IMf./stdIM;
clear IMb

%time matched filter
gamma = exp(-1/tau);
mem = max(0,gamma*IMf(:,:,end));
for t = size(IMf,3):-1:1
    IMt = IMf(:,:,t);
    nanst = isnan(IMt);
    IMt(nanst) = mem(nanst);
    IMf(:,:,t) = gamma*mem + (1-gamma)*IMt;
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%Difference of Gaussians
IMf(nans) = 0;
IMf = imgaussfilt(IMf, [sigma sigma]);
IMf = IMf - imgaussfilt(IMf, 5*[sigma sigma]);
IMf(nans) = nan;
clear nans

%nonmax suppression- find maxima
skIm = zeros(sz(1:2));
for fr = size(IMf,3)-ceil(1.5*tau):-1:2 %ceil(tau) because the filtering is uncertain in the final frames
    
    IMfr = IMf(:,:,fr); IMpre = IMf(:,:,fr-1); IMpost = IMf(:,:,fr+1);
    
    selMax = IMfr==ordfilt2(IMfr,9, ones(3));
    IMlocalMax(:,:,fr) = selMax & IMfr>IMpre & IMfr>=IMpost;
    %maxinds = find(IMfr==ordfilt2(IMfr,9, ones(3)));
    % sel = IMfr(maxinds)>0 & IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % %sel = IMpre(maxinds)<=IMfr(maxinds) & IMpost(maxinds)<=IMfr(maxinds);
    % maxinds = maxinds(sel);
    skIm(IMlocalMax(:,:,fr)) = skIm(IMlocalMax(:,:,fr)) + IMfr(IMlocalMax(:,:,fr)).^2; 
end


%summary = skewness(IMf(:,:, 1:end-3*ceil(tau)), 1,3); %.*IMgamma; 
summaryEroded = skIm;
summaryEroded(~valid) = nan;
summaryEroded(isnan(summaryEroded)) = median(summaryEroded,'all', 'omitmissing');
summaryEroded = summaryEroded - medfilt2(summaryEroded, [5 5]);
summaryEroded(~valid) = nan;

%find local maxima
peaks = summaryEroded == ordfilt2(summaryEroded, 9, ones(3)); %> circshift(summaryEroded,1,dim) &  summaryEroded > circshift(summaryEroded,-1,dim);

p = summaryEroded(peaks);
sortedP = sort(p, 'descend');
totalPix = sum(~isnan(summaryEroded(:)));

if totalPix<10 || sum(peaks(:))==0
    P.row = [];
    P.col = [];
    P.val = [];
    P.peakIM = zeros(size(summaryEroded));
    return
end

threshP = 1.5*sortedP(min(end,ceil(totalPix * params.maxSynapseDensity * (1-exp(-nTimePoints./params.alignHz./10)))));
pp = summaryEroded; pp(~peaks) = 0; pp(pp<threshP) = 0;
[rrr,ccc,vvv] = find(pp);

%upsample for superresolution
pC = []; pR = [];
for peakIx = length(vvv):-1:1
    R = summaryEroded(rrr(peakIx)+(-1:1), ccc(peakIx));
    C = summaryEroded(rrr(peakIx), ccc(peakIx)+(-1:1));

    ratioR = min(1e6,(R(2) - R(1))/(R(2) - R(3)));
    dR = (1-ratioR)/(1+ratioR)/2;
    pR(peakIx) = rrr(peakIx)-dR;

    ratioC = min(1e6,(C(2) - C(1))/(C(2) - C(3)));
    dC = (1-ratioC)/(1+ratioC)/2;
    pC(peakIx) = ccc(peakIx)-dC;
end
P.row = pR(:);
P.col = pC(:);
P.val = vvv(:);
P.peakIM = pp;

if doPlot
    figure, imagesc(summaryEroded); %hAx1 = gca;
    hold on, scatter(P.col, P.row, 20*P.val, 'r' ); %'margeredgecolor', 'r');
    figure, imagesc(summaryEroded); %hAx2 = gca;
end

end
