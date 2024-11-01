function [summaryEroded, P] = localizeSourcesSLAP2(IM, aData, params, doPlot)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
nTimePoints = size(IM,3);
tau = params.tau_s./(params.frametime*params.dsFac); %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels
baselineWindow = ceil(params.baselineWindow_Glu_s/(params.frametime*params.dsFac));
nans = isnan(IM);
IMavg = mean(IM,3, 'omitmissing');
IMgamma = sqrt(max(0,IMavg));

if nargin<4
    doPlot = false;
end

IMf= IM;
clear IM;

%time match filter
%IMf(nans) = 0;
mem = max(0,IMf(:,:,end));
gamma = exp(-1/tau);
for t = size(IMf,3):-1:1
    IMt = IMf(:,:,t);
    nanst = isnan(IMt);
    IMt(nanst) = mem(nanst);
    IMf(:,:,t) = gamma*mem + (1-gamma)*IMt;
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%log transform
B = prctile(IMf(:,:,1:100), 10, 'all') -  prctile(IMf(:,:,1:100), 1, 'all'); %estimate the mean brightness of a 'dim' pixel
IMf = log(max(IMf,0) + B); %convert to log space so linear filtering computes products

%remove all variance that looks like image movement - TESTED 10/17/24, does
%not have a significant effect
% IMf(nans) = nan;
% IMstruct = mean(IMf,3,'omitmissing'); 
% IMstruct(mean(isnan(IMf),3)>0.8) = min(IMstruct,[],'all','omitmissing');
% IMf = decorrelateMotion(IMf, IMstruct, aData, params);

%Difference of Gaussians
B = median(IMf(:,:,end-50:end), 'all', 'omitnan');%prctile(IMf(:,:,1:100), 10, 'all'); %median(IMf(:,:,end-50:end), 'all', 'omitnan');
IMf(nans) = B;
IMf = imgaussfilt(IMf, [sigma sigma]);
IMf = IMf - imgaussfilt(IMf, 5*[sigma sigma]);
IMf(nans) = nan;

%we performed filtering in the log space to perform multiplicatoins; return to non-log space
IMf = exp(IMf); 

%Highpass filter in time
IMf = IMf - smoothdata(IMf, 3, 'movmean', 2*baselineWindow, 'omitnan');  %- smoothdata(IMf, 3, 'movmedian', baselineWindow, 'omitnan'); 

%clip outliers
IMf = IMf-mean(IMf,3, 'omitmissing');
stdIM = std(IMf,0,3, 'omitmissing');
IMf = max(min(IMf, 6*stdIM, 'includemissing'), -6*stdIM, 'includemissing');

summary = skewness(IMf(:,:, 1:end-3*ceil(tau)), 1,3).*IMgamma; %remove the last few points, these can be outliers
valid = mean(nans,3)<0.33;
summaryEroded = summary;

summaryEroded(isnan(summaryEroded)) = median(summaryEroded,'all', 'omitmissing');
summaryEroded = summaryEroded - imgaussfilt(summaryEroded, 5*[sigma sigma]);
summaryEroded(~valid) = nan;
%summaryEroded(imdilate(isnan(summary), ones(3, 5))) = nan; %this removes odd phenomena at edges due to alignment, could probably be fixed by treating nans appropriately

%find local maxima
%peaks = ~imdilate(isnan(summaryEroded), ones(3));  %we won't find peaks at very edge of image
peaks = summaryEroded == ordfilt2(summaryEroded, 9, ones(3)); %> circshift(summaryEroded,1,dim) &  summaryEroded > circshift(summaryEroded,-1,dim);

%[r,c] = find(peaks);
p = summaryEroded(peaks);
sortedP = sort(p, 'descend');
totalPix = sum(~isnan(summaryEroded(:)));
threshP = 1.5*sortedP(ceil(totalPix/100 * (1-exp(-nTimePoints*params.frametime*params.dsFac/10))));
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
    %linkaxes([hAx1, hAx2]);
end

end
