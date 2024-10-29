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
B = prctile(IMf(:,:,1:100), 10, 'all'); %estimate the mean brightness of a 'dim' pixel
IMf = log(IMf + B); %convert to log space so linear filtering computes products

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


% [~,P] = islocalmax2(summaryEroded); %, 'ProminenceWindow', 21);
% P(summaryEroded<median(summaryEroded, 'all', 'omitmissing')) = 0;
% 
% totalPix = sum(~isnan(summaryEroded(:)));
% [r,c,p] = find(P);
% [sortedP, sortorder] = sort(p, 'descend');
% 
% %find a prominence cutoff that corresponds to detecting a few too many
% %sources, then multiply that by a factor
% threshP = 1.5*sortedP(ceil(totalPix/81 * (1-exp(-nTimePoints*params.frametime*params.dsFac/10))));
% P(P<threshP) = 0;

if doPlot
    figure, imagesc(summaryEroded); %hAx1 = gca;
    hold on, scatter(P.col, P.row, 20*P.val, 'r' ); %'margeredgecolor', 'r');
    figure, imagesc(summaryEroded); %hAx2 = gca;
    %linkaxes([hAx1, hAx2]);
end

end

function plotSummary(summary, IMgamma, P)
red = summary./prctile(summary(~isnan(summary)), 99);
cyan = IMgamma./prctile(IMgamma(~isnan(IMgamma)),99);
figure, imshow(cat(3, red, cyan,cyan));
hold on,
scatter(P.col,P.row,100*(P.val./mean(P.val)).^2, 'm')
end

function P = getTiledPeaks(IM, IMavg, summary, params)
valid = ~imdilate(isnan(IM), ones(5)); 
peaks = valid;
for dim = 1:3
    peaks = peaks & IM > circshift(IM,1,dim) &  IM > circshift(IM,-1,dim);
end

linInds = find(peaks(:));
vvv = IM(linInds);
[rrr,ccc,ttt] = ind2sub(size(peaks), linInds);

%TO DO:
%We may want to discard localizations on dimmer pixels as inherently less
%likely, could use IMavg for this
%The selected region already selects for pixels with a large summary value
%vvv_ = vvv.*sqrt(IMavg(sub2ind(size(IMavg), rrr,ccc))); %signals normalized to image brightness

sz = size(IM);
tilesize = params.tilesizeLoc; %64;
tilestartsR = 1:tilesize/2:(sz(1)-tilesize/4);
tileendsR = min(sz(1), tilestartsR+tilesize-1);
tilestartsC = 1:tilesize/2:(sz(2)-tilesize/4);
tileendsC = min(sz(2), tilestartsC+tilesize-1);

summaryVals = summary(sub2ind(size(summary), rrr, ccc));

keep = false(1,length(ttt)); %which events to keep
vNorm = zeros(length(ttt),1); %the event sizes, Z-scored

threshSNR = params.threshSNRloc; %5; %the desired SNR as a z-score

for rix = 1:length(tilestartsR)
    for cix = 1:length(tilestartsC)
        selStats = rrr>=tilestartsR(max(1,rix-1)) & rrr<=tileendsR(min(end,rix+1))  & ccc>=tilestartsC(max(1,cix-1)) & ccc<=tileendsC(min(end,cix+1));

        S = summary(tilestartsR(max(1,rix-1)):tileendsR(min(end,rix+1)), tilestartsC(max(1,cix-1)):tileendsC(min(end,cix+1)));
        if any(S(:))
            Sp = prctile(S(:), [1 33]);
            Sthresh =  Sp(2) + params.threshSKloc*(Sp(2)-Sp(1));
            selS = summaryVals>Sthresh;

            vals = vvv(selStats & selS);
            ptile = prctile(vals, [1 50]);
            vals = 3*(vals-ptile(2))./(ptile(2) - ptile(1));
            thresh = ptile(2) + (threshSNR/3)*(ptile(2) - ptile(1)); % The [98% confint] here is about 3 sigma, hence threshSNR/3
            selTile = rrr>=tilestartsR(rix) & rrr<=tileendsR(rix)  & ccc>=tilestartsC(cix) & ccc<=tileendsC(cix) & vvv>thresh & selS; %the events within this tile that should be kept
            keep(selTile) = true;
            vNorm(selStats & selS) = max(vNorm(selStats & selS), vals);
        end
    end
end

rrr= rrr(keep);
ccc = ccc(keep);
ttt = ttt(keep);
vvv = vNorm(keep);

%upsample for superresolution
pC = []; pR = [];
for peakIx = length(ttt):-1:1
    R = IM(rrr(peakIx)+(-1:1), ccc(peakIx), ttt(peakIx));
    C = IM(rrr(peakIx), ccc(peakIx)+(-1:1), ttt(peakIx));

    ratioR = min(1e6,(R(2) - R(1))/(R(2) - R(3)));
    dR = (1-ratioR)/(1+ratioR)/2;
    pR(peakIx) = rrr(peakIx)-dR;

    ratioC = min(1e6,(C(2) - C(1))/(C(2) - C(3)));
    dC = (1-ratioC)/(1+ratioC)/2;
    pC(peakIx) = ccc(peakIx)-dC;
end

P.row = pR;
P.col = pC;
P.t = ttt;
P.val = vvv;

end

function IM = decorrelateMotion(IM, IMavg, aData,window)
nanFrames = squeeze(all(isnan(IM),[1 2]));
IMnan = isnan(IM);
sz = size(IM);
nPCs = 20;
goodPixels = mean(~IMnan,3)>0.9;
IM2 = IM;
IM2(IMnan) = 0;
IM2 = reshape(IM2(:,:,~nanFrames), sz(1)*sz(2), []);
IM2 = IM2- mean(IM2,2);
imageGrads = cat(3, IMavg-imtranslate(IMavg,[1 0]), IMavg-imtranslate(IMavg,[0 1]), IMavg-imtranslate(IMavg,[-1 0]), IMavg-imtranslate(IMavg,[0 -1]));
imageGrads = reshape(imageGrads, sz(1)*sz(2),[]);
[U,S,V] = svds(double(IM2(goodPixels(:),:)),nPCs);
b = imageGrads(goodPixels,:)\U; %mvregress(imageGrads(goodPixels,:),U);
correction = (imageGrads*b)*S*V';
IM(:,:, ~nanFrames) = IM(:,:, ~nanFrames) - reshape(correction, sz(1), sz(2), []);
end

function P = motionPredictors (aData, window)
%predictors are
C = aData.motionDSc;
R = aData.motionDSr;
C = (C-mean(C))./sqrt(mean(C.^2,'all', "omitnan"));
R = (R-mean(R))./sqrt(mean(R.^2,"all", 'omitnan'));

dC = diff(C); dC = ([0 dC] + [dC 0])/2;
dR = diff(R); dR = ([0 dR] + [dR 0])/2;

P = cat(1,ones(size(C)), dC, dC.^2,dR, dR.^2);
%if nargin>1 && ~isempty(window)
Cs = C-smoothdata(C', window, 'movmedian')';
Rs = R-smoothdata(R', window, 'movmedian')';
P = cat(1, P, Cs, Cs.^2,Rs, Rs.^2, Cs.*Rs);
%end
if isfield(aData, 'motionDSz')
    keyboard
end
end