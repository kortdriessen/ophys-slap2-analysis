function [summary, P, params] = localizeFlashesBergamo(IM, aData, params)
%inputs:
%IM:        3D recording, X x Y x Time
%aData:     alignment metadata
tau = params.tau_s./(aData.frametime*aData.dsFac); %time constant in frames
params.tau_frames = tau;
sigma = params.sigma_px; %space constant in pixels

denoiseWindow = params.denoiseWindow_samps;
baselineWindow = ceil(params.baselineWindow_Glu_s/(aData.frametime*aData.dsFac));

%nanFrames = squeeze(all(isnan(IM),[1 2]));
nans = isnan(IM);
%nansC = convn(nans, true(1,1,7), 'same')>0;
%IM(nansC) = nan;
IMavg = mean(IM,3, 'omitmissing');
IMgamma = sqrt(max(0,IMavg));

BG = prctile(IMavg(~isnan(IMavg)), 10);

%subtract background, relevant for Bergamo only (not SLAP2)
IMf = IM-BG;
IMf(nans) = 0;

%create a Difference of Gaussians filter (mean 0) that highights spots
DoGfilt = zeros(2*ceil(6.5*sigma)+1);
DoGfilt(ceil(end/2), ceil(end/2)) = 1;
DoGfilt = imgaussfilt(DoGfilt, [sigma sigma])-imgaussfilt(DoGfilt, 3*[sigma sigma]);
DoGfilt = DoGfilt - mean(DoGfilt(:));

%Apply spatial filter
IMf = convn(IMf,DoGfilt, 'same');

%Highpass filter in time
IMf(nans) = nan;
IMf = IMf - permute(computeF0(permute(IMf, [3 1 2]), denoiseWindow, baselineWindow,2), [2 3 1]);
nans = isnan(IMf);

%remove motion-associated variance
%IMf = decorrelateMotion(IMf, aData,params.baselineWindow);

% normalize by expected poisson noise
IMavg_nans = isnan(IMavg);
IMnoise = max(IMavg,0);
IMnoise = sqrt(convn(IMnoise, abs(DoGfilt), 'same')); % total signal contributing to each measurement
IMnoise = IMnoise + prctile(IMnoise(~IMavg_nans), 50); % add a noise floor to the dimmer pixels
IMnoise(IMavg_nans) = nan;
IMf = IMf./IMnoise;

%temporal matched filter
IMf(nans) = 0;
mem = IMf(:,:,end);
gamma = exp(-1/tau);
for t = size(IMf,3):-1:1
    IMf(:,:,t) = max(0,gamma*mem) + (1-gamma)*IMf(:,:,t);
    mem = IMf(:,:,t);
end
IMf(nans) = nan;

%compute a summary image based on skewness
summary = mean(IMf.^3, 3, 'omitmissing');
summary(imdilate(isnan(summary), ones(5))) = nan; %remove noisy edges
validArea = ~imdilate(nans, ones(5));
P = getTiledPeaks(IMf, validArea, summary);

%plot a figure
plotSummary(summary, IMgamma, P);
end


function plotSummary(summary, IMgamma, P)
red = summary./prctile(summary(~isnan(summary)), 99);
cyan = IMgamma./prctile(IMgamma(~isnan(IMgamma)),99);
figure, imshow(cat(3, red, cyan,cyan));
hold on,
scatter(P.col,P.row,100*(P.val./mean(P.val)).^2, 'm')
end


% 
% function [density peaks sourceR sourceC params] = clusterLocalizations(peaks, params)
% %make an upsampled image
% sz = params.sz;
% upsample = params.upsample;
% ampThresh = min(peaks.val);
% nEventsTotal = length(peaks.t);
% kernelSigma = 1/upsample;
% rGrid = linspace(1,sz(1), upsample*sz(1)+1);
% cGrid = linspace(1,sz(2), upsample*sz(2)+1);
% [rr,cc] = ndgrid(rGrid,cGrid);
% density = zeros(size(rr));
% 
% for eIx = 1:nEventsTotal
%     rMin = peaks.row(eIx)-(8*kernelSigma);
%     rMax = peaks.row(eIx)+(8*kernelSigma);
%     cMin = peaks.col(eIx)-(8*kernelSigma);
%     cMax = peaks.col(eIx)+(8*kernelSigma);
%     rSel = rGrid>rMin & rGrid<rMax;
%     cSel = cGrid>cMin & cGrid<cMax;
% 
%     rLoc = rr(rSel,cSel); cLoc = cc(rSel,cSel);
%     A = mvnpdf([rLoc(:) cLoc(:)], [peaks.row(eIx) peaks.col(eIx)], kernelSigma.*eye(2));
%     A = (peaks.val(eIx)./max(A)).*A;
% 
%     density(rSel,cSel) = density(rSel,cSel)+reshape(A, size(rLoc));
% end
% 
% %figure, scatter(peaks.R, peaks.C);
% figure, imagesc(cGrid,rGrid,density); hold on, scatter(peaks.col, peaks.row, 'r');
% %hold on, scatter(GT.C, GT.R, 'yellow', 'marker', 'x')
% 
% deconvSigma = upsample*sigma./sqrt(ampThresh/2); % sigmaXY/sqrt(amp) is the localization precision; here we assume the max spread
% filtSize = 2*ceil(3*deconvSigma)+1;
% selRows = imdilate(any(density,2), ones(filtSize,1));
% selCols = imdilate(any(density,1), ones(1,filtSize));
% PSF = fspecial('gaussian',filtSize,deconvSigma); %prior on loclaization accuracy
% IMest = density;
% IMest(selRows,selCols) = deconvlucy(density(selRows,selCols),PSF, 50); %should replace with our own algorithm
% 
% BW = imregionalmax(IMest);
% [maxR, maxC] = find(BW);
% [V, sortorder] = sort(IMest(BW), 'descend');
% maxR = maxR(sortorder);
% maxC = maxC(sortorder);
% 
% %compute pairwise distances to cull spurious maxima
% keep = true(1,length(V));
% dMaxima = squareform(pdist([maxR maxC]));
% dMaxima(eye(size(dMaxima), 'logical')) = inf;
% for vIx = 1:length(V)
%     if ~isnan(dMaxima(vIx,vIx))
%         sel = dMaxima(vIx,:)<(upsample);
%         dMaxima(sel,:) = nan;
%         dMaxima(:,sel) = nan;
%         keep(sel) = false;
%     end
% end
% maxR = maxR(keep);
% maxC = maxC(keep);
% V = V(keep);
% k= length(V);
% sourceR = rGrid(maxR);
% sourceC = cGrid(maxC);
% 
% %Perform assignment using k-means-like approach
% weights = ones(1,k);
% assignments = zeros(1, nEventsTotal);
% done = false;
% while ~done
%     zScores = (sourceR - peaks.row').^2 + (sourceC - peaks.col').^2; %squared distance from events to centers
%     likelihoods = normpdf(zScores).*weights;
%     [~, maxInds] = max(likelihoods,[],2);
%     if all(maxInds==assignments)
%         done = true;
%         %remove extra sources
%         for ix = k:-1:1
%             keepSources(ix) = sum(assignments==ix)>=params.minEvents;
%             keepEvents(assignments==ix) = keepSources(ix);
%         end
%         peaks.row = peaks.row(keepEvents);
%         peaks.col = peaks.col(keepEvents);
%         peaks.val = peaks.val(keepEvents);
%         peaks.t = peaks.t(keepEvents);
% 
%         sourceR = sourceR(keepSources);
%         sourceC =sourceC(keepSources);
%         weights = weights(keepSources);
%         k = sum(keepSources);
% 
%         zScores = (sourceR - peaks.row').^2 + (sourceC - peaks.col').^2; %squared distance from events to centers
%         likelihoods = normpdf(zScores).*weights;
%         [~, assignments] = max(likelihoods,[],2);
%         assignProbs = likelihoods./sum(likelihoods,2);
%     else
%         assignments = maxInds;
%         for ii = 1:k
%             weights(ii) = sum(maxInds==ii);
%         end
%     end
% end
% peaks.assignments = assignments;
% peaks.assignProbs = assignProbs;
% 
% %Plot event assignments to sources
% figure('name', 'Event assignments to sources'), imagesc(cGrid,rGrid,density); hold on;
% colors = hsv(k);
% colors = colors(randperm(k),:);
% for sourceIx = 1:k
%     sel = assignProbs(:,sourceIx)>0.5;
%     scatter(sourceC(sourceIx), sourceR(sourceIx), 300, 'marker', 'x', 'markeredgecolor',colors(sourceIx,:), 'linewidth', 2); hold on;
%     scatter(peaks.col(sel), peaks.row(sel),'markeredgecolor', colors(sourceIx,:));
% end
% end

function P = getTiledPeaks(IM, valid, summary)
peaks = valid;
for dim = 1:3
    peaks = peaks & IM > circshift(IM,1,dim) &  IM > circshift(IM,-1,dim);
end

linInds = find(peaks(:));
vvv = IM(linInds);
[rrr,ccc,ttt] = ind2sub(size(peaks), linInds);

sz = size(IM);
tilesize = 64;
tilestartsR = 1:tilesize/2:(sz(1)-tilesize/4);
tileendsR = min(sz(1), tilestartsR+tilesize-1);
tilestartsC = 1:tilesize/2:(sz(2)-tilesize/4);
tileendsC = min(sz(2), tilestartsC+tilesize-1);

summaryVals = summary(sub2ind(size(summary), rrr, ccc));

keep = false(1,length(ttt));
vNorm = zeros(length(ttt),1);

for rix = 1:length(tilestartsR)
    for cix = 1:length(tilestartsC)
        selStats = rrr>=tilestartsR(max(1,rix-1)) & rrr<=tileendsR(min(end,rix+1))  & ccc>=tilestartsC(max(1,cix-1)) & ccc<=tileendsC(min(end,cix+1));

        S = summary(tilestartsR(max(1,rix-1)):tileendsR(min(end,rix+1)), tilestartsC(max(1,cix-1)):tileendsC(min(end,cix+1)));
        Sthresh = prctile(S(:), 50);
        selS = summaryVals>Sthresh;

        vals = vvv(selStats & selS);
        ptile = prctile(vals, [1 50]);
        vals = 3*(vals-ptile(2))./(ptile(2) - ptile(1));
        thresh = ptile(2) + 2*(ptile(2) - ptile(1)); % threshold is 2*[98% confint], corresponding to an SNR of ~6
        
        selTile = rrr>=tilestartsR(rix) & rrr<=tileendsR(rix)  & ccc>=tilestartsC(cix) & ccc<=tileendsC(cix) & vvv>thresh & selS;
        keep(selTile) = true;
        vNorm(selStats & selS) = max(vNorm(selStats & selS), vals);
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

function IM = decorrelateMotion(IM,aData,window)
nanFrames = squeeze(all(isnan(IM),[1 2]));
IMnan = isnan(IM);
IM(IMnan)=0;
sz = size(IM);
nPCs = 20;

P = motionPredictors(aData,window);
P = P(:, ~nanFrames);

%IM = IM-mean(IM,3, 'omitnan'); %to be safe, remove any offset

[U,S,V] = svds(reshape(IM(:,:,~nanFrames), sz(1)*sz(2),[]),nPCs);
correction = zeros(sz(1)*sz(2), size(V,1));
for pc = 1:size(S,1)
    b = mvregress(P',V(:,pc));
    correction = correction + U(:,pc)*S(pc,pc)*b(2:end)'*P(2:end,:); %P(2:end) because we don't correct DC component
end

IM(:,:,~nanFrames) = IM(:,:,~nanFrames) - reshape(correction, sz(1),sz(2),[]);
IM(IMnan) = nan;
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

function hF = visualize_comps(S, sz)
nS = size(S,2);
RGB = rand(3,nS).^2; RGB = RGB./repmat(sum(RGB,1), 3,1);
S_RGB = sqrt([S*RGB(1,:)' S*RGB(2,:)' S*RGB(3,:)']);
S_RGB = 1.5* S_RGB./max(S_RGB(:));
S_RGB = reshape(full(S_RGB), [sz 3]);
hF = figure('Name', 'NMF components'); imshow(S_RGB);
end