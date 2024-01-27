function F0 = svdF0(IM, nPCs, baselineWindow)
%compute F0 for all pixels in IM, ignoring nans, using an SVD-based
%approach. This is good if individual pixels are noisy
%time is the first dimension
if nargin<2
    nPCs = 4;
end
if nargin<3
    baselineWindow = 101;
end
origsz = size(IM);
IM = reshape(IM, origsz(1), []);
nanframes = all(isnan(IM),2);
selPix = mean(IM(~nanframes,:),1)>0.9;

IMsel = IM(~nanframes, selPix);

IMsmooth = smoothdata(IMsel,1, 'movmean', size(IMsel,1)/4, 'omitmissing');
IMsel(isnan(IMsel)) = IMsmooth(isnan(IMsel));
IMsel(isnan(IMsel)) = 0; %anything left over

[UU,~,VV] = svds(IMsel',nPCs);
selU = mean(abs(UU),1)./sqrt(mean(UU.^2,1));
selU = selU>mean(selU); %select only components that are spatially distributed

VV2 = smoothExp(VV, 'movmedian', baselineWindow);

warning("off", 'stats:statrobustfit:IterationLimit')
Y = IM(~nanframes,:);
for px = 1:size(IM,2)
    y = Y(:,px);
    selT = ~isnan(y);
    if sum(selT)>2*sum(selU)
        %use iteratively reweighted least squares to fit
        b = robustfit(VV2(selT,selU), y(selT), 'bisquare', 2.5, true);
        Y(selT, px) = [ones(sum(selT),1) VV2(selT,selU)  ]*b;
    end
end
F0 = nan(size(IM));
F0(~nanframes,:) = Y;
%F0 = computeF0(F0, denoiseWindow,baselineWindow,2);
F0 = reshape(F0, origsz);
end