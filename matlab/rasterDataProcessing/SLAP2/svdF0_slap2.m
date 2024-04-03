function F0 = svdF0_slap2(IM, nPCs, denoiseWindow, baselineWindow)
%compute F0 for all pixels in IM, ignoring nans, using an SVD-based
%approach. This is good if individual pixels are noisy
%time is the first dimension
if nargin<2
    nPCs = 8;
end
if nargin<3
    denoiseWindow = 33;
end
if nargin<4
    baselineWindow = 101;
end
origsz = size(IM);
IM = reshape(IM, origsz(1), []);
nanframes = all(isnan(IM),2);
selPix = mean(~isnan(IM(~nanframes,:)),1)>0.9;

IMsel = IM(:, selPix);

IMsmooth = smoothdata(IMsel,1, 'movmean', size(IMsel,1)/4, 'omitnan');
IMsel(isnan(IMsel)) = IMsmooth(isnan(IMsel));
IMsel(isnan(IMsel)) = 0; %anything left over


IMsel = smoothdata(IMsel, 'movmean', denoiseWindow, 'omitnan');
IMsel = smoothExp(IMsel, 'movmedian', baselineWindow);

IMsel = IMsel-mean(IMsel,1);

[UU,~,VV] = svds(IMsel',nPCs);
spaceFac = mean(abs(UU),1)./sqrt(mean(UU.^2,1));
selU = (spaceFac>0.6) | (spaceFac>mean(spaceFac)); %select only components that are spatially distributed

%VV2 = smoothExp(VV, 'movmedian', baselineWindow);

warning("off", 'stats:statrobustfit:IterationLimit')
Y = IM;
Ysmooth = smoothdata(Y, 1, 'movmean',denoiseWindow, 'omitnan');
for px = 1:size(IM,2)
    y = Y(:,px);
    selT = ~isnan(y);
    if sum(selT)>2*sum(selU) && sum(selT)>5
        %use iteratively reweighted least squares to fit
        b = regress(y(selT), [VV(selT,selU) ones(sum(selT),1)]);
        Y(selT, px) = min(Ysmooth(selT,px), [VV(selT,selU) ones(sum(selT),1)]*b);
    else
        Y(:,px) = nan; %very few timepoints valid
    end
end
% F0 = nan(size(IM));

Y(nanframes,:) = nan;
%F0 = computeF0(F0, denoiseWindow,baselineWindow,2);
F0 = reshape(Y, origsz);
end