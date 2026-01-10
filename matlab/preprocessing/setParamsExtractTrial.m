function params = setParamsExtractTrial (params)

%precompute the filter for turning Hs into H
tmp = zeros(4*ceil(params.sigma_px)+1);
tmp(ceil(end/2), ceil(end/2)) = 1;
params.Hfilter= imgaussfilt(tmp, max(0.5, 0.9*params.sigma_px), 'FilterSize',size(tmp,1));

params.validKernel = true(3,3); %for the spatial footprints, where intensity can be put surrounding the center pixels

%precompute the filters for estimating high-frequency baseline changes from
%surrounding pixels
tmp = zeros(4*params.selRadius+1);
tmp(ceil(end/2), ceil(end/2)) = 1;
tmp1= imgaussfilt(tmp, params.selRadius, 'FilterSize',4*params.selRadius+1); %tmp1 = tmp1./max(tmp1(:));
tmp2 = imgaussfilt(tmp, params.selRadius/2,  'FilterSize',4*params.selRadius+1); %tmp2 = tmp2./max(tmp2(:));
params.Bfilter(:,:,1) = max(0, tmp1-tmp2);%.^2;
mask = zeros(size(tmp)); mask(:,ceil(end/2)+1:end) = 1; mask(:, ceil(end/2)) = 0.5;
mask = mask-mean(mask(:));
params.Bfilter(:,:,2) = params.Bfilter(:,:,1).*mask;
params.Bfilter(:,:,3) = params.Bfilter(:,:,1).*mask';

%precompute temporal filter
params.k = [zeros(1, ceil(6*params.tau_full)) exp(-(0:ceil(6*params.tau_full))/params.tau_full)];
params.k = params.k./sum(params.k);

params.baselineWindow_samps =  params.baselineWindow_Glu_s .*params.analyzeHz;
params.denoiseWindow_samps = params.denoiseWindow_s .* params.analyzeHz;

if ~isfield(params, 'lambda') || isempty(params.lambda)
    % switch params.microscope
    %     case 'SLAP2'
    %         params.lambda = 10;
    %     case 'bergamo'
    %         params.lambda = 70;
    %     otherwise
            warning('params.lambda not supplied; extractTrial will estimate it independently for each recording. This may cause trial-to-trial baseline and scaling variations')
            params.lambda = [];
    % end
end

end