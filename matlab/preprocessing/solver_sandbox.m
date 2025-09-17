%working on solver for SLAP2 glutamate recordings
%TO DO
%Split out solver from simulation
%implement parallelization across subproblems
%
%Baseline estimation- baseline seems to increase over iterations, Why??
%
% NaN handling!!

clear all;
rng(0);

%Parameters
sz = [200,200];
num_sources = 100; nSources_GT = num_sources;
tau = 3;
num_pixels = prod(sz);
num_time_points = 1000;
k = [zeros(1, 6*tau) exp(-(0:6*tau)/tau)];
kernel_len = numel(k);
k_flip = fliplr(k);
alpha = exp(-1/tau);

%inference parameters
params.denoiseWindow_samps = 5;
params.baselineWindow_samps = 101;
params.sigma = 7;
params.selRadius = 7;
params.validRadius = 4;
params.k = k;
params.k_flip = k_flip;

%precompute the filter for turning Hs into H
params.Hsigma = 1; %spatial filter for sources, 2D, usually a gaussian
tmp = zeros(4*ceil(params.Hsigma)+1);
tmp(ceil(end/2), ceil(end/2)) = 1;
params.Hfilter = imgaussfilt(tmp, params.Hsigma, 'FilterSize',size(tmp,1));

%precompute the filter for estimating high-frequency baseline changes from
%surrounding pixels
tmp = zeros(4*params.selRadius+1);
tmp(ceil(end/2), ceil(end/2)) = 1;
tmp1= imgaussfilt(tmp, params.selRadius, 'FilterSize',4*params.selRadius+1); tmp1 = tmp1./max(tmp1(:));
tmp2 = imgaussfilt(tmp, params.selRadius/2,  'FilterSize',4*params.selRadius+1); tmp2 = tmp2./max(tmp2(:));
params.Bfilter = max(0, tmp1-tmp2).^2;

%simulate background
Bim = imgaussfilt(rand(sz), 2); Bim = max(0, Bim-prctile(Bim(:),20))+(0.1.*Bim);
B_GT = Bim .* reshape(1+exp(-(1:num_time_points)/num_time_points), [1,1,num_time_points]);
doMotion = true;
if doMotion
    Bmotion = smoothdata(imgaussfilt(randn([2*sz num_time_points]) , 4*params.selRadius+1), 3,"movmean",5); %fast, nonlocal noise
    Bmotion = Bmotion(ceil(sz(1)/2)+(1:sz(1)), ceil(sz(2)/2)+(1:sz(2)), :);
    Bmotion = 1+ 0.05.*Bmotion./std(Bmotion(:));
    B_GT = B_GT.*Bmotion;
    shiftR = rand(1,num_time_points)/4;
    shiftC = rand(1,num_time_points)/4;
    for frameN = 1:num_time_points
        B_GT(:,:,frameN) = imtranslate(B_GT(:,:,frameN), [shiftR(frameN) shiftC(frameN)]);
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
    Hs_GT(:,ix) = tmp(selPix);
    tmp = imgaussfilt(double(tmp), 1.5).*Bim;
    H_GT(:,ix) = tmp(selPix);
end

%simulate spikes
S_GT = (2^5) .* max(0, rand(nSources_GT,num_time_points).^8 - 0.5);
%S_GT = 0.*S_GT; %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

%normalize H to maintain scale
normFac = sqrt(sum(H_GT.^2, 1));
H_GT = H_GT./normFac;
Hs_GT = Hs_GT./normFac;
S_GT = S_GT.*normFac';

B_GT = reshape(B_GT, [],num_time_points); B_GT = B_GT(selPix(:),:);
X_GT = convn(S_GT, k, 'same');
Y_GT = B_GT + H_GT*X_GT;

%simulate observations
F = 10.*(1+10*rand(size(Y_GT,1),1)) .* (1+rand(size(Y_GT))); %F=Freshness; a weighting factor proportional to the number of independent observations averaged into a measurement
Y_obs = poissrnd(F.*Y_GT)./F; %F=Freshness

params.lambda = prctile(mean(Y_obs,2),5);

%break the problem into separable chunks based on selPix
selIdxs = nan(sz);
selIdxs(selPix) = 1:sum(selPix(:));
CC = bwconncomp(selPix);
for problemIx = 1:CC.NumObjects %10
    pxList = CC.PixelIdxList{problemIx};

    selSources = any(sub2ind(sz,sources.R,sources.C)==pxList',2);
    selPix_p = false(sz);
    selPix_p(pxList) = true;
    sources_p.R = sources.R(selSources);
    sources_p.C = sources.C(selSources);
    Y_p = Y_obs(selIdxs(pxList),:); %observations
    F_inv_p = 1./F(selIdxs(pxList),:); % inverse freshness

    %assemble ground truth
    GT.S = S_GT(selSources,:); % spikes,sources x time
    GT.X = X_GT(selSources,:); % S convolved with kernel
    GT.B = B_GT(selIdxs(pxList),:); % background
    GT.H = H_GT(selIdxs(pxList), selSources); % source images; pixels x sources
    GT.Hs = Hs_GT(selIdxs(pxList), selSources); % superres source images; pixels x sources

    %crop, to reduce memory and improve visualization
    rowSupport = [find(any(selPix_p,2),1,'first') find(any(selPix_p,2),1,'last')];
    colSupport = [find(any(selPix_p,1),1,'first') find(any(selPix_p,1),1,'last')];
    sources_p.R = sources_p.R-rowSupport(1)+1;
    sources_p.C = sources_p.C-colSupport(1)+1;
    selPix_p = selPix_p(rowSupport(1):rowSupport(2), colSupport(1):colSupport(2));

    [H_est,S_est, B_est] = extractSources(Y_p, F_inv_p, sources_p, selPix_p, params, GT);
end

keyboard
%evaluate performance


function [Hs_est,S_est, B_est, errFinal] = extractSources(Y_obs, F, sources, selPix, params, GT)
%performs source extraction by alternating gradient descent on the
%footprints and their activities,in a constrained NMF framework

%Yobs is [#selected pixels]x[#timepoints]
%F is [#selected pixels]x[#timepoints], "freshness", proportional to #measurements averaged into each sample
%selpix is a 2D logical map, containing [#selected pixels] true values, that maps Yobs to 2D space
%params

%outputs the superresolution source estimate (Hs), the superresolution
%spikes (S), and the baseline (B)

sz = size(selPix);
num_sources = numel(sources.R);
num_time_points = size(Y_obs,2);
SE = strel('disk', params.validRadius);

%Initial estimates for B, S, H
for six = num_sources:-1:1
    tmp = zeros(sz);
    tmp(sources.R(six), sources.C(six)) = 1;
    tmp1 = imgaussfilt(tmp, params.Hsigma);
    tmp2 = imdilate(logical(tmp), SE);
    Hs_est(:,six) = tmp(selPix); %initial estimate
    H_est(:,six) = tmp1(selPix); %initial estimate
    Hvalid(:,six) = tmp2(selPix); %spatial support for each source in solver
end

%initialize B
denoised = smoothdata(Y_obs,2,"movmean",params.denoiseWindow_samps,'omitmissing');
LP = splitFreq(denoised, params.baselineWindow_samps);
%medRes = median(denoised-LP,2);
B_est = LP;
typicalX = sqrt(mean((denoised(:,1:100)-LP(:,1:100)).^2,'all'))*ones(num_sources,num_time_points);

%initialize S
S_est = typicalX.*rand(num_sources,num_time_points);

%overwrite with GT for testing?
% B_est = GT.B;
% S_est = GT.S;

problemS.lb = -eps*ones(size(S_est));
problemS.ub = inf(size(S_est));

problemH.lb = -eps*ones(size(Hs_est));
problemH.ub = eps*ones(size(Hs_est));
problemH.ub(Hvalid) = inf;

nLoops = 3;
if nargin>5 %ground truth supplied
    B_err = nan(1, nLoops);
    H_err = nan(1, nLoops);
    S_err = nan(1, nLoops);
    B_err(1) = mean(B_est(:)-GT.B(:)).^2;
    H_err(1) = mean(H_est(:)-GT.H(:)).^2;
    S_err(1) = 1-corr(S_est(:),GT.S(:));
    hAx = plotGT([], S_est, Hs_est, B_est, S_err, H_err, B_err, GT, selPix, params);
end

%optimization options
opts = optimoptions('fmincon', ...
    'Algorithm','trust-region-reflective', ...   %'interior-point'
    'SpecifyObjectiveGradient',true,'SubproblemAlgorithm', 'cg', ...
    'Display','none','MaxPCGIter', 50);

doFitS = true;
doFitH = true;
doFitB = true;

%OPTIMIZE
for outerLoop = 1:nLoops
    opts.MaxIterations = 10*outerLoop;

    %SOLVE FOR S
    if doFitS
        opts.TypicalX = typicalX;
        % Objective function handle (returns [f,g, Hinfo])
        objS = @(x) objfun_S_wrapper(x, Y_obs, H_est, B_est, params.k, F, params.lambda);

        % Hessian multiply for fmincon signature (x,y,flag)
        opts.HessianMultiplyFcn = @(hinfo, v, flag) hessmult_S_wrapper(hinfo, Y_obs, H_est, B_est, params.k, F, params.lambda, v); %hessmult_S_wrapper takes (Svec, Z, H, B, k, F, lambda, v)

        % Call fmincon
        [S_est_new, lossS] = fmincon(objS, S_est, [], [], [], [], problemS.lb, problemS.ub, [], opts);
    else
        S_est_new = S_est;
    end
    %update X
    X_est_new = convn(S_est_new, params.k, 'same');

    %SOLVE FOR Hs
    if doFitH
        opts.TypicalX = max(H_est, 0.01);
        % Objective function handle (returns [f,g, Hinfo])
        objHs = @(Hs) objfun_Hs_wrapper(Hs, Y_obs, X_est_new, B_est, params.Hfilter, selPix, F, params.lambda); %(Hs_vec, Z, X, B, kk, selPix, F, lambda, v)
        opts.HessianMultiplyFcn = @(Hinfo, v, flag) hessmult_Hs_wrapper(Hinfo, Y_obs, X_est_new, B_est, params.Hfilter,selPix, F, params.lambda, v);
        % Call fmincon
        [Hs_est_new,lossH] = fmincon(objHs, Hs_est, [], [], [], [], problemH.lb, problemH.ub, [], opts);
    else
        Hs_est_new = Hs_est;
    end

    %update H
    tmp = zeros([sz num_sources]);
    tmp(repmat(selPix,1,1,num_sources)) = Hs_est_new;
    tmp = reshape(convn(tmp,params.Hfilter,'same'), numel(selPix), num_sources);
    H_est_new = tmp(selPix,:);

    %normalize new H and S
    normFac = sqrt(sum(H_est_new.^2, 1));
    H_est_new = H_est_new./normFac;
    Hs_est_new = Hs_est_new./normFac;
    S_est_new = S_est_new.*normFac';
    X_est_new = X_est_new.*normFac';

    %update with new values
    S_est = S_est_new;
    H_est = H_est_new;
    Hs_est = Hs_est_new;

    if outerLoop==nLoops
        break %we stop optimizing here, since updating B is only heuristic
    end

    %SOLVE FOR B
    if doFitB
        HX = H_est_new*X_est_new;
        B_est_new = fitB(Y_obs,HX,B_est, selPix, params);
        B_est = B_est_new;
    end

    if nargin>5 %if ground truth was supplied, make error plots
        B_err(outerLoop+1) = mean(B_est(:)-GT.B(:)).^2;
        H_err(outerLoop+1) = mean(H_est(:)-GT.H(:)).^2;
        S_err(outerLoop+1) = 1-corr(S_est(:),GT.S(:));
        plotGT(hAx, S_est, H_est, B_est, S_err, H_err, B_err, GT, selPix, params);
        errFinal = S_err(outerLoop+1);
    else
        errFinal = nan;
    end
end


end

%checkGradients
% [valid,err] = checkGradients(problemS.objective, problemS.x0);
[valid,err] = checkGradients(problemH.objective, problemH.x0);

function B = fitB(Y,HX, B0, selPix, params)
%fits a baseline (i.e. F0) to the measurements in Y, given the current
%estimate of the activity, HX

%calculate corrections to avoid inflation of spikes
meanRes = mean(Y-HX-B0,2, 'omitmissing');
HXfloor = smoothdata(ordfilt2(HX,1, ones(1,params.baselineWindow_samps)), 2,"movmean",params.baselineWindow_samps,"omitmissing");

HX1 = HX + meanRes- HXfloor; %adjust our fit by reducing the activity component; this encourages baseline to go up
B1 = Y-HX1; %the data to fit

%split into HighPass and LowPass components
[LP, HP] = splitFreq(B1, params.baselineWindow_samps);

%convert to pixel space
HPfull = zeros(size(selPix,1), size(selPix,2), size(Y,2));
LPfull = zeros(size(selPix,1), size(selPix,2), size(Y,2));
HPfull(repmat(selPix, 1, 1, size(Y,2))) = HP;
LPfull(repmat(selPix, 1, 1, size(Y,2))) = LP;

%compute the ratio of the average highpass vs lowpass signal over surrounding pixels
Bmod = convn(HPfull, params.Bfilter, 'same')./ convn(LPfull,params.Bfilter, 'same');

%convert back to matrix form
Bmod = reshape(Bmod, [], size(Bmod,3));
Bmod = Bmod(selPix(:),:);

%regressions
b = sum(HP .* Bmod, 2) ./ sum(Bmod.^2, 2);
B = LP + b.*Bmod;
end

function hAx = plotGT(hAx, S_est, H_est, B_est, S_err, H_err, B_err, GT, selPix, params)
if isempty(hAx) % SET UP PLOTS
    figure;
    hAx(1) = subplot(3,2,1);
    hAx(2) = subplot(3,2,2); plot(B_err); title(hAx(2), 'B_err History');
    hAx(3) = subplot(3,2,3); title('S')
    hAx(4) = subplot(3,2,4); plot(S_err); title('S_err History')
    hAx(5) = subplot(3,2,5); title('H')
    hAx(6) = subplot(3,2,6); plot(H_err); title('H_err History')
end

sz = size(selPix);

cla(hAx(1))
imagesc(B_est - GT.B, 'parent', hAx(1)); title(hAx(1), '{B_{est} - B_{GT}}'); colorbar(hAx(1));

plot(hAx(2), B_err); title(hAx(2), 'B error history');

cla(hAx(3));
plot(conv2(GT.S, params.k, 'same')' + (1:size(GT.S,1)), 'r','parent', hAx(3)),
hold(hAx(3), 'on'),
plot(conv2(S_est, params.k, 'same')' + (1:size(GT.S,1)), 'b', 'parent', hAx(3));
set(hAx(3), 'xlim', [0 1000]);
title(hAx(3), '{X_{est} (blue) , X_{GT} (red)} (first 1000 samples)');

plot(hAx(4), S_err); title(hAx(4), 'S error history');

cla(hAx(5))
Him = nan(sz); Him(selPix) = sum(H_est,2);
imagesc(Him, 'parent', hAx(5)); title(hAx(5), 'H');

plot(hAx(6), H_err); title(hAx(6), 'H error history');

drawnow
end



function [f, gHs, HvHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda, v)
% OBJFUN_HS  loss, gradient, and Hessian-vector product w.r.t Hs
%
% Usage:
%   [f, gHs] = objfun_Hs_new(Hs, Z, X, B, kk, F, lambda);
%   [f, gHs, HvHs] = objfun_Hs_new(Hs, Z, X, B, kk, F, lambda, VHs);
%
% Inputs:
%   Hs     : m-by-p matrix (optimization variable)
%   Z, B, F: m-by-n matrices
%   X      : p-by-n matrix
%   kk     : 2D convolution kernel
%   lambda : scalar (regularizer)
%   VHs    : optional m-by-p perturbation for Hessian-vector product
%
% Outputs:
%   f      : scalar loss
%   gHs    : m-by-p gradient (dL/dHs)
%   HvHs   : m-by-p Hessian-vector product (if VHs provided)
ns = size(Hs,2);
selPix3D = repmat(selPix, 1,1,ns);

% Forward: H from Hs
tmp = zeros([size(selPix) ns]);
tmp(repmat(selPix, 1, 1, ns)) = Hs;
tmp = convn(tmp, kk, 'same');
tmp = reshape(tmp,[numel(selPix), ns]);
H = tmp(selPix(:),:);

T = H * X + B;               % m-by-n
E = Z - T;                   % m-by-n
s = T + lambda;              % m-by-n

% Hessian-vector product (BUT NOT F and G)
if nargin >= 9 && ~isempty(v)
    for rix = size(v,2):-1:1
        % Forward perturbation
        tmp = zeros(size(selPix3D));
        tmp(selPix3D) = v(:,rix);
        tmp = convn(tmp, kk, 'same');
        V = reshape(tmp(selPix3D), size(Hs));

        % Perturbation in T
        dT = V * X;                     % m-by-n

        % coef = 2*(Z + lambda).^2 ./ (F .* (s.^3))
        coef = 2 .* (Z + lambda).^2 ./ ( F .* (s.^3) );   % m-by-n

        % dp = coef .* dT  (elementwise)
        dp = coef .* dT;                % m-by-n

        % Map back to H-space and then to Hs-space
        HvH = dp * X.';                % m-by-p

        %Backprop to Hs
        tmp = zeros(size(selPix3D));
        tmp(selPix3D) = HvH;
        tmp = convn(tmp, rot90(kk,2), 'same');
        HvHs(:,rix) = tmp(selPix3D);
    end
    f = [];
    gHs = [];
else %loss and gradient
    % Compute loss
    % Guard against zero or negative entries in F or s externally if needed.
    denom = F .* s;              % m-by-n
    f = sum( (E.^2) ./ denom , 'all' );

    % Gradient wrt T (elementwise)
    % p = - (2 E s + E.^2) ./ (F .* s.^2)
    p = - (2 .* E .* s + E.^2) ./ ( F .* (s.^2) );  % m-by-n

    % Gradient wrt H
    gH = p * X.';                % m-by-p

    % % Chain rule: gradient wrt Hs
    tmp = zeros([size(selPix) ns]);
    tmp(selPix3D) = gH;
    tmp = convn(tmp, rot90(kk,2), 'same');
    gHs = reshape(tmp(selPix3D), size(Hs));

    HvHs = [];
end
end

function [f, gHs, Hinfo] = objfun_Hs_wrapper(Hs, Z, X, B, kk, selPix, F, lambda)
[f, gHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda);
Hinfo = Hs;
end

function HvHs = hessmult_Hs_wrapper(Hs, Z, X, B, kk, selPix, F, lambda, v)
[~,~,HvHs] = objfun_Hs(Hs, Z, X, B, kk, selPix, F, lambda, v);
end

%%%%%
%%%%%

function [f, gS, HvS] = objfun_S(S, Z, H, B, k, F, lambda, v)
% OBJFUN_S  loss, gradient, Hessian-vector product w.r.t. S
%
% S : p-by-n
% Z,B,F : m-by-n
% H : m-by-p
% k : 1-D kernel (row or column) used in convn(S,k,'same')
% lambda : scalar regularizer, roughly the scale of 1 photon
% VS : optional p-by-n perturbation for Hv (same size as S)
%
% Outputs:
%  f   : scalar loss
%  gS  : p-by-n gradient dL/dS
%  HvS : p-by-n Hessian-vector product (if VS provided), else []
%
% The loss:
%   T = H * X + B,  X = convn(S, k, 'same')
%   f = sum( (Z - T).^2 ./ (F .* (T + lambda)), 'all' )

% Forward
X = convn(S, k, 'same');   % p-by-n
T = H * X + B;             % m-by-n
E = Z - T;                 % m-by-n
s = T + lambda;            % m-by-n

% Hessian-vector product branch
if nargin >= 8 && ~isempty(v)
    for rix = size(v,2):-1:1
        V = reshape(v(:,rix), size(S,2), size(S,1))';

        % Forward map of perturbation through conv: dX = convn(VS, k, 'same')
        dX = convn(V, k, 'same');    % p-by-n

        % Perturbation in T: dT = H * dX
        dT = H * dX;                  % m-by-n

        % Coefficient: coef = 2*(Z + lambda).^2 ./ (F .* s.^3)
        coef = 2 .* (Z + lambda).^2 ./ ( F .* (s.^3) );  % m-by-n

        % dp = coef .* dT
        dp = coef .* dT;              % m-by-n

        % Map back: tmp = H.' * dp  (p-by-n)
        tmp = H.' * dp;              % p-by-n

        % Back through conv adjoint: HvS = convn(tmp, kflip, 'same')
        conv_adj = convn(tmp, flip(k), 'same');  % p-by-n
        HvS(:,rix) = conv_adj(:);
    end
    f = [];
    gS = [];
else
    % Objective value
    denom = F .* s;            % m-by-n
    f = sum( (E.^2) ./ denom, 'all' );

    % Gradient wrt T (elementwise)
    % p_elem = dL/dT = - (2 E s + E.^2) ./ (F .* s.^2)
    p_elem = - (2 .* E .* s + E.^2) ./ ( F .* (s.^2) );  % m-by-n

    % Gradient wrt X: gX = H.' * p_elem   (p-by-n)
    gX = H.' * p_elem;        % p-by-n

    % Gradient wrt S: adjoint of convn => convn(gX, flip(k), 'same')
    gS = convn(gX, flip(k), 'same');  % p-by-n

    HvS = [];
end
end

function [f, gS, Hinfo] = objfun_S_wrapper(S, Z, H, B, k, F, lambda)
[f, gS] = objfun_S(S, Z, H, B, k, F, lambda); %objfun_S takes (S, Z, H, B, k, F, lambda, v)
Hinfo = S; %Hessian depends on S
end

function HvS = hessmult_S_wrapper(S, Z, H, B, k, F, lambda, v)
% Svec : current S flattened
% p = size(H,2);
% n = size(Z,2);
% S = reshape(S, p, n); %perhaps unnecessary?

[~,~,HvS] = objfun_S(S, Z, H, B, k, F, lambda, v);
end


function [LP,HP] = splitFreq(A, nSamps)
nPages= floor(size(A,2)/nSamps);
totSamps = nPages*nSamps;
t = 1:size(A,2);

a = reshape(A(:,1:totSamps), size(A,1), nSamps, nPages);
aDS = smoothdata(mean(a,2),3, 'lowess',15, 'omitmissing');
a2 = a-aDS;
sigma = std(reshape(a2, size(A,1),[]),0,2);
sel = a2>1.5*sigma;
a2(sel) = 0;
a= a2 + aDS;
aDS = smoothdata(mean(a,2),3, 'lowess',15, 'omitmissing');

tDS = (nSamps+1)/2 + (nSamps).*(0:size(aDS,3)-1);

LP = nan(size(A));
for pxIx = 1:size(LP,1)
    LP(pxIx,:) = interp1(tDS,aDS(pxIx,:)',t,'linear','extrap');
end
HP = A-LP;
end