function [F0] = computeF0(Fin, denoiseWindow, cutWindow, algoNum)
%Fin:   fluorescence traces, #timepoints is first dimension
%tau:   the timescale, in samples, of transients; algorithms are tested for
%       sparse spikes with time constant of tau samples
%cutWindow:     sustained transients longer than cutwindow tend to be absorbed
%               into F0
if nargin<4
    algoNum = 1;
end

algos = {@algo1, @algo2};
algo = algos{algoNum};
warning('off');
F0 = algo(Fin, denoiseWindow, cutWindow);
warning('on');
end

function [F0] = algo1(Fin, denoiseWindow, hullWindow)
%algorithm 1:
%   use a median filter to remove noise and spikes
%   apply a rolling convex hull on the median filtered data
%   discard regions near nans that are more sensitive to noise, and replace with
%   extrapolated smooth F0
T = size(Fin,1);
hullWindow = min(hullWindow, floor(T/4));
deltaDes = (denoiseWindow/4); % a safe spacing for computing F0, to save compute time
sampleTimes = round(linspace(1, T, ceil(T/deltaDes)+1));
nSampsInHull = ceil(hullWindow/deltaDes);

%median filter to reduce noise
F0 = medfilt1(Fin,denoiseWindow, [],1,'omitnan', 'truncate');

%minimum convex hull-like operation to remove effect of transients
origsz = size(F0);
F0 = reshape(F0, size(F0,1), []);
for cix = 1:size(F0,2)
    if all(isnan(F0(:,cix)))
        continue
    end
      
    for dix= nSampsInHull:-1:1
        xi =sampleTimes(dix:nSampsInHull:end);
        F00(:,dix) = interp1(xi, F0(xi,cix), sampleTimes);
    end
    FF = min(F00,[],2, 'omitmissing'); 

    %convert to a smooth curve, discarding untrustworthy samples
    doubt = sum(~isnan(F00),2)<ceil(nSampsInHull/2);

        F2 = FF; 
        if sum(~doubt)>2
                F2(doubt) = nan;
        end
        F2 = smoothdata(F2,1,'lowess', 2*ceil(nSampsInHull/2)+1, 'omitmissing');
        F0(:,cix) = interp1(sampleTimes,F2,1:T,"pchip");
end
F0 = reshape(F0, origsz);
F0(isnan(Fin)) = nan;
end



function F0 = algo2(Fin, tau, window)
%fits DFF quickly; assumes first dimension is time
    Ffit = Fin;
    for iter = 1:2
        F0 = smoothdata(Ffit, 1,'movmean',window);
        Ffit = min(Ffit, F0);%-thresh);
    end
    F0 = smoothdata(Ffit, 1,'movmean',window);
end

