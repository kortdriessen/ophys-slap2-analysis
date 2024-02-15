function D = matchedExpFilter(D, tau_frames)
gamma = exp(-1/tau_frames);
selNans = isnan(D);
D(selNans) = 0;
for t = size(D,2)-1:-1:1
    D(:,t) = max(0,gamma*D(:,t+1)) + (1-gamma)*D(:,t);
end
D(selNans) = nan;
end