function D = matchedExpFilter(D, tau_frames)
gamma = exp(-1/tau_frames);
selNans = isnan(D);
D(selNans) = 0;
for t = size(D,2)-1:-1:1
    W1 = gamma.*~isnan(D(:,t+1));
    W2 = (1-gamma).*~isnan(D(:,t));
    D(:,t) = (W1.*(max(0,D(:,t+1))) + W2.*D(:,t))./(W1+W2);
end
D(selNans) = nan;
end