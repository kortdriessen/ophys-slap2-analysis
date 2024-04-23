function D = matchedExpFilter(D, tau_frames)
gamma = exp(-1/tau_frames);
selNans = isnan(D);
%D(selNans) = 0;

D(:,end) = gamma.*D(:,end); %shrink the last timepoint, i.e. behave as if the boundary condition is 0

for t = size(D,2)-1:-1:1
    W1 = gamma.*~isnan(D(:,t+1));
    W2 = (1-gamma).*~isnan(D(:,t));
    D(isnan(D(:,t+1)),t+1) = 0;
    D(:,t) = (W1.*(D(:,t+1)) + W2.*D(:,t))./(W1+W2);
end
D(selNans) = nan;
end