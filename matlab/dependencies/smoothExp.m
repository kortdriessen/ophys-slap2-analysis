function sVals = smoothExp(D, mode, window)
%smooth data but do a linear fit at the start to avoid issues with edge
%effects. Useful for exponential bleaching curves

%first dimension is time
origsz = size(D);
D = reshape(D, size(D,1), []);
window = min(window, size(D,1));
t = (1:window)';
p = [t ones(size(t))];

sVals = smoothdata(D, mode, window, 'omitnan');
lVals = nan(length(t), size(D,2));
for ix = 1:size(D,2)    
    inds = find(~isnan(D(t,ix)));
    b = regress(D(inds,ix), [inds ones(size(inds))]);
    lVals(:,ix) =p*b;%csaps(inds,D(inds,ix), 1/(1 + (window.^3)/6), t);
end
weights = t./length(t);
sVals(t,:) = weights.*sVals(t,:) + (1-weights).*lVals;

mVals = smoothdata(D(1:min(end,2*window),:) - sVals(1:min(end,2*window),:), mode, window, 'omitnan');
sVals(1:window,:) = sVals(1:window,:) + mVals(1:window,:);
sVals = reshape(sVals, origsz);
end



