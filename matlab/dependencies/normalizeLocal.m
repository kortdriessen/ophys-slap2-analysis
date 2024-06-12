function IM = normalizeLocal (IM, dXY, weight, timeVary)
    %normalizes an image to its local and global max in a dXY radius
    globalMax = prctile(IM(:), 99.9);
    nans = isnan(IM);
    IM(nans) = 0;
    if nargin<4 || timeVary
        localMax = imdilate(IM, strel('disk', dXY));
    else
        localMax = imdilate(max(IM,[],3), strel('disk', dXY));
    end
    localMax = imgaussfilt(localMax, dXY, 'padding', 'symmetric');
    BG = prctile(IM(:), 33);
    IM = IM./(weight.*localMax + (1-weight).*globalMax + BG);
    IM(nans) = nan;
end
