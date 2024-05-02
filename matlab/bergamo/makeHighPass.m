function IMHP = makeHighPass(IM1)
IMHP = IM1; %sum the channels to speed alignment
IMHP(IMHP<0) = 0;
IMHP(isnan(IMHP)) = 0;
%IMHP = imgaussfilt(IMHP, [0.5 0.5])-imgaussfilt(IMHP, [4 4]);
end