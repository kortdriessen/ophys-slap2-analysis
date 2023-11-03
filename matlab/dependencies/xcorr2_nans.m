function [motion, R] = xcorr2_nans(moving, fixed, shiftsCenter, dShift)
%Kaspar Podgorski 2023
%perform a somewhat-efficient local normalized crosscorrelation for images with
%nans

%moving: moving image
%fixed: fixed image
%shiftsCenter: the center offset around which to perform a local search
%dShift: the maximum shift (scalar, in pixels) to consider on each axis around shiftsCenter

dShift = round(dShift(1)); %sanity check

SE = strel(ones(2*dShift+1)); %we will erode with this structuring element to select the set of pixels that are always non-nan within the search space

fValid = imerode(~isnan(fixed) & circshift(~isnan(moving), -shiftsCenter), SE); %valid pixels of the fixed image
fValid(1:dShift,:) = false; fValid(end-dShift+1:end,:) = false; %remove edges
fValid(:,1:dShift) = false; fValid(:,end-dShift+1:end) = false; %remove edges
mValid = circshift(fValid, shiftsCenter); %valid pixels of the moving image

F = fixed(fValid); %fixed data;
ssF = sqrt(sum(F.^2));

%correlation is sum(A.*B)./(sqrt(ssA)*sqrt(ssB)); ssB is constant though
shifts =  -dShift:dShift;
C = nan(length(shifts), length(shifts));
for drix = 1:length(shifts)
    for dcix = 1:length(shifts)
        M = moving(circshift(mValid, [shifts(drix) shifts(dcix)]));
        ssM = sum(M.^2);
        C(drix,dcix) = sum(F .* M)./sqrt(ssM);
    end
end

%find maximum of correlation map
[maxval,I] = max(C(:));
[rr,cc] = ind2sub(size(C),I);

R = maxval./ssF; %correlation coefficient

if rr>1 && rr<length(shifts) && cc>1 && cc<length(shifts)
    %perform superresolution upsampling
    ratioR = min(1e6,(C(rr,cc) - C(rr-1,cc))/(C(rr,cc) - C(rr+1,cc)));
    dR = (1-ratioR)/(1+ratioR)/2;
    
    ratioC =min(1e6, (C(rr,cc) - C(rr,cc-1))/(C(rr,cc) - C(rr,cc+1)));
    dC = (1-ratioC)/(1+ratioC)/2;

    motion = shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
else %the optimum is at an edge of search range; no superresolution
    motion = shiftsCenter' + [shifts(rr) shifts(cc)];
end

if any(isnan(motion))
    keyboard
end