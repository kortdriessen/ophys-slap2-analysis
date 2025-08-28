function [motion, R] = xcorr2_nans_weighted(frame, freshness, template, shiftsCenter, dShift)
%Kaspar Podgorski 2023
%perform a somewhat-efficient local normalized crosscorrelation for images with
%nans

%template: the template
%frame: the frame to be aligned; this has more NaNs

%shiftsCenter: the center offset around which to perform a local search
%dShift: the maximum shift (scalar, in pixels) to consider on each axis around shiftsCenter

sz= size(template);
dShift = round(dShift(1)); %sanity check

template = circshift(template, shiftsCenter);

fValid = ~isnan(frame); %valid pixels of the fixed image
[Fr0,Fc0] = find(fValid);

%correlation is sum(A.*B)./(sqrt(ssA)*sqrt(ssB)); ssB is constant though
shifts =  -dShift:dShift;
C = nan(length(shifts), length(shifts));
for drix = 1:length(shifts)
    for dcix = 1:length(shifts)
        Tr = Fr0 + shifts(drix);
        Tc = Fc0 + shifts(dcix);
        
        valid = Tr>=1 & Tr<=size(template,1) & Tc>=1 & Tc<=size(template,2);
        Fr=Fr0(valid); Fc=Fc0(valid); Tr=Tr(valid); Tc=Tc(valid);
        
        T = template(sub2ind(sz, Tr,Tc));
        sel = ~isnan(T); T = T(sel);
        indsF = sub2ind(sz, Fr(sel),Fc(sel));
        F = frame(indsF);
        Ff = freshness(indsF);

        sFT = sum(Ff.*(F-mean(F)).*(T-mean(T))); %./sum(Ff);
        sT = mean((T-mean(T)).^2); 
        sF = sum(Ff.*(F-mean(F)).^2); %./sum(Ff);
        C(drix,dcix) = sFT./sqrt(sT.*sF.*sum(Ff));
    end
end

%find maximum of correlation map
[maxval,I] = max(C(:));
[rr,cc] = ind2sub(size(C),I);

R = maxval; %correlation coefficient

if rr>1 && rr<length(shifts) && cc>1 && cc<length(shifts)
    %perform superresolution upsampling
    ratioR = min(1e6,(C(rr,cc) - C(rr-1,cc))/(C(rr,cc) - C(rr+1,cc)));
    dR = (1-ratioR)/(1+ratioR)/2;
    
    ratioC =min(1e6, (C(rr,cc) - C(rr,cc-1))/(C(rr,cc) - C(rr,cc+1)));
    dC = (1-ratioC)/(1+ratioC)/2;

    motion = -shiftsCenter' + [shifts(rr)-dR shifts(cc)-dC];
else %the optimum is at an edge of search range; no superresolution
    motion = -shiftsCenter' + [shifts(rr) shifts(cc)];
end

motion = -motion;

if any(isnan(motion))
    keyboard
end