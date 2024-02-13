function [refPlane, corrPlane, skewPlane] = alignMultiChannel(IMs, b,a)
Y = squeeze(makeHighPass(sum(IMs,3))); % data order is XYCTZV
options_rigid = NoRMCorreSetParms('d1',size(Y,1),'d2',size(Y,2),'bin_width',200,'max_shift',20,'us_fac',4,'init_batch',40, 'correct_bidir', true);
refPlane = [];
[~,shifts,~,options, col_shift] = normcorre(Y,options_rigid);
for chIx= size(IMs,3):-1:1
    tmp= apply_shifts(squeeze(IMs(:,:,chIx,:)),shifts,options,0,0,0,col_shift);
   
    refPlane(:,:,chIx) = mean(tmp,3);

    %compute correlation image
    %downsample in space
    nDS = 1;
    for DSiter = 1:nDS
        tmp = tmp(1:2:2*floor(end/2), 1:2:2*floor(end/2),:) + tmp(2:2:2*floor(end/2), 1:2:2*floor(end/2),:) + tmp(1:2:2*floor(end/2), 2:2:2*floor(end/2),:) + tmp(2:2:2*floor(end/2), 2:2:2*floor(end/2),:);
    end
    %highpass in time
    HP = permute(filtfilt(b,a,permute(double(tmp), [3 1 2])), [2 3 1]);
    ss = sum(HP.^2,3, 'omitnan');
    vertC = sum(HP .* circshift(HP, [1 0 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [1 0 0 ]));
    horzC = sum(HP .* circshift(HP, [0 1 0 0]),3, 'omitnan')./sqrt(ss.*circshift(ss, [0 1 0 ]));
    C = mean(cat(3, horzC, circshift(horzC,1,2), vertC, circshift(vertC, 1,1)),3, 'omitnan');
    C(isnan(C)) = 0;
% 
    %compute skewness image
%     IMgamma = mean(tmp,3); %average downsampled image
%     IMgamma = max(0, IMgamma - prctile(IMgamma, 1,'all')); %estimate baseline
%     IMgamma = sqrt(IMgamma);
    sk = skewness(HP,0,3);
    
    corrPlane(:,:,chIx) = C;
    skewPlane(:,:,chIx) = sk;
end
end