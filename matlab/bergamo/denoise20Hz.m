function [fnRaw, fnRegDS] = denoise20Hz(dr, fns, ds_time)
if nargin<3
    ds_time =1;  
end
dsFac = 2.^ds_time;
% [fns, dr] = uigetfile('*REGISTERED_RAW*.tif', 'multiselect', 'on');

[IM, ~, ~] = networkScanImageTiffReader(fullfile(dr, fns));

ind =strfind(fns, '_REGISTERED');
fnRawStem = fullfile(dr,fns(1:ind));
fnstem = fullfile(dr,[fns(1:ind) 'DENOISED']);
pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

load(fullfile([fnRawStem 'ALIGNMENTDATA.mat']), 'aData');
numChannels = aData.numChannels;
save(fullfile([fnstem '_ALIGNMENTDATA.mat']), 'aData'); %simply copy alignment data over to denoised file

IM = permute(reshape(IM, size(IM,1), size(IM,2), numChannels, []), [2 1 3 4]);

[rows, cols, numChannels, frames] = size(IM);

IM1 = squeeze(IM(:,:,1,:));
if numChannels == 2
    IM2 = squeeze(IM(:,:,2,:));
end

IM1 = reshape(IM1, [], size(IM1,3));
IM1nans = isnan(IM1);
IM1(IM1nans) = 0;

if numChannels == 2
    IM2 = reshape(IM2, [], size(IM2,3));
    IM2nans = isnan(IM2);
    IM2(IM2nans) = 0;
end

Fs = frames / 120;

[b,a] = butter(4,[20 22]/(Fs/2), 'bandpass');

IM1filt = filtfilt(b,a,double(IM1)')';

% F0 = computeF0(IM1',35,ceil(2/aData.frametime))';

[UU1filt,~,VV1filt] = svds(double(IM1filt),5);
selU1 = corr(abs(UU1filt),mean(IM1,2));
selU1 = selU1 > mean(selU1);
correction1 = (UU1filt(:,selU1) \ IM1filt); % (IM1-F0));
fixedIM1 = IM1 - UU1filt(:,selU1) * correction1;
fixedIM1(IM1nans) = nan;

if numChannels == 2
    % F0 = computeF0(IM2',35,ceil(2/aData.frametime))';
    IM2filt = filtfilt(b,a,double(IM2)')';
    [UU2filt,~,VV2filt] = svds(double(IM2filt),5);
    selU2 = corr(abs(UU2filt),mean(IM2,2));
    selU2 = selU2 > mean(selU2);
    correction2 = (UU2filt(:,selU2) \ IM2filt); %(IM2-F0));
    fixedIM2 = IM2 - UU2filt(:,selU2) * correction2;
    fixedIM2(IM2nans) = nan;
end

%save an original-time-resolution recording
fnwrite = [fnstem '_REGISTERED_RAW.tif'];
% fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
tiffSave = single(zeros([cols rows numChannels frames]));
for frame = 1:frames
    tiffSave(:,:,1,frame) = single(reshape(fixedIM1(:,frame),rows,cols)');
    if numChannels == 2
        tiffSave(:,:,2,frame) = single(reshape(fixedIM2(:,frame),rows,cols)');
    end
end
networkTiffWriter(reshape(single(tiffSave),cols,rows,[]), fnwrite, pixelscale);
clear('tiffSave')
[~, fnRaw, ext] = fileparts(fnwrite); fnRaw = [fnRaw ext]; 

 %save a downsampled aligned recording
fnwrite = [fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];

nDSframes = frames / dsFac;

tiffSave = single(zeros([cols rows numChannels nDSframes]));

for DSframe = 1:nDSframes
    tiffSave(:,:,1,DSframe) = single(reshape(sum(fixedIM1(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
    if numChannels == 2
        tiffSave(:,:,2,DSframe) = single(reshape(sum(fixedIM2(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
    end
end

networkTiffWriter(reshape(single(tiffSave),cols, rows, []), fnwrite, pixelscale);
clear('tiffSave');
[~, fnRegDS, ext] = fileparts(fnwrite); fnRegDS = [fnRegDS ext]; 


% figure; clf; hold on; for i = 1:size(VV,2)
% plot((1:T)/Fs,(VV(:,i)+(i-1)*.05)/.05 + 1);
% end
% xlabel('time (s)')
% hold off;
% 
% [b,a] = butter(4,[20.5 22]/(Fs/2), 'bandpass');
% wobbleComp = filtfilt(b,a, VV(:,5));

% figure(30); clf; ax0 = subplot(1,1,1); hold on;
% for i = 1:size(traces,2)
% plot((1:frames)/Fs,(traces(:,i)+(i-1)*20)/20 + 1, 'red');
% end
% xlabel('time (s)')
% hold off;
% 
% figure(31); clf;  ax1 = subplot(1,1,1); hold on;
% for i = 1:size(traces_original,2)
% plot((1:frames)/Fs,(traces_original(:,i)+(i-1)*20)/20 + 1, 'blue');
% end
% xlabel('time (s)')
% hold off;
% 
% linkaxes([ax0 ax1], 'x')

% fnstem = '\\allen\aind\scratch\ophys\Adrian\iGluSnFR testing\Michael_testingCode\sourceExtraction_tests\719673\01-23-24\scans\scan_00001_20240123_153740\scan_00001_20240123_153740_DENOISED';
% pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM
% 
% %save an original-time-resolution recording
% fnwrite = [fnstem '_REGISTERED_RAW.tif'];
% % fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
% tiffSave = single(zeros([cols rows frames * numChannels]));
% for frame = 1:frames
%     tiffSave(:,:,(frame-1)*numChannels+1) = single(reshape(IM1(:,frame),rows,cols)');
%     tiffSave(:,:,(frame-1)*numChannels+2) = single(reshape(IM2(:,frame),rows,cols)');
% end
% networkTiffWriter(single(tiffSave), fnwrite, pixelscale);
% clear('tiffSave')
% 
% dsFac = 2;
%  %save a downsampled aligned recording
% fnwrite = [fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];
% 
% nDSframes = frames / dsFac;
% 
% tiffSave = single(zeros([cols rows nDSframes*numChannels]));
% 
% for DSframe = 1:nDSframes
%     tiffSave(:,:,(DSframe-1)*numChannels+1) = single(reshape(sum(IM1(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
%     tiffSave(:,:,(DSframe-1)*numChannels+2) = single(reshape(sum(IM2(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
% end
% 
% networkTiffWriter(tiffSave, fnwrite, pixelscale);
% clear('tiffSave');

end
