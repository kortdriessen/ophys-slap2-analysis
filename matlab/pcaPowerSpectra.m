function pcaPowerSpectra()

[IM, ~, ~] = networkScanImageTiffReader('\\allen\aind\scratch\ophys\Adrian\iGluSnFR testing\Michael_testingCode\sourceExtraction_tests\719673\01-23-24\scans\scan_00001_20240123_153740\scan_00001_20240123_153740_REGISTERED_RAW.tif');

IM = permute(reshape(IM, size(IM,1), size(IM,2), 2, []), [2 1 3 4]);

[rows, cols, numChannels, frames] = size(IM);

IM = reshape(IM, [], size(IM,3));
IM(isnan(IM)) = 0;

[UU,~,VV] = svds(double(IM),5);
% Assume VV is your TxN matrix where T is the number of time points
% and N is the number of time series.
[T, N] = size(VV); % Get the dimensions of the matrix
% Frequency vector (assuming sampling rate is 1, adjust as necessary)
Fs = T/120; % Sampling frequency
f = (0:T/2-1)*(Fs/T); % Frequency vector
% Pre-allocate space for power spectra
powerSpectra = zeros(length(f), N);
% Prepare a figure outside the loop for all plots
figure;
hold on; % This allows all plots to be drawn on the same figure
% Loop through each time series
for i = 1:N
% Compute the Fourier Transform of the time series
Y = fft(VV(:,i));
% Compute the two-sided spectrum P2. Then compute the single-sided spectrum P1 based on P2 and the even-valued signal length L.
P2 = abs(Y/T);
P1 = P2(1:T/2);
P1(2:end-1) = 2*P1(2:end-1);
% Normalize the power spectrum
P1_normalized = P1 / max(P1);
% Store the normalized power spectrum in the matrix
powerSpectra(:,i) = P1_normalized;
% Plot the normalized power spectrum
plot(f, P1_normalized, 'DisplayName', ['Series ' num2str(i)]);
end
hold off; % No more plots to add
title('Normalized Power Spectra of Time Series');
xlabel('Frequency (Hz)');
ylabel('Normalized Power');
legend show; % Display a legend to identify each time

figure; clf; hold on; for i = 1:size(VV,2)
plot((1:T)/Fs,(VV(:,i)+(i-1)*.05)/.05 + 1);
end
xlabel('time (s)')
hold off;

[b,a] = butter(4,[20.5 22]/(Fs/2), 'bandpass');
wobbleComp = filtfilt(b,a, VV(:,5));

figure; clf; hold on;
plot((1:T)/Fs, VV(:,5));
plot((1:T)/Fs, VV(:,5) - wobbleComp .* (wobbleComp \ VV(:,5)));
hold off;

figure(30); clf; hold on; for i = 1:size(exptSummary.dFFls{1},1)
plot((1:T)/Fs,(exptSummary.dFFls{1}(i,:,1)+(i-1)*20)/20 + 1, 'black');
end
xlabel('time (s)')
hold off;

nanIdx = isnan(exptSummary.dFFls{1}(:,:,1));
tmpTraces = exptSummary.dFFls{1}(:,:,1);
tmpTraces(nanIdx) = 0;
fixedTraces = tmpTraces - (wobbleComp(1:2:end) \ tmpTraces')' * wobbleComp(1:2:end)';
fixedTraces(nanIdx) = nan;

figure(30); hold on; 
for i = 1:size(fixedTraces,1)
plot((1:2:T)/Fs,(fixedTraces(i,:)+(i-1)*20)/20 + 1, 'red');
end
xlabel('time (s)')
hold off;

fnstem = '\\allen\aind\scratch\ophys\Adrian\iGluSnFR testing\Michael_testingCode\sourceExtraction_tests\719673\01-23-24\scans\scan_00001_20240123_153740\scan_00001_20240123_153740_DENOISED';
pixelscale = 4e4; %PIXEL SIZE IN DOTS PER CM

%save an original-time-resolution recording
fnwrite = [fnstem '_REGISTERED_RAW.tif'];
% fTIF = Fast_BigTiff_Write(fnwrite,pixelscale,0);
tiffSave = single(zeros([cols rows frames * numChannels]));
for frame = 1:frames
    tiffSave(:,:,(frame-1)*numChannels+1) = single(reshape(IM1(:,frame),rows,cols)');
    tiffSave(:,:,(frame-1)*numChannels+2) = single(reshape(IM2(:,frame),rows,cols)');
end
networkTiffWriter(single(tiffSave), fnwrite, pixelscale);
clear('tiffSave')

dsFac = 2;
 %save a downsampled aligned recording
fnwrite = [fnstem '_REGISTERED_DOWNSAMPLED-' int2str(dsFac) 'x.tif'];

nDSframes = frames / dsFac;

tiffSave = single(zeros([cols rows nDSframes*numChannels]));

for DSframe = 1:nDSframes
    tiffSave(:,:,(DSframe-1)*numChannels+1) = single(reshape(sum(IM1(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
    tiffSave(:,:,(DSframe-1)*numChannels+2) = single(reshape(sum(IM2(:,((DSframe-1)*dsFac+1):(DSframe*dsFac)),2),rows,cols)');
end

networkTiffWriter(tiffSave, fnwrite, pixelscale);
clear('tiffSave');

end