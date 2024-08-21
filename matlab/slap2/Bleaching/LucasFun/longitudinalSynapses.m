testFile = 'D:\iGluSnFR_testing_project\Spine_Analysis\data\activity_Analysis\710234\12-18-23\scans\scan_00003_20231218_182433\scan_00003_20231218_182433_REGISTERED_DOWNSAMPLED-8x.tif';
testTiff_info = imfinfo(testFile);
nSlices = numel(testTiff_info);
imageData = NaN(testTiff_info(1).Height, testTiff_info(1).Width, nSlices, 'like', imread(testFile));
for k = 1:nSlices
    imageData(:,:,k) = imread(testFile,k);
end
%% After loading test image
%Step 1, reshape tiff stack from 3D to 2D
[h, w, s] = size(imageData);
reshapedStack = reshape(imageData, h*w, s);
% imshow(reshapedStack)
%throwing away movement artifacts outside of abnormal range
nanMask = ~isnan(reshapedStack);
figure, imshow(nanMask);
[~, firstNonNaNRow] = max(nanMask, [], 1);
moveDist = histogram(firstNonNaNRow);
figure, plot(firstNonNaNRow)
[~,idx] = sort(moveDist.Values, 'descend');
stableRange_top = max(moveDist.BinEdges(idx((1:3))));
stableRange_bot = min(moveDist.BinEdges(idx((1:3))))-50;

firstNonNaNRow(firstNonNaNRow > stableRange_top) = NaN;
firstNonNaNRow(firstNonNaNRow < stableRange_bot) = NaN;
stabilityTrace = reshapedStack(:,~isnan(firstNonNaNRow)); 
% stabilityTrace(stableIDXS) = firstNonNaNRow(stableIDXS);
figure, plot(firstNonNaNRow)
figure, imshow(stabilityTrace)
returnedSanityCheck = reshape(reshapedStack, [h,w,s]);
% figure, bar(moveDist)
display('Done')
% Now lets run clustering alg on the reshaped stack

features_omitted = create_img_features(stabilityTrace);
features_raw = create_img_features(reshapedStack);
disp('Done')
%% fitting curve to average recording to characterize noise and bleaching
PCA_the_means(stabilityTrace)
disp('done')

function featureArray = create_img_features(tiffMov)
    avgPix = mean(tiffMov, 2, 'omitnan');
    stdPix = std(tiffMov, 0,2, 'omitnan');
    regularizedMov = (tiffMov - avgPix) ./ stdPix;
    feature1 = max(tiffMov,[],2);
    feature2 = avgPix;
    feature3 = var(regularizedMov, 0,2, 'omitnan');
    feature4 = std(regularizedMov, 0,2, 'omitnan');
    featureArray =[feature1.'; feature2.'; feature3.'; feature4.'];
end

function incdices = find_closest_centroids(X, centroids)
    K = size(centroids,1);
    incdices = zeros(size(X,1));
    for mi = 1:length(size(X,1))
        distVec = zeros(K);
        for c = 1:length(K)
            distVec(c) = sum((X(mi) - centroids(c))^2);
        end
        incdices(mi) = min(distVec);
    end
end

function PCA_the_means(stabilityTrace)
    %Step 1: Take mean of zero axis for stability stack
    backgroundAvg = mean(stabilityTrace, 'omitnan');
    %Step 2: Standardize the data
    z = (backgroundAvg - mean(backgroundAvg, 'omitnan')) / std(backgroundAvg, 'omitnan');
    %Step 3: Compute Covariance between xx, xy, yx, and yy  || this will
    %give us at the very least a gauge of bleaching of the scan (itll be
    %different for different variants)
    covMat = cov(z);
    

end