function [ydata, time] = extractIntegrationROIs(smoothwindow, baselinewindow, chIdx)

if nargin<2
    smoothwindow = 20;
    baselinewindow = 5000;
end
if nargin<3
    chIdx = 1;
end

% %load a slap2 datafile
% hDF = slap2.Slap2DataFile();
% parsePlan = hDF.hDataFile.metaData.AcquisitionContainer.ParsePlan  
% %for every ROI, get the timeseries
% [deltaFOverF, dFFerr, tq] = getTimeSeries(hDF, 1, iPixels, 100)

hDataFile = slap2.util.DataFile(); %'acquisition_20230907_175045_DMD2.dat'
dt = hDataFile.metaData.linePeriod_s;
zIdx = 1;
chIdx = 1;
nROIs = numel(hDataFile.metaData.AcquisitionContainer.ROIs.rois);

labels = zeros(800,1280);
for rix = 1:nROIs
    hTrace(rix) = slap2.util.datafile.trace.Trace(hDataFile,zIdx,chIdx);
    
    pixelMask = false(800,1280);
    shapeData = hDataFile.metaData.AcquisitionContainer.ROIs.rois{rix}.shapeData;
    linInds = sub2ind(size(pixelMask), shapeData(:,1), shapeData(:,2));
    labels(linInds) = rix;
    pixelMask(linInds) = true;

    rasterPixels = pixelMask;
    integrationPixels = pixelMask;
    hTrace(rix).setPixelIdxs(rasterPixels,integrationPixels);
    windowWidth_lines = 20;
    expectedWindowWidth_lines = 5000;
    futures(rix) = hTrace(rix).processAsync(windowWidth_lines,expectedWindowWidth_lines);
end
disp('Loading data, please be patient...')
waitfor(futures)
for rix = 1:nROIs
    ydata(:,rix) = futures(rix).fetchOutputs();
end
time = (0:size(ydata,1)-1).*dt;
figure, plot(time,ydata);
xlabel('time (s)')
ylabel('intensity')

figure,
imagesc(labels); axis image;


