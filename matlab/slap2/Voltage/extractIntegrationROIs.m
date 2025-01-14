function extractIntegrationROIs(smoothwindow, baselinewindow, chIdx)

if nargin<2
    smoothwindow = 20;
    baselinewindow = 5000;
end
if nargin<3
    chIdx = 1;
end

hDataFile = slap2.util.DataFile();
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
[folder, fileName, ~] = fileparts(hDataFile.filename);
traces_fname = [fileName, '_TRACES', '.mat'];
traces_fpath = fullfile(folder, traces_fname);
ROIs_fname = [fileName, '_ROIs', '.tiff'];
ROIs_fpath = fullfile(folder, ROIs_fname);
save(traces_fpath, 'time', 'ydata');
imwrite(uint8(labels), ROIs_fpath, 'tiff');