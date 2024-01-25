dr = uigetdir();

dataDrs = dir(fullfile(dr,'scan_*'));

for idx = 1:length(dataDrs)
    disp(['Running directory ' dataDrs(idx).name])
    tmp = dir(fullfile(dataDrs(idx).folder,dataDrs(idx).name,'*_DOWNSAMPLED*'));
    downsampledFile = tmp.name;
    spAnalysis = spineAnalysis(fullfile(dataDrs(idx).folder,dataDrs(idx).name,downsampledFile));

    tmp = dir(fullfile(dataDrs(idx).folder,dataDrs(idx).name,'*ROIs.mat'));
    roisFile = tmp(1).name;
    spAnalysis.loadROIsDirect(fullfile(dataDrs(idx).folder,dataDrs(idx).name,downsampledFile,roisFile));

    try
        spAnalysis.analyze;
    catch
        msgbox(['Error analyzing: ' dataDrs(idx).name]);
    end
end