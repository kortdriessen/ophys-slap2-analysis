dr = uigetdir();

dataDrs = dir(fullfile(dr,'scan_*'));

parfor idx = 2 %1:length(dataDrs)
    disp(['Running directory ' dataDrs(idx).name])
    stripRegistrationBergamo([],fullfile(dataDrs(idx).folder,dataDrs(idx).name,[dataDrs(idx).name '.tif']));
    
    tmp = dir(fullfile(dataDrs(idx).folder,dataDrs(idx).name,[dataDrs(idx).name dataDrs(idx).name(11:end)],'*_DOWNSAMPLED*'));
    downsampledFile = tmp.name;
    summarizeBergamo_Peaks(fullfile(dataDrs(idx).folder,dataDrs(idx).name,[dataDrs(idx).name dataDrs(idx).name(11:end)]),downsampledFile);
end