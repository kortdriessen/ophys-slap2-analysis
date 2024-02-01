dr = uigetdir();

dataDrs = dir(fullfile(dr,'scan_*'));

<<<<<<< HEAD
parfor idx = 1:length(dataDrs)
=======
for idx = 1:length(dataDrs)
>>>>>>> 789268ec4e8623de3b82bf55912630cce65fe440
    disp(['Running directory ' dataDrs(idx).name])
    stripRegistrationBergamo([],fullfile(dataDrs(idx).folder,dataDrs(idx).name,[dataDrs(idx).name '.tif']));
    
    tmp = dir(fullfile(dataDrs(idx).folder,dataDrs(idx).name,dataDrs(idx).name,'*_DOWNSAMPLED*'));
    downsampledFile = tmp.name;
    summarizeBergamo_Peaks(fullfile(dataDrs(idx).folder,dataDrs(idx).name,dataDrs(idx).name),downsampledFile);
end