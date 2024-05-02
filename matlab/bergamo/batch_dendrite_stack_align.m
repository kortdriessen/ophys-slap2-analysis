dr = uigetdir();
dataDrs = dir(fullfile(dr,'zstack_*')); %select dr
drSize = size(dataDrs);
numScans = drSize(1);
disp(['Aligning ' , int2str(numScans) ,  ' Stacks'])
idx = 1;
parfor idx = 1:numScans
    fn = dataDrs(idx).name;
    dr = dataDrs(idx).folder;
    disp(fullfile(dr, fn))
    try
        alignBergamoStack_dendriteFunc(dr, fn);
        disp([dataDrs(idx).name, ' Complete'])
    catch
        disp('Bugged Out')
    end
end
