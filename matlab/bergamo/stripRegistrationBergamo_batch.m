dr = uigetdir();

dataDrs = dir(fullfile(dr,'scan_*'));
drSize = size(dataDrs);
numScans = drSize(1);
disp(['Registering ' , int2str(numScans) ,  ' scans'])


parfor idx = 1:numScans
    
    workingScan = fullfile(dataDrs(idx).folder,dataDrs(idx).name);
    disp(workingScan)
    try
        stripRegistrationBergamo(3, workingScan);
        disp([dataDrs(idx).name, ' Complete'])
    catch
        disp('Bugged Out')
    end

end
