function [data, desc, meta] = networkScanImageTiffReader(fname, localDr)

[dr, name, ext] = fileparts(fname);

if all(fname(1:2) == 'Z:') || all(fname(1:2) == '\\')
    if ~exist('localDr', 'var'); localDr = 'C:\temp'; end % 'F:\tmp_tiffIO'; end

    dr = localDr;

    if exist(dr,'dir') ~= 7; mkdir(dr); end

    try
        copyfile(fname,dr);
    catch
        disp(['Could not copy file ' fname]);
        data = nan;
        desc = nan;
        meta = nan;
    end
end


A = ScanImageTiffReader(fullfile(dr, [name ext]));
desc = A.descriptions();
meta = A.metadata;

data = single(A.data);

clear('A');

if all(fname(1:2) == 'Z:') || all(fname(1:2) == '\\')
    delete(fullfile(localDr, [name ext]));
end

end