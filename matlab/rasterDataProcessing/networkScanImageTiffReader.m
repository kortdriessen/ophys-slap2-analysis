function [data, desc, meta] = networkScanImageTiffReader(fname)

[dr, name, ext] = fileparts(fname);

if fname(1:2) == 'Z:'
    dr = 'F:\tmp_tiffReader';

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

delete(fullfile(dr, [name ext]));

end