function [IMs,meta, ImDescriptions] = readTiff(tiffFn)
hTiff= ScanImageTiffReader(tiffFn);
try
    ImDescriptions = hTiff.descriptions;
    meta = hTiff.metadata();
    IMs = hTiff.data;
catch ME
    most.idioms.safeDeleteObj(hTiff);
    ME.rethrow();
end
try
    delete(hTiff)
catch
end
% most.idioms.safeDeleteObj(hTiff);
end