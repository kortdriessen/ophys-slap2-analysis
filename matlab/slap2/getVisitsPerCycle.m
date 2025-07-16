function [numVisitsRaster, numVisitsBlock, cycleRateHz] = getVisitsPerCycle(meta)

parsePlan = meta.AcquisitionContainer.ParsePlan.acqParsePlan;
numVisitsRaster = zeros(meta.dmdPixelsPerRow,meta.dmdPixelsPerColumn, 'uint16');
numVisitsBlock = zeros(meta.dmdPixelsPerColumn, meta.dmdPixelsPerRow, 'uint16');
nPx = meta.dmdPixelsPerColumn*meta.dmdPixelsPerRow;
for lineIx = 1:numel(parsePlan) 
    pxList = parsePlan(lineIx).superPixelID+1; %+1 to convert to matlab 1-indexing
    rasterPx = pxList(pxList<=nPx);
    numVisitsRaster(rasterPx) = numVisitsRaster(rasterPx)+1;
end

numVisitsRaster = numVisitsRaster'; %transpose to convert to matlab convention
numVisitsBlock = numVisitsBlock'; %transpose to convert to matlab convention
cycleRateHz = 1/(numel(parsePlan)*meta.linePeriod_s);