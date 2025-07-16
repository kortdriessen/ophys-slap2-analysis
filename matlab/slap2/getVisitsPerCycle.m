function [numVisitsRaster, numVisitsBlock, cycleRateHz] = getVisitsPerCycle(meta)

parsePlan = meta.AcquisitionContainer.ParsePlan.acqParsePlan;
numVisitsRaster = zeros(meta.dmdPixelsPerRow,meta.dmdPixelsPerColumn, 'uint16');
numVisitsBlock = zeros(meta.dmdPixelsPerColumn, meta.dmdPixelsPerRow, 'uint16');
nPx = meta.dmdPixelsPerColumn*meta.dmdPixelsPerRow;
for lineIx = 1:numel(parsePlan) 
    pxList = parsePlan(lineIx).superPixelID;
    rasterPx = pxList(pxList<=nPx);
    numVisitsRaster(rasterPx) = numVisitsRaster(rasterPx)+1;
end

numVisitsRaster = numVisitsRaster';
numVisitsBlock = numVisitsBlock';
cycleRateHz = 1/(numel(parsePlan)*meta.linePeriod_s);