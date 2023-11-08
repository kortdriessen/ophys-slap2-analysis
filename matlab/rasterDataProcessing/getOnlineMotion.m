function [xOffset_pix,yOffset_pix,zOffset_um] = getOnlineMotion (hDataFile, lineIds)

linesPerCycle = numel(hDataFile.lineHeaderIdxs);
numCycles = hDataFile.numCycles;
totalNumLines = linesPerCycle * numCycles;

if nargin<2
    lineIds = 1:linesPerCycle:totalNumLines;
end

xOffset_pix = nan(length(lineIds),1);
yOffset_pix = nan(length(lineIds),1);
zOffset_um = nan(length(lineIds),1);


cycleIdxs = ceil(lineIds/linesPerCycle);
resIdxs = lineIds - (cycleIdxs-1)*linesPerCycle;
for idx = 1:length(lineIds)
    lineHeader = hDataFile.getLineHeader(resIdxs(idx),cycleIdxs(idx));
    xOffset_pix(idx) = lineHeader.xOffset_pix;
    yOffset_pix(idx) = lineHeader.yOffset_pix;
    zOffset_um(idx)  = lineHeader.zOffset_um;
end
