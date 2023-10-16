
function aData = loadAlignmentData(datFilename)
    aFilename = [datFilename(1:end-4) '_ALIGNMENTDATA.mat'];
    load(aFilename, 'aData');
end