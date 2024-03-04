function networkTiffWriter(mat, fnwrite, pixelscale, localDr)

[dr, name, ext] = fileparts(fnwrite);

writeDr = dr;

    if ~exist('localDr', 'var'); localDr = 'F:\tmp_tiffIO'; end
if all(fnwrite(1:2) == 'Z:') || all(fnwrite(1:2) == '\\')
  
    writeDr = localDr;
    if exist(writeDr,'dir') ~= 7; mkdir(writeDr); end
end

fTIF = Fast_BigTiff_Write(fullfile(writeDr, [name ext]),pixelscale,0);

for frame = 1:size(mat,3)
    fTIF.WriteIMG(single(mat(:,:,frame)));
end

fTIF.close;

if all(fnwrite(1:2) == 'Z:') || all(fnwrite(1:2) == '\\')
    try
        copyfile(fullfile(writeDr, [name ext]),dr);
        delete(fullfile(writeDr, [name ext]));
    catch
        disp(['Could not copy file to dir ' dr])
    end
end

end