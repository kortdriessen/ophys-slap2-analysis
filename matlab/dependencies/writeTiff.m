function writeTiff(filename,data,isTransposed,imageDescriptions)
    if iscell(data)
        nFrames = numel(data);
    elseif isnumeric(data)
        nFrames = size(data,3);
    else
        error('A needs to be cell array or numeric array')
    end

    if nargin<3 || isempty(isTransposed)
        isTransposed = false;
    end

    if nargin<4 || isempty(imageDescriptions)
        imageDescriptions = repmat({''},nFrames,1);
    end

    assert(iscellstr(imageDescriptions));
    validateattributes(imageDescriptions,{'cell'},{'numel',nFrames});
    
    hTiff = Tiff(filename,'w8');
    
    try
        for idx = 1:nFrames

            if idx>1
                hTiff.writeDirectory(); % do not write Directory on first frame, otherwise Tiff file will be invalid
            end

            if iscell(data)
                im = data{idx};
            else
                im = data(:,:,idx);
            end

            if isTransposed
                im = im';
            end

            [bitsPerSample,sampleFormat] = getDataFormat(im);

            samplesPerPixel = 1;

            hTiff.setTag('Photometric',Tiff.Photometric.MinIsBlack);
            hTiff.setTag('Compression',Tiff.Compression.None);
            hTiff.setTag('BitsPerSample',bitsPerSample);
            hTiff.setTag('SamplesPerPixel',samplesPerPixel);
            hTiff.setTag('SampleFormat',sampleFormat);
            hTiff.setTag('ImageLength',size(im,1));
            hTiff.setTag('ImageWidth',size(im,2));         
            hTiff.setTag('PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);

            bytesPerPixel = samplesPerPixel*bitsPerSample/8;
            bytesPerRow = bytesPerPixel*size(im,2);
            stripTargetSizeByte = 8000;
            rowsPerStrip = max(1,floor(stripTargetSizeByte/bytesPerRow));
            hTiff.setTag('RowsPerStrip',rowsPerStrip);

            if ~isempty(imageDescriptions{idx})
                hTiff.setTag('ImageDescription',imageDescriptions{idx});
            end
            
            hTiff.write(im);
        end
    catch ME
        hTiff.delete();
        ME.rethrow()
    end
    
    hTiff.delete();
end

function [bitsPerSample,sampleFormat] = getDataFormat(A)
    switch class(A)
        case 'logical'
            bitsPerSample = 8;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'uint8'
            bitsPerSample = 8;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'int8'
            bitsPerSample = 8;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'uint16'
            bitsPerSample = 16;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'int16'
            bitsPerSample = 16;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'uint32'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'int32'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'uint64'
            bitsPerSample = 64;
            sampleFormat = Tiff.SampleFormat.UInt;
        case 'int64'
            bitsPerSample = 64;
            sampleFormat = Tiff.SampleFormat.Int;
        case 'single'
            bitsPerSample = 32;
            sampleFormat = Tiff.SampleFormat.IEEEFP;
        case 'double'
            bitsPerSample = 64;
            sampleFormat = Tiff.SampleFormat.IEEEFP;
        otherwise
            error('Unsupported data format: %s',class(A));
    end
end