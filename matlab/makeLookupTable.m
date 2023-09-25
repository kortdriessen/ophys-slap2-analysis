function motionLUT = makeLookupTable(reference, sparseMaskInds, numFastZs, varargin)

[dmdPixelsPerColumn, dmdPixelsPerRow, numPlanes] = size(reference);

if nargin == 3
    yMotRange = -25:25;
    xMotRange = -25:25;
    zMotRange = 1:numPlanes;
else
    yMotRange = varargin{1};
    xMotRange = varargin{2};
    zMotRange = varargin{3};
end


fprintf("Calculating lookup table... ")

tic;
motionLUT = single(zeros(length(yMotRange), length(xMotRange), length(zMotRange),max(sparseMaskInds(:,2)))); % X x Y x Z x no. of superpixels
for n = 1:max(sparseMaskInds(:,2))
    openPixs = sparseMaskInds(sparseMaskInds(:,2) == n, 1);

    refPix = ceil(length(openPixs)/2);

    [r, c, d] = ind2sub([dmdPixelsPerColumn dmdPixelsPerRow numFastZs], openPixs);
    for y = 1:length(yMotRange)
        for x = 1:length(xMotRange)
            for z = 1:length(zMotRange)
                % select expected value of superpixel at (x,y,z)
                % displacement from the reference image
                motionLUT(y,x,z,n) = reference(r(refPix)+yMotRange(y), c(refPix)+xMotRange(x), d(refPix)+zMotRange(z)-1);
            end
        end
    end
end
clear('tmpMask');

fprintf('done - took %f sec\n', toc);