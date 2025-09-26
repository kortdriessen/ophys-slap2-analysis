function varargout = plotDifferentColoredROIS(imgArray, varargin)

% varargout = plotDifferentColoredROIS(imgArray, varargin)
% 
% Plots the ROIs (a 3D array) in different colors to visualize their
% spatial locations. 
%
% From Charles Lab CoDybase-MATLAB
% 2023 - Adam Charles
% 2025 - updated by Michael Xie

p = inputParser;
p.addParameter('figNo',       1);
p.addParameter('normOpt',     'single');
p.addParameter('filterSize',  [0,Inf]);
p.addParameter('aspectRatio', []);
p.addParameter('labels',      true);
p.addParameter('structuralImage',      []);
parse(p,varargin{:});
p = p.Results;

imgArray(isnan(imgArray)) = 0;
p.structuralImage(isnan(p.structuralImage)) = min(p.structuralImage,[],"all","omitmissing");
p.structuralImage = min(1, p.structuralImage / prctile(p.structuralImage(:),99.95));

numPix = sum(sum(bsxfun(@gt, imgArray, 0.1*max(max(imgArray,[],1),[],2)  ),1),2);

imgArray = imgArray(:,:,(numPix>=p.filterSize(1))&(numPix<=p.filterSize(2)));

numImgs = size(imgArray,3);
allColors = distinguishable_colors(numImgs,'k');

if ~isempty(p.structuralImage)
    imgArray = cat(3,imgArray,p.structuralImage);
    numImgs = numImgs+1;
    allColors = cat(1,allColors,[1 1 1]);
end

imgFull = zeros(size(imgArray,1),size(imgArray,2),3);
for ll = 1:numImgs
    TMPimg = imgArray(:,:,ll).*(imgArray(:,:,ll)>0);
    if ll < numImgs
        switch p.normOpt
            case 'all'; TMPimg  = TMPimg./prctile(imgArray(:),99.999);
            otherwise;  TMPimg  = TMPimg./max(TMPimg(:));
        end
    end
    imgFull = imgFull + bsxfun(@times, TMPimg,...
                                       reshape(allColors(ll,:),[1,1,3]));
end

if nargout >0; varargout{1} = imgFull;
else
    figure(p.figNo)
    imagesc(imgFull,'AlphaData',1)
    if isempty(p.aspectRatio); axis image;
    else;                      pbaspect([p.aspectRatio,1]);
    end

    axis off

    if p.labels
        for ll = 1:numImgs
            if ll < numImgs || isempty(p.structuralImage)
                loc = find(imgArray(:,:,ll) == max(reshape(imgArray(:,:,ll),1,[])),1);
                figure(p.figNo); hold on;
                [yloc,xloc] = ind2sub(size(imgArray,1:2),loc);
                text(xloc,yloc,num2str(ll),'FontWeight','bold','FontSize',15,'Color','white','HorizontalAlignment','center','VerticalAlignment','middle');
                text(xloc,yloc,num2str(ll),'FontSize',14,'Color',allColors(ll,:),'HorizontalAlignment','center','VerticalAlignment','middle');
                hold off;
            end
        end
    end
end



end
