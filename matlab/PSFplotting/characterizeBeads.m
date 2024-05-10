function characterizeBeads (im, dXY, dZ, savename)
%by convention, rotate all images so the long axis of the DMD is first
if size(im,1)<size(im,2)
    im = permute(im, [2 1 3]);
end

if nargin<3 %default XY and Z pixel sizes
    warning('Using default pixel sizes!')
    dXY = 0.25;
    dZ = 0.5;
end
sigma = [0.8/dXY 0.8/dXY 2/dZ]; %estimate gaussian parameters for the PSF; X Y Z
upsample = 2;


sz = size(im);
bthresh = 0.2; % a minimum brightness criterion (as fraction of the 90th percentile brightness) for

%3d filter the image to smooth noise; ideal filter is mirror image of PSF
imf = imgaussfilt3(im, sigma*1.5, 'padding', 0);

%find local maxima
mx =    imf>cat(1, imf(2:end,:,:), zeros(1,sz(2),sz(3))) &  imf>cat(1, zeros(1,sz(2),sz(3)), imf(1:end-1,:,:)) & ...
        imf>cat(2, zeros(sz(1),1,sz(3)), imf(:,1:end-1,:)) & imf>cat(2, imf(:,2:end,:), zeros(sz(1),1,sz(3))) & ...
        imf>cat(3, zeros(sz(1),sz(2),1), imf(:,:,1:end-1)) & imf>cat(3, imf(:,:,2:end), zeros(sz(1),sz(2),1));


imf = imf-median(imf(:));
B = imf(mx); %brightnesses of local maxima
mx = mx & (imf > prctile(B,99).*bthresh); %remove local maxima that are not very bright

%TESTING
% for testplane = 7:size(mx,3)-6
% [xc,yc] = find(mx(:,:,testplane));
% figure, imagesc(imf(:,:, testplane)); hold on, scatter(yc,xc, 'r');
% end

% scaling factor for sigma (MH changed on 2024/4/12)
factor = 2.5; % originally 3

mxc = convn(double(mx), ones(2*ceil(3.5*sigma)+1), 'same');
valid = mx & mxc<2; %select pixels for analysis that are sufficiently far from other local maxima
valid([1:ceil(factor*sigma(1)), end-ceil(factor*sigma(1))+1:end],:,:) = false; 
valid(:,[1:ceil(factor*sigma(2)), end-ceil(factor*sigma(2))+1:end],:) = false;
valid(:,:,[1:ceil(factor*sigma(3)), end-ceil(factor*sigma(3))+1:end]) = false;
valid_ixs = find(valid);

[xx,yy,zz] = ind2sub(sz, valid_ixs);

nPatches = 3;
patchIndsY = linspace(1, sz(2), nPatches+1);
patchIndsX = patchIndsY + round((size(im,1)-size(im,2))/2);

%extize the outputs
brightness = nan(nPatches);
resX = nan(nPatches);
resY = nan(nPatches);
resZ = nan(nPatches);

[bx,by,bz] = ndgrid([-ceil(factor*sigma(1)):ceil(factor*sigma(1))], [-ceil(factor*sigma(2)):ceil(factor*sigma(2))], [-ceil(factor*sigma(3)):ceil(factor*sigma(3))]);
[bxq,byq,bzq] = ndgrid(linspace(-ceil(factor*sigma(1)),ceil(factor*sigma(1)), 2*factor*upsample*ceil(sigma(1))+1), linspace(-ceil(factor*sigma(2)), ceil(factor*sigma(2)),  2*factor*upsample*ceil(sigma(2))+1), linspace(-ceil(factor*sigma(3)),ceil(factor*sigma(3)),  2*factor*upsample*ceil(sigma(3))+1));

beadAvAll = nan(2*factor*upsample*ceil(sigma(1))+1,2*factor*upsample*ceil(sigma(2))+1,2*factor*upsample*ceil(sigma(3))+1);
nbeads = length(xx);
bead_ixs = (1:nbeads)';
if nbeads<3
    warning('Too few beads in FOV!')
else
    beadX = nan(1, length(bead_ixs)); beadY = nan(1, length(bead_ixs));  beadZ = nan(1, length(bead_ixs)); beadB = nan(1, length(bead_ixs));
    beadData_q = [];
    for bead_n = length(bead_ixs):-1:1
        beadData = im(xx(bead_ixs(bead_n)) + [-ceil(factor*sigma(1)):ceil(factor*sigma(1))], yy(bead_ixs(bead_n)) + [-ceil(factor*sigma(2)):ceil(factor*sigma(2))], zz(bead_ixs(bead_n)) + [-ceil(factor*sigma(3)):ceil(factor*sigma(3))]);
    
        params = gaussfitn([bx(:) by(:) bz(:)],reshape(beadData, [],1));
        if all(abs(params{3})<1) %if the center is centered
            beadX(bead_n) = params{4}(1,1);
            beadY(bead_n) = params{4}(2,2);
            beadZ(bead_n) = params{4}(3,3);
            beadB(bead_n) = params{2};
        end
        
        %upsample
        beadData_q(:,:,:,bead_n) = interp3(by,bx, bz,beadData, byq+params{3}(2),bxq+params{3}(1),bzq+params{3}(3), 'spline');
    end
    
    %select the smaller beads to ignore big junk
    selb = (beadX.^2+beadY.^2+beadZ.^2) < prctile((beadX.^2+beadY.^2+beadZ.^2), 40); 
    if any(selb)
        beadAvAll = mean(beadData_q(:,:,:,selb),4);
    end
end

%characterize beads in a given region
for patchX = nPatches:-1:1
    for patchY = nPatches:-1:1
        %select beads in this patch
        sel =  xx >= patchIndsX(patchX) & xx < patchIndsX(patchX+1) & yy >= patchIndsY(patchY) & yy < patchIndsY(patchY+1); 
        
        beadAv{patchX,patchY} = nan(2*factor*upsample*ceil(sigma(1))+1,2*factor*upsample*ceil(sigma(2))+1,2*factor*upsample*ceil(sigma(3))+1);

        nbeads = sum(sel);
        bead_ixs = find(sel);
        if nbeads<3
            warning('Too few beads in a region!')
            continue
        end
        
        beadX = nan(1, length(bead_ixs)); beadY = nan(1, length(bead_ixs));  beadZ = nan(1, length(bead_ixs)); beadB = nan(1, length(bead_ixs));
        beadData_q = [];
        for bead_n = length(bead_ixs):-1:1
            beadData = im(xx(bead_ixs(bead_n)) + [-ceil(factor*sigma(1)):ceil(factor*sigma(1))], yy(bead_ixs(bead_n)) + [-ceil(factor*sigma(2)):ceil(factor*sigma(2))], zz(bead_ixs(bead_n)) + [-ceil(factor*sigma(3)):ceil(factor*sigma(3))]);

            params = gaussfitn([bx(:) by(:) bz(:)],reshape(beadData, [],1));
            if all(abs(params{3})<1) %if the center is centered
                beadX(bead_n) = params{4}(1,1);
                beadY(bead_n) = params{4}(2,2);
                beadZ(bead_n) = params{4}(3,3);
                beadB(bead_n) = params{2};
            end
            
            %upsample
            beadData_q(:,:,:,bead_n) = interp3(by,bx, bz,beadData, byq+params{3}(2),bxq+params{3}(1),bzq+params{3}(3), 'spline');
        end
        
        %select the smaller beads to ignore big junk
        selb = (beadX.^2+beadY.^2+beadZ.^2) < prctile((beadX.^2+beadY.^2+beadZ.^2), 40);
        brightness(patchX,patchY) = mean(beadB(selb), 'omitnan');
        resX(patchX,patchY) = mean(beadX(selb), 'omitnan');
        resY(patchX,patchY) = mean(beadY(selb), 'omitnan');
        resZ(patchX,patchY) = mean(beadZ(selb), 'omitnan');
        
        if any(selb)
            beadAv{patchX,patchY} = mean(beadData_q(:,:,:,selb),4);
        end
    end
end

% figure, imagesc(1./resX, [0 max(1./resX(:))]); axis image; axis off; colorbar; title('Vertical Resolution');
% figure, imagesc(1./resY,[0 max(1./resY(:))]); axis image; axis off; colorbar; title('Horizontal Resolution');
% figure, imagesc(resZ, [0 max(1./resZ(:))]); axis image; axis off; colorbar; title('Axial Resolution');
% figure, imagesc(brightness, [0 max(brightness(:))]); axis image; axis off; colorbar; title('Brightness');

hF1 = fwhmPlots(beadAv, upsample./[dXY dXY dZ]);

avImXY = cell2mat(beadAv);
avImXZ = cell2mat(cellfun(@(x) permute(x, [1 3 2]), beadAv, 'UniformOutput',false));
avImYZ = cell2mat(cellfun(@(x) permute(x, [2 3 1]), beadAv, 'UniformOutput',false));
hF2 = resolutionPlots(avImXY(:,:, ceil((end+1)/2))', nPatches, 'XY projections', [dXY*upsample dXY*upsample]);
hF3 = resolutionPlots(avImXZ(:,:, ceil((end+1)/2))', nPatches, 'XZ projections', [dXY*upsample dZ*upsample]);
hF4 = resolutionPlots(avImYZ(:,:, ceil((end+1)/2))', nPatches, 'YZ projections', [dXY*upsample dZ*upsample]);

% display reference stack (red) with selected beads highlighted (blue)
figure;
imBW = (im-min(im(:)))/(max(im(:))-min(im(:)));
highlights = convn(valid, ones(2*ceil(3.5*sigma)+1), 'same');
highlightedStack = 255*cat(4,imBW,zeros(size(imBW)),highlights);
imshow3D(permute(highlightedStack,[2 1 3 4]));
% draw patch grid
hold on;
for idx = 1:length(patchIndsY)
    plot([patchIndsX(1) patchIndsX(end)],[patchIndsY(idx) patchIndsY(idx)],'g')
end
for idx = 1:length(patchIndsX)
    plot([patchIndsX(idx) patchIndsX(idx)],[patchIndsY(1) patchIndsY(end)],'g')
end

if nargin>3
    disp('Saving figures...')
    savefig(hF1, [savename '_FWHMs.fig'], 'compact');
    savefig(hF2, [savename '_XY.fig'], 'compact');
    savefig(hF3, [savename '_XZ.fig'], 'compact');
    savefig(hF4, [savename '_YZ.fig'], 'compact');
    save([savename '_AVG_PSF'],'beadAvAll');
end
end

function hF = fwhmPlots(beadAv, pixPerUm)
hF = figure;
colors = hsv(3);
for patchR = size(beadAv,1):-1:1
    for patchC = size(beadAv,2):-1:1
        ind = sub2ind([size(beadAv,1), size(beadAv,2)], patchR,patchC);
        ax(ind) = subplot(size(beadAv,2), size(beadAv,1), ind);
        for AX = 1:3
            t = shiftdim(beadAv{patchR,patchC}, AX-1);
            t = mean(mean(t(:, ceil((end+1)/2)+[-2:2], ceil((end+1)/2)+[-2:2]), 2),3) ;
            width(AX) = fwhm(t(2:end-1))/pixPerUm(AX);
            plot(((0:length(t)-1) - (length(t)/2))/pixPerUm(AX), t, 'color', colors(AX,:), 'linewidth', 2); hold on;
        end
        title(['FWHMs: ' num2str(width,2) ' um'])
        xlabel('um'); ylabel('brightness')
    end
end
set(ax, 'box', 'off', 'tickdir', 'out')
linkaxes(ax);
end

function width = fwhm(y)
x = 1:length(y);
width = nan;
y = y-min(y);
y = y./max(y);
N = length(y);
lev50 = 0.5;
i = find(y>0.5,1, 'first');
if i ==1
    return
end
interp = (lev50-y(i-1)) / (y(i)-y(i-1));
tlead = x(i-1) + interp;
i = find(y>0.5,1, 'last')+1;
if i <=N
    interp = (lev50-y(i-1)) / (y(i)-y(i-1));
    ttrail = x(i-1) + interp*(x(i)-x(i-1));
    width = ttrail - tlead;
end
end

function fwhm = diffractionLimitXY(NA)
lambda = 1.03;
fwhm = 2 * (0.32*lambda./(sqrt(2).*NA));
fwhm(NA>0.7) = 2 * (0.325*lambda./(sqrt(2).*NA(NA>0.7).^0.91));
end

function fwhm = diffractionLimitZ(NA)
lambda = 1.03;
fwhm = 2 * (0.532*lambda./(sqrt(2))) .* (1./(1.33-sqrt(1.33.^2 - NA.^2)));
end


    function hF = resolutionPlots(tmpIM, nPatches, name, pixPerUm)
        %IMAGES
        hF = figure;
        tmpIM = sqrt(tmpIM); %apply gamma=0.5 correction
        imshow(tmpIM, [0 max(tmpIM(:))], 'initialmag', 800);
        title(name)
        hold on;
        dR = size(tmpIM,2)./nPatches;
        dC = size(tmpIM,1)./nPatches;
        for patchX = 1:nPatches-1
            plot(patchX.*dR.*[1 1], [0 size(tmpIM,1)+1] , 'r', 'linewidth', 2)
        end
        for patchY = 1:nPatches-1
            plot([0 size(tmpIM,2)+1],patchY.*dC.*[1 1],  'r', 'linewidth', 2)
        end
        label = sprintf('%sm', '\mu'); % micrometer
        hScalebarX = scalebar(gca, 'x', 1, label,'Location', 'southeast', 'ConversionFactor', pixPerUm(1));
        hScalebarX.Color = 'w';
        hScalebarY = scalebar(gca, 'y', 1, '', 'Location', 'southeast', 'ConversionFactor', pixPerUm(2));
        hScalebarY.Color = 'w';
        hcb = colorbar; hcb.Label.String = 'sqrt(brightness) (au)';
        % hcb = colorbar; hcb.Label.String = 'brightness (au)';
        hcb.Position = [0.9 0.25 0.02 0.5];
    end
