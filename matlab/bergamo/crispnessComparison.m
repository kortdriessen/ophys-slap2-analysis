drs = {'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\750098\2024-09-24\scans\scan_00001_20240924_110500',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\750098\2024-09-24\scans\scan_00002_20240924_111556',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\750098\2024-09-24\scans\scan_00006_20240924_120622',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\758796\2024-10-01\scans\scan_00002_20241001_100847',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\758796\2024-10-01\scans\scan_00006_20241001_105830',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\758796\2024-10-01\scans\scan_00008_20241001_112221',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\758796\2024-10-01\scans\scan_00010_20241001_114510',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00001_20240801_110557',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00003_20240801_113412',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00006_20240801_121111'};

% drs = {'/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/750098_01',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/750098_02',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/750098_06',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/758796_02',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/758796_06',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/758796_08',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/758796_10',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/743713_01',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/743713_03',
%     '/root/capsule/data/iGluSnFR4f_dendrite_passiveDriftingGratings_Registered/caiman/743713_06'};

ngs_reg = nan(size(drs));
ngs_orig = nan(size(drs));

for idx = 1:length(drs)
    disp(['Running file ' num2str(idx)])

    basename = drs{idx}(regexp(drs{idx}, 'scan_[0-9]+'):end);
    
    mov_reg_fn = fullfile(drs{idx},[basename '_REGISTERED_RAW.tif']);
    [mov_reg, ~, ~] = networkScanImageTiffReader(mov_reg_fn);
    
    mov_orig_fn = fullfile(drs{idx},[basename '.tif']);
    [mov_orig, ~, ~] = networkScanImageTiffReader(mov_orig_fn);
    
    meanIM_reg = mean(mov_reg(:,:,2:2:end),3,'omitnan');
    meanIM_orig = mean(mov_orig(:,5:end,2:2:end),3,'omitnan');
    
    fullFrame_orig = padarray(meanIM_orig, floor((size(meanIM_reg) - size(meanIM_orig))/2),nan,'both');
    fullFrame_orig = padarray(fullFrame_orig, size(meanIM_reg) - size(fullFrame_orig),nan,'pre');
    
    templateShifts = xcorr2_nans(fullFrame_orig,meanIM_reg,[0;0],25);
    
    nanPix = isnan(meanIM_reg);
    meanIM_reg(nanPix) = 0;
    
    meanIM_reg_shifted = imtranslate(meanIM_reg,round(templateShifts(2:-1:1)));
    meanIM_reg_shifted(imtranslate(nanPix,round(templateShifts(2:-1:1)))) = nan;
    
    meanIM_reg_crop = meanIM_reg_shifted(any(~isnan(fullFrame_orig),2),any(~isnan(fullFrame_orig),1));
    
    figure; imshow(imfuse(meanIM_reg_crop / prctile(meanIM_reg_crop(:),99.5),meanIM_orig / prctile(meanIM_orig(:),99.5)))
    
    [gx,gy] = gradient(meanIM_orig);
    ngs_orig(idx) = norm(sqrt(gx.^2+gy.^2),'fro');
    
    [gx,gy] = gradient(meanIM_reg_crop);
    ngs_reg(idx) = norm(sqrt(gx.^2+gy.^2),'fro');
end

%%
[H,P,CI,STATS]=ttest(ngs_orig,ngs_reg,'Tail','left');
[P,H,STATS]=signrank(ngs_orig,ngs_reg,'Tail','left');

if P < 0.001
    sig_reg = '***';
elseif P < 0.01
    sig_reg = '**';
elseif P < 0.05
    sig_reg = '*';
else 
    sig_reg = 'n.s.';
end

figure(234); hold on;
plot([ones(size(ngs_orig)),2*ones(size(ngs_reg))]',[ngs_orig,ngs_reg]','.-','Color',ones(1,3)*0.5)
errorbar(1:2,[mean(ngs_orig);mean(ngs_reg)],[std(ngs_orig);std(ngs_reg)]/sqrt(length(drs)),'r.-','MarkerSize',10,'LineWidth',1.5)

ylimits = ylim;
maxPoint = max([ngs_orig;ngs_reg]);
plot([1;2],ones(2,1)*(ylimits(2)+maxPoint)/2,'-k','LineWidth',2)
text(1.5,(ylimits(2)+maxPoint)/2,sig_reg,'FontSize',16,'HorizontalAlignment','center','VerticalAlignment','baseline')

hold off;
ylabel('Crispness Metric')
xlim([0,3])
xticks([1,2])
xticklabels({'Original','Strip Registration'})
