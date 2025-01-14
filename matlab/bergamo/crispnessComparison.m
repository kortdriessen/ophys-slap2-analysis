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
ngs_suite2p = nan(size(drs));
ngs_caiman = nan(size(drs));
ngs_orig = nan(size(drs));

CM_reg = cell(size(drs));
CM_suite2p = cell(size(drs));
CM_caiman = cell(size(drs));
CM_orig = cell(size(drs));

for idx = 1:length(drs)
    disp(['Running file ' num2str(idx)])

    basename = drs{idx}(regexp(drs{idx}, 'scan_[0-9]+'):end);

    mov_orig_fn = fullfile(drs{idx},[basename '.tif']);
    [mov_orig, ~, ~] = networkScanImageTiffReader(mov_orig_fn);
    mov_orig = mov_orig(:,5:end,2:2:end);

    load(fullfile(drs{idx},[basename '_ALIGNMENTDATA.mat']))

    motionC_caiman = -h5read(fullfile(drs{idx},[basename '_CAIMAN_ADATA.h5']),'/C');
    motionC_caiman = motionC_caiman - mean(motionC_caiman) + mean(aData.motionC);
    motionR_caiman = -h5read(fullfile(drs{idx},[basename '_CAIMAN_ADATA.h5']),'/R');
    motionR_caiman = motionR_caiman - mean(motionR_caiman) + mean(aData.motionR);

    motionC_suite2p = h5read(fullfile(drs{idx},[basename '_SUITE2P_ADATA.h5']),'/R');
    motionC_suite2p = motionC_suite2p(2:2:end) - mean(motionC_suite2p(2:2:end)) + mean(aData.motionC);
    motionR_suite2p = h5read(fullfile(drs{idx},[basename '_SUITE2P_ADATA.h5']),'/C');
    motionR_suite2p = motionR_suite2p(2:2:end) - mean(motionR_suite2p(2:2:end)) + mean(aData.motionR);

    nframes = min([length(aData.motionC),length(motionC_caiman),length(motionC_suite2p)]);
    
    mov_reg = nan([size(mov_orig,1:2) nframes]);
    mov_caiman = nan([size(mov_orig,1:2) nframes]);
    mov_suite2p = nan([size(mov_orig,1:2) nframes]);

    for fix = 1:nframes
        mov_reg(:,:,fix) = imtranslate(mov_orig(:,:,fix),-[aData.motionR(fix) aData.motionC(fix)],'FillValues',nan);
        mov_caiman(:,:,fix) = imtranslate(mov_orig(:,:,fix),-[motionR_caiman(fix) motionC_caiman(fix)],'FillValues',nan);
        mov_suite2p(:,:,fix) = imtranslate(mov_orig(:,:,fix),-[motionR_suite2p(fix) motionC_suite2p(fix)],'FillValues',nan);
    end

    meanIM_reg = mean(mov_reg,3,'omitnan');
    meanIM_suite2p = mean(mov_suite2p,3,'omitnan');
    meanIM_caiman = mean(mov_caiman,3,'omitnan');
    meanIM_orig = mean(mov_orig,3,'omitnan');

    figure; subplot(1,4,1); imagesc(meanIM_orig); axis image; colormap gray; %imshow(imfuse(meanIM_reg_crop / prctile(meanIM_reg_crop(:),99.5),meanIM_orig / prctile(meanIM_orig(:),99.5)))
    subplot(1,4,2); imagesc(meanIM_reg); axis image; colormap gray;
    subplot(1,4,3); imagesc(meanIM_caiman); axis image; colormap gray;
    subplot(1,4,4); imagesc(meanIM_suite2p); axis image; colormap gray;
    drawnow;

    [gx,gy] = gradient(meanIM_orig);
    ngs_orig(idx) = norm(sqrt(gx.^2+gy.^2),'fro');
    
    [gx,gy] = gradient(meanIM_reg);
    ngs_reg(idx) = norm(sqrt(gx.^2+gy.^2),'fro');

    [gx,gy] = gradient(meanIM_caiman);
    ngs_caiman(idx) = norm(sqrt(gx.^2+gy.^2),'fro');

    [gx,gy] = gradient(meanIM_suite2p);
    ngs_suite2p(idx) = norm(sqrt(gx.^2+gy.^2),'fro');
    
    CM_orig{idx} = nan(nframes,1);
    CM_reg{idx} = nan(nframes,1);
    CM_caiman{idx} = nan(nframes,1);
    CM_suite2p{idx} = nan(nframes,1);
    
    mov_orig = reshape(mov_orig(:,:,1:nframes),[],nframes);
    mov_reg = reshape(mov_reg,[],nframes);
    mov_caiman = reshape(mov_caiman,[],nframes);
    mov_suite2p = reshape(mov_suite2p,[],nframes);

    for fix = 1:nframes
        CM_orig{idx}(fix) = corr(mov_orig(:,fix),meanIM_orig(:),'rows','pairwise');
        CM_orig{idx}(fix) = corr(mov_reg(:,fix),meanIM_reg(:),'rows','pairwise');
        CM_orig{idx}(fix) = corr(mov_caiman(:,fix),meanIM_caiman(:),'rows','pairwise');
        CM_orig{idx}(fix) = corr(mov_suite2p(:,fix),meanIM_suite2p(:),'rows','pairwise');
    end
end

%%
figure(233); hold on;
pct_crisp_inc_reg = (ngs_reg ./ ngs_orig - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_reg)),1,'Strip Registration');
boxchart(model_labels,pct_crisp_inc_reg,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_reg,'o','filled','b');

pct_crisp_inc_suite2p = (ngs_suite2p ./ ngs_orig - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_suite2p)),1,'suite2p');
boxchart(model_labels,pct_crisp_inc_suite2p,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_suite2p,'o','filled','g');

pct_crisp_inc_caiman = (ngs_caiman ./ ngs_orig - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_caiman)),1,'CaImAn');
boxchart(model_labels,pct_crisp_inc_caiman,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_caiman,'o','filled','r');

% set(gca,'YDir','reverse')

yline(0,'--')

ylims = ylim;
ylim([-100, ylims(2)])

ylabel('% increase in crispness from original')
% xlim([0.5,1.5])
% xticks([1])
% xticklabels({'Strip Registration'})


%%
figure(234); hold on;
pct_crisp_inc_reg = (ngs_reg ./ ngs_orig - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_reg)),1,'Unregistered');
boxchart(model_labels,pct_crisp_inc_reg,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_reg,'o','filled','b');

pct_crisp_inc_suite2p = (ngs_reg ./ ngs_suite2p - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_suite2p)),1,'suite2p');
boxchart(model_labels,pct_crisp_inc_suite2p,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_suite2p,'o','filled','g');

pct_crisp_inc_caiman = (ngs_reg ./ ngs_caiman - 1) * 100;
model_labels = categorical(ones(size(pct_crisp_inc_caiman)),1,'CaImAn');
boxchart(model_labels,pct_crisp_inc_caiman,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,pct_crisp_inc_caiman,'o','filled','r');

% set(gca,'YDir','reverse')

yline(0,'--')

ylims = ylim;
ylim([-100, ylims(2)])

ylabel('% increase in crispness of strip registration from []')
% xlim([0.5,1.5])
% xticks([1])
% xticklabels({'Strip Registration'})

%%
figure(235); hold on;
mCM_reg = cellfun(@(x) mean(x,'all'), CM_reg);
model_labels = categorical(ones(size(mCM_reg)),1,'Strip Registration');
boxchart(model_labels,mCM_reg,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_reg,'o','filled','b');

mCM_suite2p = cellfun(@(x) mean(x,'all'), CM_suite2p);
model_labels = categorical(ones(size(mCM_suite2p)),1,'suite2p');
boxchart(model_labels,mCM_suite2p,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_suite2p,'o','filled','g');

mCM_caiman = cellfun(@(x) mean(x,'all'), CM_caiman);
model_labels = categorical(ones(size(mCM_caiman)),1,'CaImAn');
boxchart(model_labels,mCM_caiman,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_caiman,'o','filled','r');

% set(gca,'YDir','reverse')

yline(0,'--')

ylims = ylim;
% ylim([-100, ylims(2)])

ylabel('mean correlation to the mean image')
% xlim([0.5,1.5])
% xticks([1])
% xticklabels({'Strip Registration'})

%%
figure(236); hold on;
mCM_reg = cellfun(@(x) sum(1-x), CM_reg);
model_labels = categorical(ones(size(mCM_reg)),1,'Strip Registration');
boxchart(model_labels,mCM_reg,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_reg,'o','filled','b');

mCM_suite2p = cellfun(@(x) sum(1-x), CM_suite2p);
model_labels = categorical(ones(size(mCM_suite2p)),1,'suite2p');
boxchart(model_labels,mCM_suite2p,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_suite2p,'o','filled','g');

mCM_caiman = cellfun(@(x) sum(1-x), CM_caiman);
model_labels = categorical(ones(size(mCM_caiman)),1,'CaImAn');
boxchart(model_labels,mCM_caiman,'Notch','on','BoxFaceColor','none','BoxEdgeColor','k');
swarmchart(model_labels,mCM_caiman,'o','filled','r');

% set(gca,'YDir','reverse')

yline(0,'--')

ylims = ylim;
% ylim([-100, ylims(2)])

ylabel('integrated 1-CM')
% xlim([0.5,1.5])
% xticks([1])
% xticklabels({'Strip Registration'})

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
