function evaluateSummary(exptSummary, GT)

params = exptSummary.params;
denoiseWindow = params.denoiseWindow_samps*params.dsFac+1;
baselineWindow = ceil(params.baselineWindow_Glu_s/params.frametime);

%first, match sources from GT to extracted
[rr,cc] = ndgrid(1:size(exptSummary.footprints,1), 1:size(exptSummary.footprints,2));
for sourceIx = 1:size(exptSummary.footprints,3)
    sel = exptSummary.footprints(:,:,sourceIx,1);
    sel(isnan(sel)) = 0;
    est(sourceIx,1) = sum(rr.*sel, 'all')./sum(sel, 'all');
    est(sourceIx,2) = sum(cc.*sel, 'all')./sum(sel, 'all');
end

ROIs = permute(GT{1}.ROIs, [2 1 3]);
[rr,cc] = ndgrid(1:size(ROIs,1), 1:size(ROIs,2));
for sourceIx = 1:size(ROIs,3)
    sel = ROIs(:,:,sourceIx);
    gt(sourceIx,1) = sum(rr.*sel, 'all')./sum(sel, 'all') - exptSummary.perTrialAlignmentOffsets(1,1);
    gt(sourceIx,2) = sum(cc.*sel, 'all')./sum(sel, 'all') - exptSummary.perTrialAlignmentOffsets(2,1);
end

for iter = 1:3
    D = pdist2(est,gt,'euclidean');
    [~, closestGTforEachEst] = min(D,[],2);
    [~, closestEstforEachGT] = min(D,[],1);
    dR = median(est(closestEstforEachGT,1) - gt(:,1));
    dC =median(est(closestEstforEachGT,2) - gt(:,2));
    gt = gt + [dR dC];
end
D = pdist2(est,gt,'euclidean');
[~, closestGTforEachEst] = min(D,[],2);
[~, closestEstforEachGT] = min(D,[],1);

%Scatter the loclaizations vs ground truth
figure, imagesc(exptSummary.meanIM);
hold on, scatter(gt(:,2), gt(:,1), 'g', 'filled'); hold on, scatter(est(:,2), est(:,1), 'm', 'filled')
for ix = 1:size(gt,1)
    plot([gt(ix,2) est(closestEstforEachGT(ix),2)], [gt(ix,1) est(closestEstforEachGT(ix),1)], 'g', 'linewidth', 2)
end

%plot traces
for trial = 1:length(GT)
    %compute the correlation of gt to estimated acivity at each GT site
    ESTtmp = exptSummary.dFF{trial}(closestEstforEachGT,:);
    selFrames = ~all(isnan(ESTtmp),1);
    ESTtmp(isnan(ESTtmp)) = 0;
    C_est = corr(GT{trial}.activity(:,selFrames)', ESTtmp(:, selFrames)');

    ROIF = GT{trial}.ROI_activity;
    ROIF0 = computeF0(ROIF',denoiseWindow, baselineWindow, 1)';
    GTtmp = (ROIF-ROIF0)./ROIF0;
    GTtmp(isnan(GTtmp)) = 0;
    C_roi = corr(GT{trial}.activity(:,selFrames)', GTtmp(:,selFrames)');

    figure,
    subplot(1,2,1);
    imagesc(C_roi.^2, [0 1]);
    xlabel('optimal linear ROI'); ylabel('GT');
    subplot(1,2,2);
    imagesc(C_est.^2, [0 1]);
    xlabel('CNMF'); ylabel('GT');
    cb = colorbar;
    ylabel(cb, 'R^2')

    GTnorm = GT{trial}.activity./sqrt(mean(GT{trial}.activity.^2, 2));
    ESTnorm = exptSummary.dFF{trial}./sqrt(mean(exptSummary.dFF{trial}.^2, 2, 'omitnan'));
    ESTnorm = ESTnorm(closestEstforEachGT,:);
    ROInorm = GTtmp./sqrt(mean(GTtmp.^2,2));

    figure,
    h1 =plot(GTnorm' + (1:size(GTnorm,1))*5, 'g', 'linewidth', 2);
    hold on,
    h3 = plot(ROInorm' + (1:size(ROInorm,1))*5, 'color', [0.8 0.8 0.8]);
    h2 = plot(ESTnorm' + (1:size(ESTnorm,1))*5, 'm', 'linewidth', 2);

    legend ([h1(1) h2(1) h3(1)], {'Ground Truth', 'Estimate', 'Optimal linear ROI'})
end