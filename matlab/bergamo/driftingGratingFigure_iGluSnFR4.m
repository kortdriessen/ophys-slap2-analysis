%% generate data tables from raw movies

files8880 = {'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00001_20240912_143208\scan_00001_20240912_143208_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00003_20240912_150010\scan_00003_20240912_150010_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00006_20240912_155432\scan_00006_20240912_155432_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00007_20240912_160302\scan_00007_20240912_160302_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00002_20240816_114130\scan_00002_20240816_114130_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00003_20240816_115328\scan_00003_20240816_115328_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00007_20240816_131021\scan_00007_20240816_131021_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00006_20240816_130139\scan_00006_20240816_130139_DENOISED_REGISTERED_RAW.tif'};

filesGlu3 = {'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755823\2024-09-11\scans\scan_00010_20240911_135849\scan_00010_20240911_135849_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755823\2024-09-11\scans\scan_00009_20240911_135003\scan_00009_20240911_135003_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755823\2024-09-11\scans\scan_00008_20240911_134125\scan_00008_20240911_134125_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755823\2024-09-11\scans\scan_00007_20240911_133234\scan_00007_20240911_133234_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755826\2024-09-18\scans\scan_00008_20240918_112104\scan_00008_20240918_112104_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755826\2024-09-18\scans\scan_00005_20240918_104519\scan_00005_20240918_104519_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755826\2024-09-18\scans\scan_00006_20240918_105528\scan_00006_20240918_105528_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\755826\2024-09-18\scans\scan_00002_20240918_100740\scan_00002_20240918_100740_DENOISED_REGISTERED_RAW.tif'};

files9601 = {'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00001_20240816_142643\scan_00001_20240816_142643_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00004_20240816_153225\scan_00004_20240816_153225_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00005_20240816_154521\scan_00005_20240816_154521_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00006_20240816_155419\scan_00006_20240816_155419_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00001_20240801_110557\scan_00001_20240801_110557_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00003_20240801_113412\scan_00003_20240801_113412_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00004_20240801_114829\scan_00004_20240801_114829_DENOISED_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00006_20240801_121111\scan_00006_20240801_121111_DENOISED_REGISTERED_RAW.tif'};

HPfreqs = 0;

% results data tables: each row represents a pair of pixels in a movie,
% columns represent variables such as distance between pixels, xval_cov,
% file index, pixel identities, etc.
results8880 = makeXValCovPlot(files8880, HPfreqs);
bg_results8880 = makeXValCovPlot(files8880, HPfreqs, 1);

resultsGlu3 = makeXValCovPlot(filesGlu3, HPfreqs);
bg_resultsGlu3 = makeXValCovPlot(filesGlu3, HPfreqs, 1);

results9601 = makeXValCovPlot(files9601, HPfreqs);
bg_results9601 = makeXValCovPlot(files9601, HPfreqs, 1);

save('Z:\scratch\ophys\Michael\figures\iglusnfr4 paper\xval_cov_results_20250125.mat','-v7.3');

%% plot xval_xcov (both normalized to max vs unnormalized) vs. distance between pixels
colors = [0 0 1; 0.6350 0.0780 0.1840; 0.4660 0.6740 0.1880; 0 0.4470 0.7410; 1 0 0; 0 1 0];

figure(123); clf;
figure(124); clf;

glu3_norm = plotResults(resultsGlu3, 'iGluSnFR3 v857', colors(1,:), 123, 124);
glu4f_norm = plotResults(results9601, 'iGluSnFR4 v9601', colors(2,:), 123, 124);
glu4s_norm = plotResults(results8880, 'iGluSnFR4 v8880', colors(3,:), 123, 124);

plotResults(bg_resultsGlu3, 'iGluSnFR3 v857 (background)', colors(4,:), 123, 124, glu3_norm);
plotResults(bg_results9601, 'iGluSnFR4 v9601 (background)', colors(5,:), 123, 124, glu4f_norm);
plotResults(bg_results8880, 'iGluSnFR4 v8880 (background)', colors(6,:), 123, 124, glu4s_norm);

figure(123);
xlabel('distance (um)')
ylabel('xval xcov')

figure(124);
xlabel('distance (um)')
ylabel('xval xcov (normalized to max of foreground)')

%% bootstrap to get 1/E^2 decay distance and test for difference between variants
rng(2,"twister");

fileIDsGlu3 = unique(resultsGlu3.fileIx);
fileIDs9601 = unique(results9601.fileIx);
fileIDs8880 = unique(results8880.fileIx);

file_mean_vals_glu3 = {};
file_mean_vals_9601 = {};
file_mean_vals_8880 = {};

for i = 1:length(fileIDsGlu3)
    % Compute unique distances
    [unique_distances_glu3, file_mean_vals_glu3{i}, ~, ~, ~, ...
        ~, ~, ~] = calculateStats(resultsGlu3(resultsGlu3.fileIx == fileIDsGlu3(i),:), 0);
end
file_mean_vals_glu3 = cell2mat(file_mean_vals_glu3);

for i = 1:length(fileIDs8880)
    % Compute unique distances
    [unique_distances_8880, file_mean_vals_8880{i}, ~, ~, ~, ...
        ~, ~, ~] = calculateStats(results8880(results8880.fileIx == fileIDs8880(i),:), 0);
end
file_mean_vals_8880 = cell2mat(file_mean_vals_8880);

for i = 1:length(fileIDs9601)
    % Compute unique distances
    [unique_distances_9601, file_mean_vals_9601{i}, ~, ~, ~, ...
        ~, ~, ~] = calculateStats(results9601(results9601.fileIx == fileIDs9601(i),:), 0);
end
file_mean_vals_9601 = cell2mat(file_mean_vals_9601);

decay_dist_glu3 = findDecayByE2(2.^unique_distances_glu3, mean(file_mean_vals_glu3,2));
decay_dist_9601 = findDecayByE2(2.^unique_distances_9601, mean(file_mean_vals_9601,2));
decay_dist_8880 = findDecayByE2(2.^unique_distances_8880, mean(file_mean_vals_8880,2));

B = 1000;

bootstrap_decay_glu3 = nan(B,1);
bootstrap_decay_9601 = nan(B,1);
bootstrap_decay_8880 = nan(B,1);

for ix = 1:B
    trialIdxsGlu3 = randi(length(filesGlu3),1,length(filesGlu3));
    trialIdxs9601 = randi(length(files9601),1,length(files9601));
    trialIdxs8880 = randi(length(files8880),1,length(files8880));

    bootstrap_decay_glu3(ix) = findDecayByE2(2.^unique_distances_glu3, mean(file_mean_vals_glu3(:,trialIdxsGlu3),2));
    bootstrap_decay_9601(ix) = findDecayByE2(2.^unique_distances_9601, mean(file_mean_vals_9601(:,trialIdxs9601),2));
    bootstrap_decay_8880(ix) = findDecayByE2(2.^unique_distances_8880, mean(file_mean_vals_8880(:,trialIdxs8880),2));
end

figure; hold on;
errorbar(1:3,[decay_dist_glu3; decay_dist_8880; decay_dist_9601], ...
    [decay_dist_glu3; decay_dist_8880; decay_dist_9601] - [prctile(bootstrap_decay_glu3,2.5); prctile(bootstrap_decay_8880,2.5); prctile(bootstrap_decay_9601,2.5)], ...
    [prctile(bootstrap_decay_glu3,97.5); prctile(bootstrap_decay_8880,97.5); prctile(bootstrap_decay_9601,97.5)] - [decay_dist_glu3; decay_dist_8880; decay_dist_9601], '.')
xlim([0.5,3.5])
xticks(1:3)
xticklabels({'857','8880','9601'})

p_8880_9601 = sum((bootstrap_decay_8880 - bootstrap_decay_9601) < 0) / B;
p_glu3_9601 = sum((bootstrap_decay_glu3 - bootstrap_decay_9601) < 0) / B;
p_glu3_8880 = sum((bootstrap_decay_glu3 - bootstrap_decay_8880) < 0) / B;

ylimits = ylim;
maxPoint = max([prctile(bootstrap_decay_glu3,97.5); prctile(bootstrap_decay_8880,97.5); prctile(bootstrap_decay_9601,97.5)]);

if abs(ylimits(2)-maxPoint) < diff(ylimits)/5
    ylim([ylimits(1) ylimits(2)+diff(ylimits)/5]);
    ylimits = ylim;
end

plot([1;2],ones(2,1)*((ylimits(2)-maxPoint)/4+maxPoint),'-k','LineWidth',1)
text(1.5,((ylimits(2)-maxPoint)/4+maxPoint),sprintf("p = %0.3f",p_glu3_8880),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','baseline')
plot([2;3],ones(2,1)*((ylimits(2)-maxPoint)/2+maxPoint),'-k','LineWidth',1)
text(2.5,((ylimits(2)-maxPoint)/2+maxPoint),sprintf("p = %0.3f",p_8880_9601),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','baseline')
plot([1;3],ones(2,1)*((ylimits(2)-maxPoint)/4*3+maxPoint),'-k','LineWidth',1)
text(2,((ylimits(2)-maxPoint)/4*3+maxPoint),sprintf("p = %0.3f",p_glu3_9601),'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','baseline')

%%
function decay_dist = findDecayByE2(X, Y)
difference = Y - max(Y)*exp(-2);

idx = find(difference(1:end-1) .* difference(2:end) < 0);
idx = idx(1);

decay_dist = X(idx) - difference(idx) .* (X(idx+1) - X(idx)) ./ (difference(idx+1) - difference(idx));

end

%%
function norm_factor = plotResults(results, label, color, figRaw, figNorm, norm_factor)
    if nargin < 6
        norm_factor = 0;
    end

    fileIDs = unique(results.fileIx);
    file_mean_vals = {};

    for i = 1:length(fileIDs)
        % Compute unique distances
        [unique_distances, file_mean_vals{i}, norm_mean_vals, upper_bound, lower_bound, ...
            norm_upper_bound, norm_lower_bound, ~] = calculateStats(results(results.fileIx == fileIDs(i),:), norm_factor);
    end
    file_mean_vals = cell2mat(file_mean_vals);
    if norm_factor == 0
        norm_factor = mean(file_mean_vals(1,:));
    end

    % Plot raw covariance
    figure(figRaw);
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], ...
         [mean(file_mean_vals,2)+std(file_mean_vals,[],2)./sqrt(size(file_mean_vals,2)); flipud(mean(file_mean_vals,2)-std(file_mean_vals,[],2)./sqrt(size(file_mean_vals,2)))], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(2.^unique_distances, mean(file_mean_vals,2), 'Color', color, 'LineWidth', 2, 'DisplayName', label); % Add line to legend
    hold off;

    % Plot normalized covariance
    figure(figNorm);
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], ...
         [mean(file_mean_vals,2)+std(file_mean_vals,[],2)./sqrt(size(file_mean_vals,2)); flipud(mean(file_mean_vals,2)-std(file_mean_vals,[],2)./sqrt(size(file_mean_vals,2)))]./norm_factor, ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(2.^unique_distances, mean(file_mean_vals,2)./norm_factor, 'Color', color, 'LineWidth', 2, 'DisplayName', label); % Add line to legend
    hold off;

    % Add legend dynamically after all plots
    figure(figRaw);
    legend('show', 'Interpreter', 'none');

    figure(figNorm);
    legend('show', 'Interpreter', 'none');
end

function [unique_distances, mean_vals, norm_mean_vals, upper_bound, lower_bound, ...
          norm_upper_bound, norm_lower_bound, norm_factor] = calculateStats(results, norm_factor, mask)
    if nargin < 2
        norm_factor = 0;
    end
    if nargin < 3
        mask = true(size(results.distance));
    end

    % Compute unique distances and group stats
    [unique_distances, ~, idx] = unique(round(log2(results.distance(mask)) * 25) / 25);
    mean_vals = accumarray([repmat(idx,[size(results.xval_cov,2) 1]),repelem(reshape(1:size(results.xval_cov,2),[],1),size(idx,1),1)], ...
        reshape(results.xval_cov(find(mask),:),[],1), [], @mean);
    if norm_factor == 0
        norm_factor = max(mean_vals);
    end
    norm_mean_vals = mean_vals ./ norm_factor;

    std_vals = accumarray([repmat(idx,[size(results.xval_cov,2) 1]),repelem(reshape(1:size(results.xval_cov,2),[],1),size(idx,1),1)], ...
        reshape(results.xval_cov(find(mask),:),[],1), [], @std);
    count_vals = accumarray([repmat(idx,[size(results.xval_cov,2) 1]),repelem(reshape(1:size(results.xval_cov,2),[],1),size(idx,1),1)], ...
        reshape(results.xval_cov(find(mask),:),[],1), [], @numel);
    sem_vals = std_vals ./ sqrt(count_vals);

    upper_bound = mean_vals + sem_vals;
    lower_bound = mean_vals - sem_vals;

    norm_sem_vals = sem_vals ./ norm_factor;
    norm_upper_bound = norm_mean_vals + norm_sem_vals;
    norm_lower_bound = norm_mean_vals - norm_sem_vals;
end