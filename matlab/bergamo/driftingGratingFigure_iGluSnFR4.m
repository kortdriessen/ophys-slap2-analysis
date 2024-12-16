files8880 = {'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00001_20240912_143208\scan_00001_20240912_143208_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00003_20240912_150010\scan_00003_20240912_150010_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00006_20240912_155432\scan_00006_20240912_155432_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\754588\2024-09-12\scans\scan_00007_20240912_160302\scan_00007_20240912_160302_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00001_20240816_113139\scan_00001_20240816_113139_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00003_20240816_115328\scan_00003_20240816_115328_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00005_20240816_125104\scan_00005_20240816_125104_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750099\08-16-24\scans\scan_00006_20240816_130139\scan_00006_20240816_130139_REGISTERED_RAW.tif'};

files8360 = {'Z:\scratch\ophys\Michael\visual stim characterization\747775\08-16-24\scans\scan_00001_20240816_164019\scan_00001_20240816_164019_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\747775\08-16-24\scans\scan_00002_20240816_165246\scan_00002_20240816_165246_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\747775\08-16-24\scans\scan_00004_20240816_171100\scan_00004_20240816_171100_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\747775\08-16-24\scans\scan_00006_20240816_173309\scan_00006_20240816_173309_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\753085\2024-08-29\scans\scan_00002_20240829_122649\scan_00002_20240829_122649_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\753085\2024-08-29\scans\scan_00004_20240829_124950\scan_00004_20240829_124950_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\753085\2024-08-29\scans\scan_00008_20240829_140301\scan_00008_20240829_140301_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\753085\2024-08-29\scans\scan_00009_20240829_141207\scan_00009_20240829_141207_REGISTERED_RAW.tif'};

files9601 = {'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00001_20240816_142643\scan_00001_20240816_142643_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00004_20240816_153225\scan_00004_20240816_153225_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00005_20240816_154521\scan_00005_20240816_154521_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\visual stim characterization\750098\scans\scan_00006_20240816_155419\scan_00006_20240816_155419_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00001_20240801_110557\scan_00001_20240801_110557_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00003_20240801_113412\scan_00003_20240801_113412_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00004_20240801_114829\scan_00004_20240801_114829_REGISTERED_RAW.tif',
    'Z:\scratch\ophys\Michael\GluSourceExtractionValidation\743713\20240801\scans\scan_00006_20240801_121111\scan_00006_20240801_121111_REGISTERED_RAW.tif'};

results8880 = makeXValCovPlot(files8880);
bg_results8880 = makeXValCovPlot(files8880,1);

results8360 = makeXValCovPlot(files8360);
bg_results8360 = makeXValCovPlot(files8360,1);

results9601 = makeXValCovPlot(files9601);
bg_results9601 = makeXValCovPlot(files9601,1);


%%
colors = distinguishable_colors(6);

figure(123); clf;
figure(124); clf;
plotResults(results8880, 'v8880', colors(1,:), 123, 124);
plotResults(results8360, 'v8360', colors(2,:), 123, 124);
plotResults(results9601, 'v9601', colors(3,:), 123, 124);

plotResults(bg_results8880, 'v8880', colors(4,:), 123, 124);
plotResults(bg_results8360, 'v8360', colors(5,:), 123, 124);
plotResults(bg_results9601, 'v9601', colors(6,:), 123, 124);

%%
figure(223); clf;
figure(224); clf;
figure(225); clf;
plotByFOV(results9601, files9601, 'v9601 (by FOV)', 223);
plotByFOV(results8880, files8880, 'v8880 (by FOV)', 224);
plotByFOV(results8360, files8360, 'v8360 (by FOV)', 225);

%%
function plotResults(results, label, color, figRaw, figNorm)
    % Compute unique distances
    [unique_distances, mean_vals, norm_mean_vals, upper_bound, lower_bound, ...
        norm_upper_bound, norm_lower_bound] = calculateStats(results);

    % Plot raw covariance
    figure(figRaw);
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], ...
         [upper_bound; flipud(lower_bound)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(2.^unique_distances, mean_vals, 'Color', color, 'LineWidth', 2, 'DisplayName', label); % Add line to legend
    hold off;

    % Plot normalized covariance
    figure(figNorm);
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], ...
         [norm_upper_bound; flipud(norm_lower_bound)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(2.^unique_distances, norm_mean_vals, 'Color', color, 'LineWidth', 2, 'DisplayName', label); % Add line to legend
    hold off;

    % Add legend dynamically after all plots
    figure(figRaw);
    legend('show', 'Interpreter', 'none');

    figure(figNorm);
    legend('show', 'Interpreter', 'none');
end

function plotByFOV(results, files, titleText, figNum)
    colors = distinguishable_colors(length(files));
    figure(figNum);
    hold on;

    for fileNum = 1:length(files)
        % Filter results by file index
        idx = results.fileIx == fileNum;
        [unique_distances, mean_vals, ~, upper_bound, lower_bound] = calculateStats(results, idx);

        % Plot data for each FOV
        fill(2.^[unique_distances; flipud(unique_distances)], ...
             [upper_bound; flipud(lower_bound)], ...
             colors(fileNum, :), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '');
        plot(2.^unique_distances, mean_vals, 'Color', colors(fileNum, :), 'LineWidth', 2);
    end

    title(titleText);
    hold off;
end

function [unique_distances, mean_vals, norm_mean_vals, upper_bound, lower_bound, ...
          norm_upper_bound, norm_lower_bound] = calculateStats(results, mask)
    if nargin < 2
        mask = true(size(results.distance));
    end

    % Compute unique distances and group stats
    [unique_distances, ~, idx] = unique(round(log2(results.distance(mask)) * 50) / 50);
    mean_vals = accumarray(idx, results.xval_cov(mask), [], @mean);
    norm_mean_vals = mean_vals ./ max(mean_vals);

    std_vals = accumarray(idx, results.xval_cov(mask), [], @std);
    count_vals = accumarray(idx, results.xval_cov(mask), [], @numel);
    sem_vals = std_vals ./ sqrt(count_vals);

    upper_bound = mean_vals + sem_vals;
    lower_bound = mean_vals - sem_vals;

    norm_sem_vals = sem_vals ./ max(mean_vals);
    norm_upper_bound = norm_mean_vals + norm_sem_vals;
    norm_lower_bound = norm_mean_vals - norm_sem_vals;
end

function fillPlot(unique_distances, upper_bound, lower_bound, color)
    % Helper function to fill plots with transparency
    fill(2.^[unique_distances; flipud(unique_distances)], ...
         [upper_bound; flipud(lower_bound)], ...
         color, 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'DisplayName', '');
end