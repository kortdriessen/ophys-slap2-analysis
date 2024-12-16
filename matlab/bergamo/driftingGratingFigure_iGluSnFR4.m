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

results8360 = makeXValCovPlot(files8360);

results9601 = makeXValCovPlot(files9601);


%%
[unique_distances, ~, idx] = unique(round(log2(results8880.distance)*50)/50);

mean_vals = accumarray(idx, results8880.xval_cov, [], @mean);
norm_mean_vals = mean_vals ./ max(mean_vals);

std_vals = accumarray(idx, results8880.xval_cov, [], @std);
norm_std_vals = accumarray(idx, results8880.xval_cov ./ max(mean_vals), [], @std);

count_vals = accumarray(idx, results8880.xval_cov, [], @numel);

sem_vals = std_vals ./ sqrt(count_vals);
norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

upper_bound = mean_vals + sem_vals;
lower_bound = mean_vals - sem_vals;
norm_upper_bound = norm_mean_vals + norm_sem_vals;
norm_lower_bound = norm_mean_vals - norm_sem_vals;

figure(123)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, mean_vals, 'b', 'LineWidth', 2);
xlabel('Distance');
ylabel('Covariance');
box on;
hold off;

figure(124)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [norm_upper_bound; flipud(norm_lower_bound)], ...
'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, norm_mean_vals, 'b', 'LineWidth', 2);
xlabel('Distance');
ylabel('Normalized Covariance');
box on;
hold off;

[unique_distances, ~, idx] = unique(round(log2(results8360.distance)*50)/50);
mean_vals = accumarray(idx, results8360.xval_cov, [], @mean);
norm_mean_vals = mean_vals ./ max(mean_vals);

std_vals = accumarray(idx, results8360.xval_cov, [], @std);
norm_std_vals = accumarray(idx, results8360.xval_cov ./ max(mean_vals), [], @std);

count_vals = accumarray(idx, results8360.xval_cov, [], @numel);

sem_vals = std_vals ./ sqrt(count_vals);
norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

upper_bound = mean_vals + sem_vals;
lower_bound = mean_vals - sem_vals;
norm_upper_bound = norm_mean_vals + norm_sem_vals;
norm_lower_bound = norm_mean_vals - norm_sem_vals;

figure(123)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, mean_vals, 'm', 'LineWidth', 2);
box on;
hold off;

figure(124)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [norm_upper_bound; flipud(norm_lower_bound)], ...
'm', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, norm_mean_vals, 'm', 'LineWidth', 2);
box on;
hold off;

[unique_distances, ~, idx] = unique(round(log2(results9601.distance)*50)/50);
mean_vals = accumarray(idx, results9601.xval_cov, [], @mean);
norm_mean_vals = mean_vals ./ max(mean_vals);

std_vals = accumarray(idx, results9601.xval_cov, [], @std);
norm_std_vals = accumarray(idx, results9601.xval_cov ./ max(mean_vals), [], @std);

count_vals = accumarray(idx, results9601.xval_cov, [], @numel);

sem_vals = std_vals ./ sqrt(count_vals);
norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

upper_bound = mean_vals + sem_vals;
lower_bound = mean_vals - sem_vals;
norm_upper_bound = norm_mean_vals + norm_sem_vals;
norm_lower_bound = norm_mean_vals - norm_sem_vals;

figure(123)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, mean_vals, 'r', 'LineWidth', 2);
box on;
hold off;

legend('','v8880','','v8360','','v9601')

figure(124)
hold on;
fill(2.^[unique_distances; flipud(unique_distances)], [norm_upper_bound; flipud(norm_lower_bound)], ...
'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
plot(2.^unique_distances, norm_mean_vals, 'r', 'LineWidth', 2);
box on;
hold off;

legend('','v8880','','v8360','','v9601')

%%
colors = distinguishable_colors(length(files9601));
for fileNum = 1:length(files9601)
    [unique_distances, ~, idx] = unique(round(log2(results9601.distance(results9601.fileIx==fileNum))*50)/50);
    mean_vals = accumarray(idx, results9601.xval_cov(results9601.fileIx==fileNum), [], @mean);
    norm_mean_vals = mean_vals ./ max(mean_vals);

    std_vals = accumarray(idx, results9601.xval_cov(results9601.fileIx==fileNum), [], @std);
    norm_std_vals = accumarray(idx, results9601.xval_cov(results9601.fileIx==fileNum) ./ max(mean_vals), [], @std);

    count_vals = accumarray(idx, results9601.xval_cov(results9601.fileIx==fileNum), [], @numel);

    sem_vals = std_vals ./ sqrt(count_vals);
    norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

    upper_bound = mean_vals + sem_vals;
    lower_bound = mean_vals - sem_vals;
    norm_upper_bound = norm_mean_vals + norm_sem_vals;
    norm_lower_bound = norm_mean_vals - norm_sem_vals;

    figure(223)
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
        colors(fileNum,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(2.^unique_distances, mean_vals, 'Color',colors(fileNum,:), 'LineWidth', 2);
    box on;
    hold off;
end

title('v9601, by FOV (8 total)')

%%
colors = distinguishable_colors(length(files8880));
for fileNum = 1:length(files8880)
    [unique_distances, ~, idx] = unique(round(log2(results8880.distance(results8880.fileIx==fileNum))*50)/50);
    mean_vals = accumarray(idx, results8880.xval_cov(results8880.fileIx==fileNum), [], @mean);
    norm_mean_vals = mean_vals ./ max(mean_vals);

    std_vals = accumarray(idx, results8880.xval_cov(results8880.fileIx==fileNum), [], @std);
    norm_std_vals = accumarray(idx, results8880.xval_cov(results8880.fileIx==fileNum) ./ max(mean_vals), [], @std);

    count_vals = accumarray(idx, results8880.xval_cov(results8880.fileIx==fileNum), [], @numel);

    sem_vals = std_vals ./ sqrt(count_vals);
    norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

    upper_bound = mean_vals + sem_vals;
    lower_bound = mean_vals - sem_vals;
    norm_upper_bound = norm_mean_vals + norm_sem_vals;
    norm_lower_bound = norm_mean_vals - norm_sem_vals;

    figure(224)
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
        colors(fileNum,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(2.^unique_distances, mean_vals, 'Color',colors(fileNum,:), 'LineWidth', 2);
    box on;
    hold off;
end

title('v8880, by FOV (8 total)')

%%
colors = distinguishable_colors(length(files8360));
for fileNum = 1:length(files8360)
    [unique_distances, ~, idx] = unique(round(log2(results8360.distance(results8360.fileIx==fileNum))*50)/50);
    mean_vals = accumarray(idx, results8360.xval_cov(results8360.fileIx==fileNum), [], @mean);
    norm_mean_vals = mean_vals ./ max(mean_vals);

    std_vals = accumarray(idx, results8360.xval_cov(results8360.fileIx==fileNum), [], @std);
    norm_std_vals = accumarray(idx, results8360.xval_cov(results8360.fileIx==fileNum) ./ max(mean_vals), [], @std);

    count_vals = accumarray(idx, results8360.xval_cov(results8360.fileIx==fileNum), [], @numel);

    sem_vals = std_vals ./ sqrt(count_vals);
    norm_sem_vals = norm_std_vals ./ sqrt(count_vals);

    upper_bound = mean_vals + sem_vals;
    lower_bound = mean_vals - sem_vals;
    norm_upper_bound = norm_mean_vals + norm_sem_vals;
    norm_lower_bound = norm_mean_vals - norm_sem_vals;

    figure(225)
    hold on;
    fill(2.^[unique_distances; flipud(unique_distances)], [upper_bound; flipud(lower_bound)], ...
        colors(fileNum,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    plot(2.^unique_distances, mean_vals, 'Color',colors(fileNum,:), 'LineWidth', 2);
    box on;
    hold off;
end

title('v8360, by FOV (8 total)')

