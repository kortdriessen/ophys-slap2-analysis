function [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, likelihood_means, log_means, ySearch, xSearch, zSearch, channels, robust)

validSPs = ~any(isnan(data),2);
sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,channels,:);
sub_log_means = log_means(ySearch,xSearch,zSearch,channels,:);

scalingFactorTable = ones(size(sub_likelihood_means,1:4));

for chIx = 1:length(channels)
    % Slice once for all y, x, z values to reduce indexing inside the loop
    expectedData = squeeze(sub_likelihood_means(:,:,:,chIx,validSPs)); % Size [y, x, z, validSPs]
    
    % Replicate tmpData to match the expectedData shape
    tmpData = reshape(data(validSPs, chIx), [1, 1, 1, length(validSPs)]);
    
    % Compute numerator and denominator using matrix operations
    numerator = sum(tmpData .* ~isnan(expectedData), 4);
    denominator = sum(expectedData, 4, 'omitnan');
    
    % Store the result in the scalingFactorTable
    scalingFactorTable(:,:,:,chIx) = numerator ./ denominator;
end

scaled_likelihood_means = scalingFactorTable .* sub_likelihood_means;
scaled_log_means = log(scalingFactorTable) + sub_log_means;

if robust
    logLikelihoodTable = sum(max(-b, min(b,(reshape(data(validSPs,:)',[1,1,1,length(channels),sum(validSPs)]) - scaled_likelihood_means) ./ sqrt(scaled_likelihood_means))).^2,[4 5],"omitnan");
    logLikelihoodTable = -abs(logLikelihoodTable(:));
    logLikelihoodTable(mean(isnan(scaled_likelihood_means),[4 5]) == 1) = nan;

else
    logLikelihoodTable = sum(reshape(data',[1,1,1,size(data')]) .* scaled_log_means - scaled_likelihood_means,[4 5],"omitnan");
    logLikelihoodTable(mean(isnan(scaled_likelihood_means),[4 5]) == 1) = nan;
end

% prevent Inf from being the max log likelihood
logLikelihoodTable(logLikelihoodTable >= Inf) = -1e10;