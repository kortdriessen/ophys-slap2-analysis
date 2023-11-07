function [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, likelihood_means, log_means, ySearch, xSearch, zSearch, robust)

sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,:);
sub_log_means = log_means(ySearch,xSearch,zSearch,:);

scalingFactorTable = ones(size(sub_likelihood_means,1:3));

for y = 1:length(ySearch)
    for x = 1:length(xSearch)
        for z = 1:length(zSearch)
            expectedData = squeeze(sub_likelihood_means(y,x,z,:));
            scalingFactorTable(y,x,z) = sum(data(~isnan(expectedData))) ./ sum(expectedData,"omitnan");
        end
    end
end

scaled_likelihood_means = scalingFactorTable .* sub_likelihood_means;
scaled_log_means = log(scalingFactorTable) + sub_log_means;

if robust
    logLikelihoodTable = sum(max(-b, min(b,(reshape(data,[1,1,1,length(data)]) - scaled_likelihood_means) ./ sqrt(scaled_likelihood_means))).^2,4,"omitnan");
    logLikelihoodTable = -abs(logLikelihoodTable(:));
else
    logLikelihoodTable = sum(reshape(data,[1,1,1,length(data)]) .* scaled_log_means - scaled_likelihood_means,4,"omitnan");
end

% prevent Inf from being the max log likelihood
logLikelihoodTable(logLikelihoodTable >= Inf) = -1e10;