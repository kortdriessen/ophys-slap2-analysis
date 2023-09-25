function [logLikelihoodTable, scalingFactorTable] = poissonLogLikelihoodTable(data, likelihood_means, log_means, ySearch, xSearch, zSearch, robust)

sub_likelihood_means = likelihood_means(ySearch,xSearch,zSearch,:);
sub_log_means = log_means(ySearch,xSearch,zSearch,:);

scalingFactorTable = sum(data) ./ sum(sub_likelihood_means,4);

scaled_likelihood_means = scalingFactorTable .* sub_likelihood_means;
scaled_log_means = log(scalingFactorTable) + sub_log_means;

if robust
    logLikelihoodTable = sum(max(-b, min(b,(reshape(data,[1,1,1,length(data)]) - scaled_likelihood_means) ./ sqrt(scaled_likelihood_means))).^2,4);
    logLikelihoodTable = -abs(logLikelihoodTable(:));
else
    logLikelihoodTable = sum(reshape(data,[1,1,1,length(data)]) .* scaled_log_means - scaled_likelihood_means,4);
end

% prevent Inf from being the max log likelihood
logLikelihoodTable(logLikelihoodTable >= Inf) = -1e10;