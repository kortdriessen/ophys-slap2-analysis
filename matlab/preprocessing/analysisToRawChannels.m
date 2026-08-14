function dataRaw = analysisToRawChannels(dataAnalysis, analysisChannelOrder, channelDimension)
%ANALYSISTORAWCHANNELS Map model-ordered channels to physical raw order.

arguments
    dataAnalysis
    analysisChannelOrder (1,:) double {mustBeInteger,mustBePositive}
    channelDimension (1,1) double {mustBeInteger,mustBePositive}
end

numChannels = numel(analysisChannelOrder);
validateChannelPermutation(analysisChannelOrder, numChannels);
assert(size(dataAnalysis, channelDimension) == numChannels, ...
    'slap2:channelCountMismatch', ...
    'Input channel dimension contains %d channels; expected %d.', ...
    size(dataAnalysis, channelDimension), numChannels);

dataRaw = zeros(size(dataAnalysis), 'like', dataAnalysis);
subscripts = repmat({':'}, 1, max(ndims(dataAnalysis), channelDimension));
subscripts{channelDimension} = analysisChannelOrder;
dataRaw(subscripts{:}) = dataAnalysis;
end
