function dataAnalysis = rawToAnalysisChannels(dataRaw, analysisChannelOrder, channelDimension)
%RAWTOANALYSISCHANNELS Permute canonical raw channels into model order.

arguments
    dataRaw
    analysisChannelOrder (1,:) double {mustBeInteger,mustBePositive}
    channelDimension (1,1) double {mustBeInteger,mustBePositive}
end

numChannels = numel(analysisChannelOrder);
validateChannelPermutation(analysisChannelOrder, numChannels);
assert(size(dataRaw, channelDimension) == numChannels, ...
    'slap2:channelCountMismatch', ...
    'Input channel dimension contains %d channels; expected %d.', ...
    size(dataRaw, channelDimension), numChannels);

subscripts = repmat({':'}, 1, max(ndims(dataRaw), channelDimension));
subscripts{channelDimension} = analysisChannelOrder;
dataAnalysis = dataRaw(subscripts{:});
end
