function validateChannelPermutation(channelOrder, numChannels)
%VALIDATECHANNELPERMUTATION Require every raw channel exactly once.

arguments
    channelOrder (1,:) double {mustBeInteger,mustBePositive}
    numChannels (1,1) double {mustBeInteger,mustBePositive}
end

assert(numel(channelOrder) == numChannels && ...
    isequal(sort(channelOrder), 1:numChannels), ...
    'slap2:invalidChannelPermutation', ...
    'Channel order must contain every channel in 1:%d exactly once.', numChannels);
end
