function contract = channelContract(numChannels, activityChannel)
%CHANNELCONTRACT Validate and describe raw and internal channel ordering.

arguments
    numChannels (1,1) double {mustBeInteger,mustBePositive}
    activityChannel (1,1) double {mustBeInteger,mustBePositive}
end

assert(activityChannel <= numChannels, 'slap2:invalidActivityChannel', ...
    'activityChannel (%d) must be in the range 1:%d.', activityChannel, numChannels);

contract.rawChannelOrder = 1:numChannels;
contract.analysisChannelOrder = [activityChannel, ...
    contract.rawChannelOrder(contract.rawChannelOrder ~= activityChannel)];

assert(numel(contract.analysisChannelOrder) == numChannels && ...
    isequal(sort(contract.analysisChannelOrder), contract.rawChannelOrder), ...
    'slap2:invalidAnalysisChannelOrder', ...
    'analysisChannelOrder must contain every physical raw channel exactly once.');
end
