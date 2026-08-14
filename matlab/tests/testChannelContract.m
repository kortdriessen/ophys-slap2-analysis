function tests = testChannelContract
%TESTCHANNELCONTRACT Exact synthetic tests for raw and analysis channel axes.
tests = functiontests(localfunctions);
end

function testSingleChannelIdentity(testCase)
contract = channelContract(1, 1);
verifyEqual(testCase, contract.rawChannelOrder, 1);
verifyEqual(testCase, contract.analysisChannelOrder, 1);
tagged = reshape(1:6, [2 3 1]);
verifyEqual(testCase, rawToAnalysisChannels(tagged, 1, 3), tagged);
verifyEqual(testCase, analysisToRawChannels(tagged, 1, 3), tagged);
end

function testActivityChannelOneIdentity(testCase)
contract = channelContract(2, 1);
verifyEqual(testCase, contract.analysisChannelOrder, [1 2]);
end

function testActivityChannelTwoPermutation(testCase)
contract = channelContract(2, 2);
verifyEqual(testCase, contract.analysisChannelOrder, [2 1]);
raw = zeros(2, 3, 2, 4);
raw(:,:,1,:) = 11;
raw(:,:,2,:) = 22;
analysis = rawToAnalysisChannels(raw, contract.analysisChannelOrder, 3);
verifyEqual(testCase, analysis(:,:,1,:), raw(:,:,2,:));
verifyEqual(testCase, analysis(:,:,2,:), raw(:,:,1,:));
verifyEqual(testCase, ...
    analysisToRawChannels(analysis, contract.analysisChannelOrder, 3), raw);
end

function testOutputRemapping(testCase)
analysis = zeros(3, 5, 2);
analysis(:,:,1) = 2002;
analysis(:,:,2) = 1001;
raw = analysisToRawChannels(analysis, [2 1], 3);
verifyEqual(testCase, raw(:,:,1), repmat(1001, 3, 5));
verifyEqual(testCase, raw(:,:,2), repmat(2002, 3, 5));
end

function testInvalidOrdersRejected(testCase)
verifyError(testCase, @()validateChannelPermutation([1 1], 2), ...
    'slap2:invalidChannelPermutation');
verifyError(testCase, @()validateChannelPermutation(1, 2), ...
    'slap2:invalidChannelPermutation');
verifyError(testCase, @()validateChannelPermutation([1 3], 2), ...
    'slap2:invalidChannelPermutation');
verifyError(testCase, @()channelContract(2, 3), ...
    'slap2:invalidActivityChannel');
end

function testDistinctBaselineWindows(testCase)
params = minimalExtractParams();
params.baselineWindow_Glu_s = 2;
params.baselineWindow_Ca_s = 7;
result = setParamsExtractTrial(params);
verifyEqual(testCase, result.baselineWindow_Glu_samps, 400);
verifyEqual(testCase, result.baselineWindow_Ca_samps, 1400);
verifyEqual(testCase, result.analysisBaselineWindows_samps, [400 1400]);

gluChanged = params;
gluChanged.baselineWindow_Glu_s = 3;
gluResult = setParamsExtractTrial(gluChanged);
verifyEqual(testCase, gluResult.baselineWindow_Glu_samps, 600);
verifyEqual(testCase, gluResult.baselineWindow_Ca_samps, 1400);

caChanged = params;
caChanged.baselineWindow_Ca_s = 9;
caResult = setParamsExtractTrial(caChanged);
verifyEqual(testCase, caResult.baselineWindow_Glu_samps, 400);
verifyEqual(testCase, caResult.baselineWindow_Ca_samps, 1800);
end

function testWisconsinSchema(testCase)
schema = buildWiscChannelSchema(2, 2, 0.03, 0.15);
verifyEqual(testCase, schema.rawChannelNumbers, [1 2]);
verifyEqual(testCase, schema.analysisChannelOrder, [2 1]);
verifyEqual(testCase, schema.kernelTauByRawChannel_s, [0.15 0.03]);
verifyEqual(testCase, schema.savedTraceOrder, 'raw');
verifyError(testCase, @()buildWiscChannelSchema(2, 1, 0.03, 0.15), ...
    'slap2:wiscSchemaActivityChannel');
end

function params = minimalExtractParams()
params.sigma_px = 1.33;
params.selRadius = 3;
params.tau_full = 6;
params.tau2_s = 0.15;
params.analyzeHz = 200;
params.baselineWindow_Glu_s = 4;
params.baselineWindow_Ca_s = 4;
params.denoiseWindow_s = 0.2;
params.numChannels = 2;
params.microscope = 'SLAP2';
params.lambda = 25;
end
