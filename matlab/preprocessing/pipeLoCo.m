function pipeLoCo(dr)
disp(['Script started at: ' char(datetime("now"))])

params.microscope = { '''SLAP2''' , '''bergamo'''};          tooltips.scope = 'SLAP2 or bergamo';
params.includeIntegrationROIs = false; tooltips.includeIntegrationROIs = 'Use integration ROIs for trace extraction?';
params.sigma_px = 1.33;          tooltips.sigma_px = 'Estimated radius of the PSF (gaussian sigma)';
params.nmfIter = 2;              tooltips.nmfIter = 'number of iterations of NMF refinement';
params.dXY = 3;                  tooltips.dXY = 'how large sources can be (radius), pixels';
params.lambda = [];              tooltips.lambda = 'regularizer; roughly the single-photon amplitude. Leave empty to use default/estimate from data.';
params.denoiseWindow_s = 0.2;   tooltips.denoiseWindow_s= 'the timescale on which signals can be smoothed when denoising, seconds';
params.baselineWindow_Glu_s = 4; tooltips.baselineWindow_Glu_s= 'timescale for calculating F0 in glutamate channel, seconds';
params.baselineWindow_Ca_s = 4;  tooltips.baselineWindow_Ca_s= 'timescale for calculating F0 in calcium channel, seconds';
params.activityChannel = 2;      tooltips.activityChannel = 'the channel of the original tiff image that contains the glutamate signal';
params.tau_s = 0.03;             tooltips.tau_s = 'decay time constant of glutamate signal';
params.tau2_s = 0.15;            tooltips.tau2_s = 'decay time constant of 2nd channel signal at synapses (usually spine calcium)';
params.maxSynapseDensity = 0.01; tooltips.maxSynapseDensity = 'maximum synapses per pixel';
params.poissBasedStdIM = 0;      tooltips.poissBasedStdIM = 'use Poisson model to estimate stdIM';
params.VIF = 1;                  tooltips.VIF = 'variance inflation factor for stdIM estimate (Poiss-based only)';
params.peakth = 3.5;             tooltips.peakth = 'peak identification threshold (actIM z-score)';
params.peakFuncOpt = 2;             tooltips.peakFuncOpt = 'peak fitting function (1=gaussian, 2=binned gaussian)';
params.actImHeteroscedasticNoise = 1;             tooltips.actImHeteroscedasticNoise = 'noise in actIM modeled as heteroscedastic';
params.peakBufferSize = 0;             tooltips.peakBufferSize = 'number of pixels to mask around each peak';
params.nParallelWorkers = 12;    tooltips.nWorkers = 'number of parallel workers';
params.drawUserRois = true;     tooltips.drawUserRois = 'pop up a GUI to annotate user ROIs?';  
params.motionThresh = 2.5;       tooltips.motionThresh = 'decrease this to be more stringent on motion correction when censoring frames';
params.analyzeHz = 200;          tooltips.analyzeHz = 'frame rate used for analysis (SLAP2 only)';
params.nanThresh = 0.33;         tooltips.nanThresh = 'Max fraction of samples that can be NaN for including a pixel in analysis';
params.discardInitial_s = 0;     tooltips.discardInitial_s = 'time in seconds to remove from analysis at the start of each trial, to accound for warmup';
params.operator = 'KD';       tooltips.operator = 'person running the analysis';
params.makeJSON = false;             tooltips.makeJSON = 'run python script to create processing.json';

%summarize
sParams.batchSize = 100;
summarize_LoCo_ISOLATED(dr, sParams);