function processLoCo(dr)
disp(['Script started at: ' char(datetime("now"))])
if ~nargin
    dr = uigetdir; %neuron folder where scans are, not project folder
end

sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams = setParams('summarize_LoCo', sParams, true);

%summarize
sParams.batchSize = 100;
summarize_LoCo_ISOLATED(dr, sParams);
