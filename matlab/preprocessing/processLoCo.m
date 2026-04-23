function processLoCo(dr)
disp(['Script started at: ' char(datetime("now"))])
if ~nargin
    %dr = '/run/user/1329238735/gvfs/smb-share:server=tononi-nas,share=slap_mi/slap_mi/data/alnair/exp_4/loc_M/acq_1';
    dr = uigetdir; %neuron folder where scans are, not project folder
end

sParams.microscope = "SLAP2";
sParams.drawUserRois = true;
sParams = setParams('summarize_LoCo', sParams, true);

%summarize
sParams.batchSize = 100;
summarize_LoCo_ISOLATED(dr, sParams);

