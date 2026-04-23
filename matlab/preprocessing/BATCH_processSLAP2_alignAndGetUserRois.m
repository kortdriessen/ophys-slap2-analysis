directories_to_process = {
    '/data/regora/exp_1/loc_K/acq_1',

};

for i = 1:numel(directories_to_process)
    processSLAP2_alignAndGetUserRois(directories_to_process{i});
end