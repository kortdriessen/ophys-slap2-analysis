directories_to_process = {
    %'/data/raw_mirror/miaplacidus/exp_1/loc_D/acq_1',
    '/data/raw_mirror/miaplacidus/exp_2/loc_D/acq_1',
    '/data/raw_mirror/bellatrix/exp_10/loc_J/acq_1',
};

for i = 1:numel(directories_to_process)
    processSLAP2_alignAndGetUserRois(directories_to_process{i});
end