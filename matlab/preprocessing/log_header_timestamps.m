%% WISC PROCESS TIME INFORMAITON - written for LINUX!!!
function log_header_timestamps(data_dir)
    % first get the metadatafile and save the fpgaTimeReference
    meta_list = dir([data_dir filesep '*.meta']);
    meta_path = fullfile(meta_list(1).folder, meta_list(1).name);
    meta = load(meta_path, "-mat");
    fpgaTimeRef = char(datetime(meta.fpgaTimeReference,'Format','MM-dd-yyyy HH:mm:ss.SSSSSS'));
    txt_name = fullfile(data_dir, 'fpgaTimeReference.txt');
    fid = fopen(txt_name, 'w');
    fprintf(fid, '%s', fpgaTimeRef);
    fclose(fid);
    
    % now get the start_timestamp of each data file
    dat_files = dir([data_dir filesep '*.dat']);
    
    for dmd_ix = 1:2
        isDMDix = cellfun(@(x)(~isempty(x)), strfind({dat_files.name}, ['DMD', num2str(dmd_ix)]));
        dat_files_dmd = dat_files(isDMDix);
        path = fullfile([dat_files_dmd(1).folder filesep dat_files_dmd(1).name]);
        hdat = slap2.Slap2DataFile(path);
        multi_data_files = hdat.hMultiDataFiles.hDataFiles;
        for md_file = 1:length(multi_data_files)
            header = multi_data_files(1, md_file).getLineHeader(1, 1);
            timestamp = header.timestamp_date;
            true_start = char(datetime(timestamp,'Format','MM-dd-yyyy HH:mm:ss.SSSSSS'));
            txt_id = ['datStartTimes', num2str(dmd_ix), '.txt'];
            txt_name = fullfile(data_dir, txt_id);
            fid = fopen(txt_name, 'a');
            fprintf(fid, '%s', true_start);
            fprintf(fid, '\r\n');
            fclose(fid);
        end
    end

end