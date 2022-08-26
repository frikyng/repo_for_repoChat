TOP_FOLDER = 'D:\Curated Data\';
[~, expe_folders] = list_expe_folders(TOP_FOLDER);

for el = fliplr(1:numel(expe_folders))
    [~, df] = list_data_folders(expe_folders{el});
    data_folder = df{2};
    h = load_header(data_folder);
    
    try
    a_s = arboreal_scan(data_folder);
    a_s.prepare_tree([TOP_FOLDER,default_batch_params_filename]);
    a_s.plot_dist_tree();  
    catch err
        data_folder
        err
        warning('error detected, an code excution pas paused. Hit F5 to resume')
        keyboard()
    end
end


meta_batch_process_ribbon_scan
obj = arboreal_scan_experiment(expe_folders{end});
