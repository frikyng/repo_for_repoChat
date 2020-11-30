opengl software


[export_folder, data_folders_per_exp, setting_file_path] = prepare_meta_analysis('C:\Users\vanto\Desktop\new_meta_june_2020_v5');

%% Now, do the analysis expe-by-expe
for expe = 1:numel(data_folders_per_exp) 
   [all_traces_per_bin, all_ROI_ids_per_bin, blacklist, all_original_folders] = load_analyses(data_folders_per_exp{expe});
   
   data_folder = all_original_folders{1};   
   [~,data_folder,expe_folder, ~,~,~,tag] = parse_paths(data_folder, true);
   ref = imread([expe_folder, '/CONSENSUS_MEDIAN_STD_MOSAIC.png']);
   mask             = show_cell_mask(data_folder,'','',true, ref);
end