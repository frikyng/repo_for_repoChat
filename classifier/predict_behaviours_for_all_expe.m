%% Use this script to run ML code on all extracted arboreal_scan_experiments

%% List extracted directories
top_export_folder   = 'C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans_final_new_zscored_compressed';
fold                = dir([top_export_folder,'/*-*-*_exp_*']);
fold                = fold([fold.isdir]);

%% Process
tic
out = {};
for fold_idx = numel(fold):-1:1
    fname = dir([fold(fold_idx).folder, '/', fold(fold_idx).name,'/20*-*-*_exp*.mat']);
    if ~isempty(fname)
        load([fname(1).folder, '/', fname(1).name]);
        out{fold_idx} = predict_behaviours(obj, true, 'linear','subtracted_peaks',{'encoder'}, '', 'optimize_hyper', true, 'optimization_method', 'manual');
    end
end
toc





