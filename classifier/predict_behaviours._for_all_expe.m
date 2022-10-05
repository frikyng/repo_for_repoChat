%% Use this script to build all arboreal_scan_experiments objects 
%% automatically, rescale peaks, extract events and save the updated object.

%% List extracted directories
top_export_folder   = 'C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans_thin_mask\extracted_arboreal_scans';
fold                = dir([top_export_folder,'/*-*-*_exp_*']);
fold                = fold([fold.isdir]);

%% Set analysis settings
process             = true;
full_process = true; % if false, we use the custom list of preocessing set (see "else" condition in the loop)

%% Process
tic
errors = {};
for fold_idx = numel(fold):-1:1
    fname = dir([fold(fold_idx).folder, '/', fold(fold_idx).name,'/20*-*-*_exp*.mat']);
    if ~isempty(fname)
        load([fname(1).folder, '/', fname(1).name]);        
        out{fold_idx} = classify_behaviours(obj);
    end
end
toc





