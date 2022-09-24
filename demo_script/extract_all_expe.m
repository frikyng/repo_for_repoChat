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
errors = {};
for fold_idx = numel(fold):-1:1
	try
        if numel(dir([fold(fold_idx).folder, '/', fold(fold_idx).name,'/'])) == 2
            errors{fold_idx} = 'empty';
        else
            current_fold = [fold(fold_idx).folder, '/', fold(fold_idx).name,'/'];
            obj = arboreal_scan_experiment(current_fold, false);
            
            %% Full processing
            if process && full_process
            	obj.process({'depth',50},[ceil(1/nanmedian(expe.timescale.sr)), 0]); %% Add grouping settings here
            elseif process
            %% Manual step by step rocessing
                obj.rescaling_method = 'peaks_trials';
                obj.prepare_binning({'distance',50});
                obj.find_events();
                obj.rescale_traces();
                obj.set_median_traces();
                obj.compute_similarity();
            end
            
            close all
            obj.save(true)
        end
    catch err
        errors{fold_idx} = err;
    end
end