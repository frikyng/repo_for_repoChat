

top_export_folder = 'C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans\extracted_arboreal_scans';
fold = dir([top_export_folder,'/*-*-*_exp_*']);
fold = fold([fold.isdir]);
errors = {};
to_redo = find(~cellfun(@isempty, errors));
for idx = numel(fold):-1:1 %71 --> double cell, 16 failing
    try
        %         expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
        %%         if contains(errors{idx}.message, 'value')
        %             [~,~,expe_f] = parse_paths(expe.ref.data_folder, true);
        %             meta_batch_extract_arboreal_scan(expe_f,'', '','C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans');
        %         end
        
        file = dir([fold(idx).folder,'/',fold(idx).name,'/*-*-*_exp_*.mat']);
        expe = load([file(1).folder,'/',file(1).name]);
        expe = expe.obj;
        expe.rendering = 1;
        expe.detrend     = 0;  
%         expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
        expe.auto_save_analysis = 0;
        expe.auto_save_figures = 0;
        expe.bad_ROI_list = 'unset'
        %expe.filter_win = [20,0]
        expe.process({'distance',100},[ceil(1/nanmedian(expe.timescale.sr)), 0]); %% Add grouping settings here
        errors{idx} = [];
    catch err
        errors{idx} = err;
        1
    end
end