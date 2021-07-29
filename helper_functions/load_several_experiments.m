function data = load_several_experiments(ROIs, data_folders, use_mask)
    if nargin < 3 || isempty(use_mask) || (islogical(use_mask) && use_mask)
        mask_name = 'auto_mask.mat';
        mask_method = 'user_defined';
    elseif (ischar(use_mask) || isstring(use_mask)) && isfile(use_mask)
        mask_name = use_mask;
        mask_method = 'user_defined';
    elseif ischar(use_mask) && ~isfile(use_mask)
        error_box('indicated mask file doesnt exist anymore')
    else
        mask_name = '';
        mask_method = 'user_defined';
    end
        
    p = analysis_params(   'data_type'      , 'raw',...
                           'tracing_source' , 'swc'         , ...
                           'signal_channel' , 2             ,...
                           'registration'   , true          ,...
                           'registration_reload','auto_offsets.mat',...
                           'smoothing'      , [10,0]        ,...
                           'mask_method'    , mask_method   ,...
                           'mask_value'     , mask_name     ,...
                           'compression'    , '2D'          ,...
                           'rendering_mode' , 'individual'      ,...
                           'rendering'      , false         ,...
                           'ROIs'           , ROIs); 
    

%     p.mask_method = 'none'
%     p.mask_value = '';
    if ~isdir(data_folders{1})
        error_box(['You are trying to reload the original data, but the indicated path ( ', data_folders{1}, ' does not exist anymore'],1)
        return
    end
    br = get_branch_id_from_ROI(data_folders{1}, 0, ROIs);
    ROIs = get_branch_id_from_ROI(data_folders{1}, br, 0);
    p.ROIs = ROIs;


    if numel(data_folders) == 1
        p.rendering_mode = 'mosaic';
        p.rendering = true;
        load_ribbon_scan(p, 'data_folder', data_folders{1})
    else    
        all_data = {};
        all_durations = {};
        fprintf('Please wait while loading all the recordings for the selected ROI\n')
        [~, header] = import_ROIs(p, 'data_folder', data_folders{1}, 'repeats', 1);
        parfor f_idx = 1:numel(data_folders)
            f                   = data_folders{f_idx};
            info                = get_recordings_info(f);
            h                   = load_header(f);
            info                = merge_params_obj(info, h); % this update some values for import ROI
            current_data        = [];
            current_duration    = 0;
            for trial = info.repeats
                    %[results, header] = import_ROIs(p, 'data_folder', f{1}, 'repeats', trial);
                    results     = [];
                for roi = ROIs % if multiple ROIs from a branch
                    fname       = dir([info.data_folder '/RibbonScan_ROI_',sprintf('%04d',roi),'_repeat_',sprintf('%04d',trial),'*.mat']); %load list of all the files (excluded concatenation and averages)
                    fname       = parse_paths([fname.folder,'/',fname.name]);
                    results     = cat(1, results, import_ROI(fname, merge_params_obj(p, info)));                
                        %[results, header] = import_ROIs(p, 'data_folder', f{1}, 'repeats', trial);
                        %[temp, ~, ~, params] = load_experiment(p, 'data_folder', f{1});
                end
                current_data    = cat(4, current_data, results);
                current_duration= current_duration + info.estimated_trial_duration;        
            end
            if iscell(current_data)
                current_data = current_data{1}
            end
            all_data{f_idx} = current_data;
            all_durations{f_idx} = current_duration;
        end
        
        
        data = cat(4, all_data{:});
        data = data(:,:,:,:,2);
        header.ROIs = ROIs;
        header.trial_duration = sum(cell2mat(all_durations));
        header.duration = sum(cell2mat(all_durations));
        header.time_range = [0, header.estimated_total_duration];
        header.timepoints = size(data, 4);
        header.points_per_s = header.timepoints/header.duration;
        header.data_folder = data_folders{1};
        header.timescale = linspace(0,header.duration,header.timepoints);
        header.data_type = 'concatenated';
        header.original_repeats = 'all';
        
        Ribbon_viewer(data, header);
    end
    
    
    for f = data_folders
    	cleanup_processed_data(f{1});
    end
end