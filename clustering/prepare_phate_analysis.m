function [obj, source_signal, signal_indices] = prepare_phate_analysis(path_or_obj, use_hd_data, time_filter, type_of_trace)
    if nargin < 1 || isempty(path_or_obj)
        path_or_obj = ''; % i.e. current matlab folder. It must then be an extracted arboreal scan folder
    end
    if nargin < 2 || isempty(use_hd_data)
        use_hd_data = false; % HD is slower
    end
    if nargin < 3 || isempty(time_filter)
        time_filter = 0; % 
    end
    if nargin < 4 || isempty(type_of_trace)
        type_of_trace = 'raw'; % 
    end
    
    if ~any(contains(type_of_trace, {'raw','rescaled', 'subtracted'}))
        error('type_of_trace must contain "raw" or "rescaled" or "subtracted" ')
    end
    
    %% Get object if it needs loading/building
    if ischar(path_or_obj)
        [~,~,ext] = fileparts(path_or_obj);
        if strcmp(ext, '.mat')
            load('C:\Users\vanto\Documents\MATLAB\extracted_arboreal_scans 2\arboreal_scans_thin_mask.mat');
        else
            obj = arboreal_scan_experiment('', use_hd_data);
        end
    else
        obj = path_or_obj;
    end

    %% Rebuild object with HD data in it 
    %obj = arboreal_scan_experiment('C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans_2\extracted_arboreal_scans\2019-10-31_exp_1',use_hd_data)
    %obj = arboreal_scan_experiment('',use_hd_data)
    %obj.breakpoints = {'12-59-34'};

    %% ####### run once #######
    %% run once to get event time from low D
    %% Define wether to use HD data or LD data
    obj.use_hd_data = false;

    %% Quick processing to have proper event detection. This is required for filtering based on activity
    obj.rescaling_method = 'peaks_trials';
    obj.prepare_binning({'distance',50});
    obj.find_events();
    obj.rescale_traces();
    obj.set_median_traces();
    obj.compute_similarity();
    close all

    %% Define wether to use HD data or LD data, if using HD,
    if use_hd_data
        obj.use_hd_data = true;
        try
            obj.rescale_traces(); %% run twice .. to fix
        end
        obj.rescale_traces(); %% run once
    end
    close all


    %% Time filtering of traces if required. 
    %% ## ! ## This introduces temporal correlation
    obj.filter_win = [time_filter, 0];

    %% If using all voxels, remove bad_ROIs_list field because it is designed for full segments
    if obj.use_hd_data    
        obj.bad_ROI_list = [];
    end

    %% Select signal source (RAW or rescaled traces or shuffled)
    if contains(type_of_trace, 'raw')
        source_signal = obj.extracted_traces_conc;
    elseif contains(type_of_trace, 'rescaled')
        source_signal = obj.rescaled_traces;
    elseif contains(type_of_trace, 'subtracted')
        source_signal = obj.rescaled_traces - obj.global_median_raw;
    end

    if contains(type_of_trace, 'shuffle')
        %% randomize/permute/shuffle all data points in time series for each ROI separately (i.e. destroy all correlation)
        % see rand_shuffle_PHATE_loadings.m
    end
    
    %% Don't keep any Inf values, if any (only happens in HD case, occasionally)
    source_signal(isinf(source_signal)) = NaN;

    %% Flag ROIs that have NaN vaues at one point as they may mess up later computations
    bad_ROI_list                    = find(any(isnan(source_signal),1));
    bad_ROI_list(bad_ROI_list > obj.n_ROIs) = [];
    signal_indices                  = true(1, obj.n_ROIs); %% ROIs or voxels, depending on the data source
    signal_indices(bad_ROI_list)    = false;
    signal_indices                  = find(signal_indices);    
    
    source_signal = double(source_signal(:, 1:obj.n_ROIs)); 
    
    %% Filter out bad ROIs/voxels
    source_signal(:, bad_ROI_list) = NaN; 
    
end