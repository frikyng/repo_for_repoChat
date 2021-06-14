function [extracted_traces, all_original_folders, all_sr, all_infos] = load_analyses(data_folders_current_exp, tracing_source)
    if nargin < 3 || isempty(tracing_source)
        tracing_source = 'swc';
    end

    extracted_traces        = {}; % N recordings cell array of M analysis bins cell array of [Timepoints x ROI] matrices
    all_sr                  = {}; % N FLOAT of sampling rates
    all_infos               = {}; % N recordings cell array of analysis_params, without data

    %% Load extracted data
    blacklist   = [];
    for data_f_idx = 1:numel(data_folders_current_exp)
        paths = parse_paths([data_folders_current_exp(data_f_idx).folder, '/', data_folders_current_exp(data_f_idx).name]);
        source              = dir([paths,'/**/*.mat']);
        if ~isempty(source)
            extracted_data                  = load([source.folder, '/', source.name]);
            extracted_data                  = extracted_data.obj;

%             %% Reselect ROI subset
%             if strcmp(tracing_source, 'swc')
%                params = analysis_params('source'          ,extracted_data.data_folder,...
%                                         'tracing_source'  , 'swc'         );    
%                [~, ROIs] = get_branch_id_from_tracing_source(params);
%             end

            extracted_traces{data_f_idx}                = extracted_data.data;
            error_box('the following fieldnames have changed. check and clean before using')
            extracted_data.data                         = [];
            extracted_data.analysis_params.full_data    = [];
            extracted_data.analysis_params.simple_data  = [];
            
            all_infos{data_f_idx}           = extracted_data;
            all_sr{data_f_idx}              = diff(extracted_data.analysis_params.timescale(1:2));  
        else
            blacklist                       = [blacklist, data_f_idx];
        end
    end
    

    %% Get rid of any bad file
    if ~isempty(blacklist)
        error_box('Curation and/or data extraction not clean. Consider reviewing this folder')
    end
    
    extracted_traces      = extracted_traces(~cellfun(@isempty , extracted_traces));
    all_sr                  = all_sr(~cellfun(@isempty , all_sr));    
    all_original_folders    = arrayfun(@(x) x.data_folder, [all_infos{:}], 'UniformOutput' ,false);
end

