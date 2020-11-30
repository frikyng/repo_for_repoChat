function [fixed_tree, soma_location, ROIs_list, listing, batch_params] = rebuild_tree(source, setting_file_path, show_distance_tree)
    if nargin < 3 || isempty(show_distance_tree)
        show_distance_tree = false;
    end
    if ischar(setting_file_path)
        rebuilt_expe_name = fileparts(source.folder);
        [~, rebuilt_expe_name] = fileparts(rebuilt_expe_name);
        rebuilt_expe_name = strrep(rebuilt_expe_name,'_exp_','/experiment_');
        batch_params = get_paths_from_settings(setting_file_path,'','',rebuilt_expe_name);
        correct_data_folder = contains(batch_params.data_folder, source.name(end-7:end));
        batch_params    = structfun(@(x) x(correct_data_folder), batch_params, 'UniformOutput', false);
    else
        % batch_params = structfun(@(x) x(1), batch_params,'UniformOutput', false); if more than 1
        batch_params =setting_file_path; %you need to pass 1 batch prms
    end

    [fixed_tree, ROI_distance_fig, soma_location] = fix_tree(batch_params.data_folder, batch_params.primary_branches, batch_params.manual_reconnections, batch_params.excluded_branches, batch_params.pia_soma, show_distance_tree);

    %[fixed_tree, ROI_distance_fig, soma_location] = fix_tree(batch_params.data_folder{1}, batch_params.primary_branches{1}, batch_params.manual_reconnections{1}, batch_params.excluded_branches{1}, show_distance_tree);

    %% 3 step process to get ROIs that are not excluded
    branch_idx = get_branch_id_from_tracing_source(struct('data_folder',batch_params.data_folder,'tracing_source','swc','ROIs',0,'branch_ids',0));
    branch_idx = branch_idx(~ismember(branch_idx, batch_params.excluded_branches)); % non exluded branches
    ROIs_list = get_branch_id_from_ROI(batch_params.data_folder, branch_idx, 0);

    [~, ~, listing] = get_branch_id_from_ROI(batch_params.data_folder, 0, ROIs_list, ROIs_list);
    listing = vertcat(listing{:});
end