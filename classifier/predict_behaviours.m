%% Classify behaviour for a given experiment, for all behaviours
%% use expe 07-26 to fix gaps in behaviour

%% Cell 2019-09-17_exp_1 anticorelated to running

function [out, data, ROI_groups, meanvalue] = predict_behaviours(obj, use_classifier, method, type_of_trace, behaviour_list, ROI_groups, n_iter, rand_ROI_groups, varargin)
    if nargin < 1 || isempty(obj)
        obj = ''; 
    end
    if nargin < 2 || isempty(use_classifier)
        use_classifier = true; % if false, use regression learner
    end
    if nargin < 3 || isempty(method)
        method = 'svm'; % i.e. fitcsvm vs fitclinear etc
    end
    if nargin < 4 || isempty(type_of_trace)
        type_of_trace   = 'subtracted_peaks'; % ['subtracted' OR 'rescaled' OR 'raw'] AND ['peaks' or '']. eg 'subtracted_peaks' , or 'raw'
    end
    build_beh           = true;
    if nargin < 5 || isempty(behaviour_list)
        % pass      
    elseif isstruct(behaviour_list)
        raw_behaviours = behaviour_list.raw_behaviours;
        beh_thr        = behaviour_list.beh_thr;
        formatted_behaviour_list = behaviour_list.formatted_behaviour_list;
        behaviour_list = behaviour_list.original_behaviour_list;
        build_beh      = false;
    end
    single_matrix_input = false;
    if nargin < 6 || isempty(ROI_groups)
        ROI_groups   = num2cell(obj.ref.indices.valid_swc_rois);
        single_matrix_input = true;
    elseif ~iscell(ROI_groups)
        ROI_groups   = num2cell(ROI_groups);    % does this force groups to be cell array regardless of input type?
        single_matrix_input = true;
    end
    if iscolumn(ROI_groups)
        ROI_groups = ROI_groups';
    end
    if nargin < 7 || isempty(n_iter)
        n_iter   = 1; 
    end
    if nargin < 8 || isempty(rand_ROI_groups)
        rand_ROI_groups   = false; % false is no randomization, true is randperm (this is to compare info content in a group of (e.g. PHATE discovered) ROIs with randomly selected ROIs (from anywhere from the tree), esuring the same # of ROIs in each)
    end
    if nargin < 9 || isempty(varargin)
        ml_parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        ml_parameters = DEFAULT_CLASSIFIER_OPTION(varargin);
    end
    
    use_hd_data             = false;
    time_filter             = 0;
    
    %% Make sure preprocessing was done correctly
    rendering       = obj.rendering;
    obj.rendering   = false;

    [obj, source_signal, ~, timepoints, lag]    = prepare_phate_analysis(obj, use_hd_data, time_filter, type_of_trace);
    
    beh_timepoints                              = timepoints + lag; % if lag == 0 then points are the same
    timepoints(beh_timepoints < 1)              = [];
    beh_timepoints(beh_timepoints < 1)          = [];
    timepoints(beh_timepoints > obj.tp)         = [];
    beh_timepoints(beh_timepoints > obj.tp)     = [];
    
    % holder = source_signal(:,1); figure;plot(source_signal(:,1));hold on; scatter(timepoints, holder(timepoints),'r','filled')
    % figure();hist(reshape(source_signal(timepoints,:),[],1),100); % TOMMY uncomment to see distribution of predictors
    obj.rendering   = rendering;
    
    %% Get all behaviours
    if ~build_beh
        %% If we already had the extracted behaviours, at this stage we just update the formatted name 
    	behaviour_list = formatted_behaviour_list;
    else
        [~, raw_behaviours, beh_thr, behaviour_list] = prepare_behaviour_data(obj, behaviour_list);
    end
    
    %% Filter timepoints
    processed_behaviours = raw_behaviours(:, beh_timepoints);
    if use_classifier
        processed_behaviours = logical(processed_behaviours > beh_thr');
    end

    %% Get the signal for the selected timepoints and ROIs
    %Valid_ROIs      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); % remove excluded branches AND bad_ROIs based on correlation)
    invalid_ROIs_logical    = ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list);
    invalid_ROIs            = obj.ref.indices.valid_swc_rois(invalid_ROIs_logical); % remove excluded branches AND bad_ROIs based on correlation)
    source_signal(:, invalid_ROIs) = NaN;
    data                    = NaN(numel(ROI_groups), numel(timepoints));
    for gp_idx = 1:numel(ROI_groups)
        ROI_groups{gp_idx}(ismember(ROI_groups{gp_idx}, obj.bad_ROI_list)) = [];
        data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups{gp_idx}),2)';        
    end

    fully_invalid_group = cellfun(@isempty, ROI_groups)' | all(isnan(data),2);
    ROI_groups(fully_invalid_group) = [];
    data(fully_invalid_group, :) = [];

    %     %% Correlation of the different ROIs with each other
    %     figure();imagesc(corr(data'))
    % 
    %     %% Correlation of the Ca2+ values per event
    %     figure();imagesc(corr(data))


    %     %% Correlation of the different behaviours with each other
    %     figure();imagesc(corr(processed_behaviours'))
    % 
    %     %% Correlation of the behaviours values per event
    %     figure();imagesc(corr(processed_behaviours))

    %     Soma_ROIs   = find(ismember(1:size(data, 1), obj.ref.indices.somatic_ROIs)); % somatic ROIs, but ignoring the NaNs

    
    %% FYI find(all(isnan(source_signal))) should be empty, or you have observation with only NaN
    
    All_ROIs    = 1:size(data, 1);

    if ~rand_ROI_groups
        for iter = 1:n_iter               
            out{iter}               = train_and_test(data, processed_behaviours, timepoints, All_ROIs, method, behaviour_list, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            out{iter}.used_ROIs     = ROI_groups;
        end
    else  
        %% we build alternative randomized groups. 
        % 0/false does no randomization
        % 1 randomize groups using all valid ROIs
        % -1 randomize groups using all valid ROIs, excluding the ROIs listed in the groups. Groups are sized matched
        if rand_ROI_groups == 1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        elseif rand_ROI_groups == -1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois((~ismember(obj.ref.indices.valid_swc_rois, horzcat(ROI_groups{:}))) & (~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list))); %list of valid ROIs excluding ROI groups
        end
        % when you pass all ROIs (so all cells in ROI_groups size == 1), randomization for leftover ROIs makes no sense so we just randomize ROIs
        if isempty(rand_ROI_pool)
            disp('WARNING : No ROI left for randomization')
            rand_ROI_groups     = 1;
            rand_ROI_pool       = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        end
        for iter = 1:n_iter
            try
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(rand_ROI_pool(randperm(numel(rand_ROI_pool),numel(ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) rand_ROI_pool(randperm(numel(rand_ROI_pool),x)), cellfun(@numel, ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            catch
                disp('WARNING : Not enough non-used ROIs available. Using all ROIs instead for randomization')
                ok = find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list));
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(ok(randperm(numel(ok),numel(ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) ok(randperm(numel(ok),x)), cellfun(@numel, ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            end
            for gp_idx = 1:numel(ROI_groups_rdm)
                data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups_rdm{gp_idx}),2)';        
            end   
            out{iter}               = train_and_test(data, processed_behaviours, timepoints, 1:numel(ROI_groups_rdm), method, behaviour_list, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            out{iter}.used_ROIs     = ROI_groups_rdm;
        end
    end

    [meanvalue,~, fig_handle] = bar_chart(out, 'beh_type','','','',ml_parameters.rendering);
    if ml_parameters.rendering
        title(ml_parameters.title)
    end
    if ml_parameters.savefig
        if islogical(ml_parameters.savefig)
            save_myfig(fig_handle,ml_parameters.title,{'png','pdf'})
        elseif ischar(ml_parameters.savefig)
            if contains(ml_parameters.savefig, '.mat')
                ml_parameters.savefig = parse_paths(fileparts(ml_parameters.savefig));
            end
            save_myfig(fig_handle,[ml_parameters.savefig, ml_parameters.title],{'png','pdf'})
        end
    end    
end