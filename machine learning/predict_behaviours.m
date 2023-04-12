%% Classify behaviour for a given experiment, for all behaviours
%% use expe 07-26 to fix gaps in behaviour

%% Cell 2019-09-17_exp_1 anticorelated to running

% * obj is a preprocessed arboreal_scan_experiment object. You need to have
% extracted peaks and resclaed traces, which is sensitive to factors such
% as detrending and smoothing
% * predictors_type defines the source data used for the predictor, as
% computed in the prepare_analysis() function. For example, you can define
% - if you want to use the whole trace or just the peaks, or any other
%   subset of timepoints
% - If you want to introduce a lag with the observations
% - The rescaled data or the median subtracted traces
% * observations is typically a single behaviour or a list of behvaiour
% names. Timepoints are the same than the one set in predictors_type,
% although it can have a temporal lag.
% * predictor_ROI_groups is the list of predictors to use. By default, all
% ROIs of the trees are used independently. If you pass an array, all ROIs 
% of the array are used indpendently. If you pass a cell array, each cell 
% form one predictor. If a cell contains more than one value, all ROis 
% listed within the cell are averaged.
% * varargin are either a machine_learning_params or set of valid
% name/argument pairs as described in machine_learning_params.m These
% options control the machine learning process

%% See also, prepare_analysis, bar_chart, train_and_test

function [results, data, predictor_ROI_groups, meanvalue, stats, values] = predict_behaviours(obj, predictors_type, observations, predictor_ROI_groups, varargin)
    if nargin < 1 || isempty(obj)
        obj = ''; 
    end
    if nargin < 2 || isempty(predictors_type)
        predictors_type   = 'subtracted_peaks'; % ['subtracted' OR 'rescaled' OR 'raw'] AND ['peaks' or '']. eg 'subtracted_peaks' , or 'raw'
    end
    build_beh           = true;
    if nargin < 3 || isempty(observations)
        % pass      
    elseif isstruct(observations)
        raw_behaviours = observations.raw_behaviours;
        beh_thr        = observations.beh_thr;
        formatted_behaviour_list = observations.formatted_behaviour_list;
        observations = observations.original_behaviour_list;
        build_beh      = false;
    end
    single_matrix_input = false;
    if nargin < 4 || isempty(predictor_ROI_groups)
        predictor_ROI_groups   = num2cell(obj.ref.indices.valid_swc_rois);
        single_matrix_input = true;
    elseif ~iscell(predictor_ROI_groups)
        predictor_ROI_groups   = num2cell(predictor_ROI_groups);    % does this force groups to be cell array regardless of input type?
        single_matrix_input = true;
    end
    if iscolumn(predictor_ROI_groups)
        predictor_ROI_groups = predictor_ROI_groups';
    end
    if nargin < 5 || isempty(varargin)
        ml_parameters = machine_learning_params;
    else
        ml_parameters = machine_learning_params(varargin);
    end

    
    use_hd_data             = false;
    
    %% Make sure preprocessing was done correctly
    rendering       = obj.rendering;
    obj.rendering   = false;

    [obj, source_signal, ~, timepoints, lag]    = prepare_analysis(obj, use_hd_data, [], predictors_type);
    
    beh_timepoints                              = timepoints + lag; % if lag == 0 then points are the same
    timepoints(beh_timepoints < 1)              = [];
    beh_timepoints(beh_timepoints < 1)          = [];
    timepoints(beh_timepoints > obj.tp)         = [];
    beh_timepoints(beh_timepoints > obj.tp)     = [];
    
    % holder = source_signal(:,1); figure;plot(source_signal(:,1));hold on; scatter(timepoints, holder(timepoints),'r','filled')
    % figure();hist(reshape(source_signal(timepoints,:),[],1),100); % TOMMY uncomment to see distribution of predictors
    obj.rendering   = rendering;
    
    %% If required add shufffing
    if ml_parameters.add_shuffle
        shuffled_beh    = strcat(observations, '_shuffled');
        observations  = reshape([observations; shuffled_beh],[],1)';
    end
    
    %% Get all behaviours
    if ~build_beh
        %% If we already had the extracted behaviours, at this stage we just update the formatted name 
    	observations = formatted_behaviour_list;
    else
        [~, raw_behaviours, beh_thr, observations] = prepare_behaviour_data(obj, observations);
    end
    
    %% Filter timepoints
    processed_behaviours = raw_behaviours(:, beh_timepoints);
    if ml_parameters.use_classifier
        processed_behaviours = logical(abs(processed_behaviours) > beh_thr');
    end
    
    %% Find groups that are exact duplicates because that would make 2 times the same predictor
    for el = 1:(numel(predictor_ROI_groups)-1)
        for el2 = (el+1):numel(predictor_ROI_groups)
            if numel(predictor_ROI_groups{el}) == numel(predictor_ROI_groups{el2}) && all(predictor_ROI_groups{el} == predictor_ROI_groups{el2})
                predictor_ROI_groups{el2} = [];
            end
        end
    end

    %% Get the signal for the selected timepoints and ROIs
    invalid_ROIs_logical    = ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list);
    invalid_ROIs            = obj.ref.indices.valid_swc_rois(invalid_ROIs_logical); % remove excluded branches AND bad_ROIs based on correlation)
    source_signal(:, invalid_ROIs) = NaN;
    data                    = NaN(numel(predictor_ROI_groups), numel(timepoints));
    for gp_idx = 1:numel(predictor_ROI_groups)
        predictor_ROI_groups{gp_idx}(ismember(predictor_ROI_groups{gp_idx}, obj.bad_ROI_list)) = [];
        data(gp_idx, :) =  nanmean(source_signal(timepoints, predictor_ROI_groups{gp_idx}),2)';        
    end
    
    %figure();imagesc(tril(corr([data',processed_behaviours'])',-1))

    fully_invalid_group = cellfun(@isempty, predictor_ROI_groups)' | all(isnan(data),2);
    predictor_ROI_groups(fully_invalid_group) = [];
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

    if ~ml_parameters.randomize_ROIs
        for iter = 1:ml_parameters.N_iter
            fprintf(['Iteration : ',num2str(iter),'\n'])
            results{iter}               = train_and_test(data, processed_behaviours, timepoints, All_ROIs, observations, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            results{iter}.used_ROIs     = predictor_ROI_groups;
        end
    else  
        %% we build alternative randomized groups. 
        % 0/false does no randomization
        % 1 randomize groups using all valid ROIs
        % -1 randomize groups using all valid ROIs, excluding the ROIs listed in the groups. Groups are sized matched
        if ml_parameters.randomize_ROIs == 1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        elseif ml_parameters.randomize_ROIs == -1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois((~ismember(obj.ref.indices.valid_swc_rois, horzcat(predictor_ROI_groups{:}))) & (~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list))); %list of valid ROIs excluding ROI groups
        end
        % when you pass all ROIs (so all cells in ROI_groups size == 1), randomization for leftover ROIs makes no sense so we just randomize ROIs
        if isempty(rand_ROI_pool)
            disp('WARNING : No ROI left for randomization')
            ml_parameters.randomize_ROIs        = 1;
            rand_ROI_pool                       = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        end
        for iter = 1:ml_parameters.N_iter   
            fprintf(['Iteration : ',num2str(iter),'\n'])
            try
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(rand_ROI_pool(randperm(numel(rand_ROI_pool),numel(predictor_ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) rand_ROI_pool(randperm(numel(rand_ROI_pool),x)), cellfun(@numel, predictor_ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            catch
                disp('WARNING : Not enough non-used ROIs available. Using all ROIs instead for randomization')
                ok = find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list));
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(ok(randperm(numel(ok),numel(predictor_ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) ok(randperm(numel(ok),x)), cellfun(@numel, predictor_ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            end
            for gp_idx = 1:numel(ROI_groups_rdm)
                data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups_rdm{gp_idx}),2)';        
            end   
            results{iter}               = train_and_test(data, processed_behaviours, timepoints, 1:numel(ROI_groups_rdm), observations, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            results{iter}.used_ROIs     = ROI_groups_rdm;
        end
    end

    [meanvalue,~, fig_handle, stats, values] = bar_chart(results, '','','','',ml_parameters.rendering, true, ml_parameters);
    if ml_parameters.rendering
        title(ml_parameters.title)
    end
    if iscell(ml_parameters.title)
        ml_parameters.title = [ml_parameters.title{:}];
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
    if ml_parameters.save
        if contains(ml_parameters.save, '.mat')
            ml_parameters.save = parse_paths(fileparts(ml_parameters.save));
        end
        indexes = predictor_ROI_groups;
        save([ml_parameters.save, ml_parameters.title], 'results', 'indexes', 'stats', 'values', '-v7.3')
    end
end