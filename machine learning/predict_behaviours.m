%% Machine learning function to see if the neuron signal predcit behaviours
% -------------------------------------------------------------------------
%% Syntax:
%   [results, data, predictor_ROI_groups, mean_score, stats, ind_scores] =
%               predict_behaviours(obj, predictors_type, observations,
%                                   predictor_ROI_groups, varargin)
%
% -------------------------------------------------------------------------
%% Inputs:
%   obj(arboreal_scan_experiment object, or CHAR) - Optional - Default is
%           empty
%       obj can be :
%       * An arboreal_scan_experiment object (Recommended)
%       * A STR PATH to an extracted arboreal_scan_experiment
%       * an empty value, in whic case it will look for the current folder
%
%   predictors_type(CHAR or CHAR ARRAY) - Optional - Default is
%       subtracted_peaks
%       * Use any valid behaviour / signal filter syntax, as described in
%       prepare_analysis.m. This will determine the signal timepoints that
%       will be used (if you filter for specific behaviours, type of
%       activity or introduce any lag) and will determine the type of
%       signal used (raw traces, rescaled, meidan-subtracted etc...)
%
%   observations(CHAR or CHAR ARRAY or STRUCT) - Optional - Default is 
%           obj.behaviours.types; 
%       This corresponds to the behaviours that will be predicted. 
%       * If you provide the name of the behaviour, it can be a single 
%       behaviour (eg. 'encoder') or a list of behaviours (eg {'trigger', 
%       'encoder'}). Each element must correspond to an valid behaviour 
%       that can be obtained from obj.behaviours.types. behaviours are then
%       extracted using prepare_behaviour_data()
%       * You can directly provide a structure with the different outputs
%       of prepare_behaviour_data(), as there is a lot of overhead calling 
%       this function. If you run different instanciations of the machine
%       learning, but if the behaviours are always the same, you can build
%       your own input. See extra Notes and Demo_FINAL_SERIES_OF_ANALYSIS
%
%
%   predictor_ROI_groups([] or 1xN MATRIX or 1xM CELL ARRAY of 1xN
%       MATRICES) - Optional - Default is obj.ref.indices.valid_swc_rois.
%       Defines the predictors indices to use
%       * If you pass [] or a MATRIX, each value listed will be used as an
%       independent predictor.
%       * If you pass a cell array, each cell will be used as an
%       independent predictor. If a cell has more than one value, the
%       corrsponding predictors are averaged
%
%   varargin(machine_learning_params object AND/OR {'Argument',Value} 
%       pairs):  
%       Any pair of 'argument' and value. For details and default values,
%       see machine_learning_params
% -------------------------------------------------------------------------
%% Outputs:
%
%	results (1xN CELL ARRAY of MACHINE LEARNING OUPTUT STRUCTURE) :
%       * Structure containing one cell per model training (as defined by
%       mlp.N_Iter. Each cell contains a structure with the following
%       fields. Fields indicated in parenthesis are empty if
%       mlp.saving_level
%           - out.mlp % a copy of the machine learning settings
%           - out.(calcium) % full input predictors (all tp)
%           - out.(full_beh) % full input observations (all tp)
%           - out.(predictors) % input predictors (selected tp)
%           - out.observation % input observation (selected tp) 
%           - out.timepoints % selected tp
%           - out.train_range_idx % points used for training (empty if
%                                   only kfold)
%           - out.test_range_idx %  points used for testing
%           - out.prediction %  predicted observation
%           - out.beh_type %  behaviour name
%           - out.shuffling_type %  type of data shuffling if any
%           - out.score %  all computed scores
%           - out.model %  the ml model
%
%   data(N x T NUMERIC)
%       The predictor data used for the training
%
%   predictor_ROI_groups(1xN INT)
%       The sets of ROIs used for the training (after some correction
%       compared to the input) Note that in the case of randomization, you
%       need to look at obj.used_ROI for the correct model iteration
%
%    mean_score(1 x N)
%       The average score for all iteration, for each N condition, as 
%       provided by the bar_chart output
%
%    stats(STRUCT)
%       a structure with the statistic results across conditions
%
%    ind_scores(M x N)
%       The score for all M iteration, for each N condition, as 
%       provided by the bar_chart output
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * Generating observations can be long if you call this function many time.
%   To prevent this, you can call 
%       extracted_beh   = prepare_behaviour_data(obj, behaviours);
%   This will generate a structure can be passed as an input. Timepoints
%   filters and lag are still aplied on top of it. This method could also 
%   be used to pass custom behaviours that are not in obj.behaviours.types
%
% -------------------------------------------------------------------------
%% Examples:
%
% * See Demo_FINAL_SERIES_OF_ANALYSIS
%
% * obj is a preprocessed arboreal_scan_experiment object. You need to have
%   extracted peaks and rescaled traces, which is sensitive to factors such
%   as detrending and smoothing. If you change detrending an smoothing, you
%   will have to re-extract the peaks and rescale again
%
% * observations is typically a single behaviour or a list of behaviour
%   names. Timepoints are the same than the one set in predictors_type,
%   although it can have a temporal lag between the two.
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%    
% -------------------------------------------------------------------------
%% Revision Date:
%   17-04-2023
% -------------------------------------------------------------------------
%% See also:
%    prepare_analysis, train_and_test, bar_chart


% use expe 07-26 to fix gaps in behaviour
% Cell 2019-09-17_exp_1 anticorelated to running

function [results, data, predictor_ROI_groups, mean_score, stats, ind_scores] = predict_behaviours(obj, predictors_type, observations_labels, predictor_ROI_groups, varargin)
    if nargin < 1 || isempty(obj) || ischar(obj)
        %% Get object if it needs loading/building
        if ischar(obj)
        	obj = arboreal_scan_experiment(obj, false);
        else
            obj = arboreal_scan_experiment(pwd, false);
        end
    elseif ~isa(obj, 'arboreal_scan_experiment')
    	disp_info('Input should be an erboreal_scan_experiment or a path pointing to an arboreal_scan_experiment', 4);
    end
    if nargin < 2 || isempty(predictors_type)
        predictors_type   = 'subtracted_peaks'; % ['subtracted' OR 'rescaled' OR 'raw'] AND ['peaks' or '']. eg 'subtracted_peaks' , or 'raw'
    end
    build_beh           = true;
    if nargin < 3 || isempty(observations_labels)
        % pass      
    elseif isstruct(observations_labels)
        raw_behaviours              = observations_labels.raw_behaviours;
        beh_thr                     = observations_labels.beh_thr;
        formatted_behaviour_list    = observations_labels.formatted_behaviour_list;
        observations_labels                = observations_labels.original_behaviour_list;
        build_beh                   = false;
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

    %% For now, no HD DATA
    use_hd_data             = false;
    debug                   = false;
    
    %% Make sure preprocessing was done correctly and adjust the predictors
    rendering                                   = obj.rendering;
    obj.rendering                               = false;
    [obj, source_signal, ~, timepoints, lag]    = prepare_analysis(obj, use_hd_data, [], predictors_type);
    
    %% Introduce time lag if required
    beh_timepoints                              = timepoints + lag; % if lag == 0 then points are the same
    timepoints(beh_timepoints < 1)              = [];
    beh_timepoints(beh_timepoints < 1)          = [];
    timepoints(beh_timepoints > obj.tp)         = [];
    beh_timepoints(beh_timepoints > obj.tp)     = [];
    
    if debug
        holder = source_signal(:,1); figure;plot(source_signal(:,1));hold on; scatter(timepoints, holder(timepoints),'r','filled')
        figure();hist(reshape(source_signal(timepoints,:),[],1),100); % TOMMY uncomment to see distribution of predictors
    end
    
    %% Restore rendering setting
    obj.rendering       = rendering;
    
    %% If required add shufffing
    if ml_parameters.add_shuffle
        shuffled_beh            = strcat(observations_labels, '_shuffled');
        observations_labels     = reshape([observations_labels; shuffled_beh],[],1)';
    end
    
    %% Get observations. If not provided, compute the values from behaviours
    if ~build_beh
        %% If we already had the extracted behaviours, at this stage we just update the formatted name 
    	observations_labels    = formatted_behaviour_list;
    else
        [~, raw_behaviours, beh_thr, observations_labels] = prepare_behaviour_data(obj, observations_labels);
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
    
    if debug
    	figure();imagesc(tril(corr([data',processed_behaviours'])',-1));
    end

    %% Remove invalid ROi groups
    fully_invalid_group                         = cellfun(@isempty, predictor_ROI_groups)' | all(isnan(data),2);
    predictor_ROI_groups(fully_invalid_group)   = [];
    data(fully_invalid_group, :)                = [];

    if debug
        % Correlation of the different ROIs with each other
        figure();imagesc(corr(data'));    
        % Correlation of the Ca2+ values per event
        figure();imagesc(corr(data));
        % Correlation of the different behaviours with each other
        figure();imagesc(corr(processed_behaviours'));   
        % Correlation of the behaviours values per event
        figure();imagesc(corr(processed_behaviours));
    end

    %% Control block if you have a single ROI
    %     data = repmat(data, 50, 1);
    %     data = data + randn(size(data))*50;
    
    %% FYI find(all(isnan(source_signal))) should be empty, or you have observation with only NaN    
    All_ROIs    = 1:size(data, 1);

    %% Now run training on regular OR randomized data
    if ~ml_parameters.randomize_ROIs
        %% Regular case, run N_iter
        for iter = 1:ml_parameters.N_iter
            fprintf(['Iteration : ',num2str(iter),'\n'])
            results{iter}               = train_and_test(data, processed_behaviours, timepoints, All_ROIs, observations_labels, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            results{iter}.used_ROIs     = predictor_ROI_groups;
        end
    else  
        %% We build alternative randomized groups. 
        % 0/false does no randomization
        % 1 randomize groups using all valid ROIs
        % -1 randomize groups using all valid ROIs, excluding the ROIs listed in the groups. Groups are sized matched
        if ml_parameters.randomize_ROIs == 1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        elseif ml_parameters.randomize_ROIs == -1
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois((~ismember(obj.ref.indices.valid_swc_rois, horzcat(predictor_ROI_groups{:}))) & (~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list))); %list of valid ROIs excluding ROI groups
        end
        
        %% When you pass all ROIs (so all cells in ROI_groups size == 1), randomization for leftover ROIs makes no sense so we just randomize ROIs
        if isempty(rand_ROI_pool)
            disp('WARNING : No ROI left for randomization')
            ml_parameters.randomize_ROIs        = 1;
            rand_ROI_pool                       = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); %list of all valid ROIs 
        end
        
        %% Now run model for each iteration, with different randomization every time using the pool of ROis that are authorized for randomization
        for iter = 1:ml_parameters.N_iter   
            fprintf(['Iteration : ',num2str(iter),'\n'])
            try
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(rand_ROI_pool(randperm(numel(rand_ROI_pool),numel(predictor_ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) rand_ROI_pool(randperm(numel(rand_ROI_pool),x)), cellfun(@numel, predictor_ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            catch
                disp_info('WARNING : Not enough non-used ROIs available. Using all ROIs instead for randomization', 2)
                ok      = find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list));
                if single_matrix_input % when input was matrix
                    ROI_groups_rdm     = num2cell(ok(randperm(numel(ok),numel(predictor_ROI_groups))));
                else
                    ROI_groups_rdm     = cellfun(@(x) ok(randperm(numel(ok),x)), cellfun(@numel, predictor_ROI_groups, 'UniformOutput', false), 'UniformOutput', false); % size matched groups from ROI_pool
                end
            end
            
            %% Update group info
            for gp_idx = 1:numel(ROI_groups_rdm)
                data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups_rdm{gp_idx}),2)';        
            end   
            
            %% Train and test model
            results{iter}               = train_and_test(data, processed_behaviours, timepoints, 1:numel(ROI_groups_rdm), observations_labels, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
            results{iter}.used_ROIs     = ROI_groups_rdm;
        end
    end

    %% Display results and/or compute stats
    [mean_score,~, fig_handle, stats, ind_scores] = bar_chart(results, '','','','',ml_parameters.rendering, true, ml_parameters);
    if ml_parameters.rendering
        title(ml_parameters.title)
    end
    if iscell(ml_parameters.title)
        ml_parameters.title = [ml_parameters.title{:}];
    end
    
    %% Save figure if required
    if ml_parameters.savefig
        if islogical(ml_parameters.savefig)
            save_myfig(fig_handle,ml_parameters.title,{'png','pdf'})
        elseif ischar(ml_parameters.savefig)
            if contains(ml_parameters.savefig, '.mat')
                ml_parameters.savefig = parse_paths(fileparts(ml_parameters.savefig));
            end
            save_myfig(fig_handle,[ml_parameters.savefig, '\', ml_parameters.title],{'png','pdf'})
        end
    end 
    
    %% Save result if required
    if ml_parameters.save
        if contains(ml_parameters.save, '.mat')
            ml_parameters.save = parse_paths(fileparts(ml_parameters.save));
        end
        indexes = predictor_ROI_groups;
        save([ml_parameters.save,'\',  ml_parameters.title], 'results', 'indexes', 'stats', 'ind_scores', '-v7.3')
    end
end