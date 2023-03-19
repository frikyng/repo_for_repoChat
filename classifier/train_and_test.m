%% Given some predictors and observation classes, train a classfier with a 
%% part of the points, and predict the remaining data

%% Note : hold out fraction is controlled by the HOLDOUT() global setting

function out = train_and_test(predictor_data, observation_data, timepoints, roi_subset, beh_types, raw_behaviour, calcium_ref, ml_parameters)
    if nargin < 3 || isempty(timepoints)
        timepoints      = 1:size(observation_data, 2);
        raw_behaviour   = [];
        calcium_ref     = [];
    end
    if nargin < 4 || isempty(roi_subset) %% Define the ROIs we want to use for the training   
        roi_subset    = 1:size(predictor_data, 1);   % cl_4; %45:5:70 in 2019-09-26_exp_1
    end
    if nargin < 5 || isempty(beh_types) %% Define the  names of the behaviour, to display on the plots   
        beh_types    = strcat('beh', strsplit(num2str(1:size(observation_data, 1)),' '));  
    end
    if nargin < 6 || isempty(raw_behaviour) %% Original behaviour trace, for context  
        raw_behaviour    = NaN(numel(beh_types), 1);  
    end
    if nargin < 7 || isempty(calcium_ref) %% Original calcium trace, for context  
        calcium_ref    = [];  
    end
    if nargin < 8 || isempty(ml_parameters)
        ml_parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        ml_parameters = DEFAULT_CLASSIFIER_OPTION(ml_parameters);
    end
    
    USABLE_BEH                      = ~all(isnan(observation_data),2);
    INVALID_TP                      = any(isnan(observation_data(USABLE_BEH,:)),1);
    observation_data(:,INVALID_TP)  = [];
    predictor_data(:,INVALID_TP)    = [];
    timepoints(INVALID_TP)          = [];
    
    if any(find(all(isnan(predictor_data'))))
        warning('Some predictors contains only NaN. This should be filtered out before caling this function');
    end
    
    if ischar(ml_parameters.shuffling)
        ml_parameters.shuffling = {ml_parameters.shuffling};
    end
    
    %% Time Shuffle Observations if required. 
    shuffle_observation_data = [];
    if any(contains(ml_parameters.shuffling, 'behaviours','IgnoreCase',true))
        if ml_parameters.obs_shuf_block_sz == 1
            shuffle_observation_data = observation_data(:,randperm(size(observation_data,2)));
        else
            beh_part       = block_shuffle(timepoints, round(ml_parameters.obs_shuf_block_sz), size(raw_behaviour, 2), 1);
            shuffle_observation_data = observation_data(:,beh_part.x_test);
        end
    end
      
    %% Time or Space Shuffle predictors if required. 
    if any(contains(ml_parameters.shuffling, 'events','IgnoreCase',true))
        predictor_data = predictor_data(:,randperm(size(predictor_data,2)));
    end
    if any(contains(ml_parameters.shuffling, 'ROIs','IgnoreCase',true))
        for tp = 1:size(predictor_data,2)
            predictor_data(:,tp) = predictor_data(randperm(size(predictor_data,1)), tp);
        end
    end
    
    %% Get a random set of timepoints for training vs testing
    if ml_parameters.block_shuffling
        partition       = block_shuffle(timepoints, round(ml_parameters.block_shuffling), size(raw_behaviour, 2), ml_parameters.holdout);
    elseif ~ml_parameters.holdout
        partition       = [];
    else
        partition       = cvpartition(numel(timepoints), 'HoldOut', ml_parameters.holdout);
    end
    
    if any(contains(ml_parameters.shuffling, 'behaviours','IgnoreCase',true))
        shuffled_beh  = strcat(beh_types, '_shuffled');
        beh_types = reshape([beh_types; shuffled_beh],[],1)';       
    end

    score       = [];
    out         = {};
    for mdl_idx = 1:numel(beh_types)
        type                = strrep(beh_types{mdl_idx},' ','_');
        type_corrected      = strrep(type,'\_','_');
        beh_idx = mdl_idx;
        if any(contains(ml_parameters.shuffling, 'behaviours','IgnoreCase',true))
            beh_idx = ceil(mdl_idx/2);
        end
        if contains(type, 'shuffle')
            current_obs         = shuffle_observation_data(beh_idx,:);
        else
            current_obs         = observation_data(beh_idx,:);
        end
        current_obs         = normalize(current_obs);
        %         assignin('base',    ['var',type_corrected], observation_data(el,:))    % why 3 rows of the same?
        %         assignin('base',    'Predictors',           predictor_data(roi_subset,:))

        %% Compute cost matrix
        if islogical(current_obs)
            ratio   = sum(current_obs)/numel(current_obs);
            cost    = [0,ratio;1-ratio,0];
            % assignin('base',    ['cost_',type_corrected],   cost)
        else
            cost = [];
        end   
        
        if all(isnan(current_obs))
            warning([type_corrected, ' has only NaNs']);
            out.calcium         = calcium_ref;
            out.bin_beh{mdl_idx}     = current_obs;
            out.peak_tp{mdl_idx}     = timepoints;
            out.train_range{mdl_idx} = [];  
            out.prediction{mdl_idx}  = [];
            out.full_beh{mdl_idx}    = raw_behaviour(mdl_idx,:);
            out.beh_type{mdl_idx}    = type_corrected;
            out.score{mdl_idx}       = NaN(1,4);
            out.model{mdl_idx}       = {};
            continue
        end
        
        %% Predict behaviour
        [y_predict, y_test, cross_val_score, x_test, x_train, y_train, model] = prediction(predictor_data(roi_subset,:), current_obs, partition, cost, merge_params_obj(ml_parameters, struct('behaviour',type_corrected)));

        %% Get accuracy score
        if isempty(y_test) && ml_parameters.kFold > 1
            %  score(el,:) = repmat(kfoldLoss(model)*100,1,4); % reveals the fraction of predictions that were incorrect, i.e. (1 - accuracy)
            ml_parameters.rendering = 1;
            temp = kfoldLoss(model, 'Mode','individual', 'LossFun', @pearson_correlation_coefficient)*100; 
            score(mdl_idx,:) = repmat(nanmean(temp), 1, 4);
        else
            [score(mdl_idx,1), score(mdl_idx,2), score(mdl_idx,3), score(mdl_idx,4)] = get_classifier_score(y_test, y_predict);
        end

        %% Plot training result %% QQ RAW BEHAVIOUR IS NEVER SHUFFLED SO IT WONT LOOK RIGHT
        if ml_parameters.rendering >= 2              
            if ~isempty(raw_behaviour)
                current_raw_behaviour = raw_behaviour(beh_idx,:);
            end
            plot_prediction(calcium_ref, current_obs, timepoints, partition, y_predict, current_raw_behaviour , type, score(mdl_idx,:), mdl_idx);
        end

        %% Store result
        out.calcium         = calcium_ref;
        out.bin_beh{mdl_idx}     = current_obs;
        out.peak_tp{mdl_idx}     = timepoints;
        if ~isempty(y_test) && ~ml_parameters.block_shuffling
            out.train_range{mdl_idx} = training(partition);
        elseif ml_parameters.block_shuffling
            out.train_range{mdl_idx} = partition.x_test;
        else
            out.train_range{mdl_idx} = [];
        end
        out.prediction{mdl_idx}  = y_predict;
        out.full_beh{mdl_idx}    = raw_behaviour(beh_idx,:);
        out.beh_type{mdl_idx}    = type_corrected;
        out.shuffling_type{mdl_idx}= ml_parameters.shuffling;
        out.score{mdl_idx}       = score(mdl_idx,:);
        out.model{mdl_idx}       = model;
    end
    if numel(beh_types) > 1 && ml_parameters.rendering >= 1
        labels = reordercats(categorical(beh_types),beh_types);
        %figure(1000);clf();bar(labels, score(:,4));
    else
       % fprintf(['Prediction score is ',num2str(out.score{el}(end)),'\n'])
    end
    %arrangefigures; drawnow
end

function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
    score = corr(y_true, y_pred); % ,'Type','Spearman'
end