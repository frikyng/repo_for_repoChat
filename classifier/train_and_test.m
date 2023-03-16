%% Given some predictors and observation classes, train a classfier with a 
%% part of the points, and predict the remaining data

%% Note : hold out fraction is controlled by the HOLDOUT() global setting

function out = train_and_test(predictor_data, observation_data, timepoints, roi_subset, method, beh_types, raw_behaviour, calcium_ref, parameters)
    if nargin < 3 || isempty(timepoints)
        timepoints      = 1:size(observation_data, 2);
        raw_behaviour   = [];
        calcium_ref     = [];
    end
    if nargin < 4 || isempty(roi_subset) %% Define the ROIs we want to use for the training   
        roi_subset    = 1:size(predictor_data, 1);   % cl_4; %45:5:70 in 2019-09-26_exp_1
    end
    if nargin < 5 || isempty(method) %% defines clasfier or learner method 
        method    = 'svm';  % i.e. fitcsvm vs fitclinear etc
    end    
    if nargin < 6 || isempty(beh_types) %% Define the  names of the behaviour, to display on the plots   
        beh_types    = strcat('beh', strsplit(num2str(1:size(observation_data, 1)),' '));  
    end
    if nargin < 7 || isempty(raw_behaviour) %% Original behaviour trace, for context  
        raw_behaviour    = NaN(numel(beh_types), 1);  
    end
    if nargin < 8 || isempty(calcium_ref) %% Original calcium trace, for context  
        calcium_ref    = [];  
    end
    if nargin < 9 || isempty(parameters)
        parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        parameters = DEFAULT_CLASSIFIER_OPTION(parameters);
    end
    
    USABLE_BEH                      = ~all(isnan(observation_data),2);
    INVALID_TP                      = any(isnan(observation_data(USABLE_BEH,:)),1);
    observation_data(:,INVALID_TP)  = [];
    predictor_data(:,INVALID_TP)    = [];
    timepoints(INVALID_TP)          = [];
    
    if any(find(all(isnan(predictor_data'))))
        warning('Some predictors contains only NaN. This should be filtered out before caling this function');
    end
    
    %% Shuffle predictor if required.     
    if strcmpi(parameters.shuffling, 'events') || strcmpi(parameters.shuffling, 'both')
        predictor_data = predictor_data(:,randperm(size(predictor_data,2)));
    end
    if strcmpi(parameters.shuffling, 'ROIs') || strcmpi(parameters.shuffling, 'both')
        for tp = 1:size(predictor_data,2)
            predictor_data(:,tp) = predictor_data(randperm(size(predictor_data,1)), tp);
        end
    end
    
    %% Get a random set of timepoints for training vs testing
    if parameters.block_shuffling
        partition       = block_shuffle(timepoints, round(parameters.block_shuffling), size(raw_behaviour, 2), parameters.holdout);
    elseif ~parameters.holdout
        partition       = [];
    else
        partition       = cvpartition(numel(timepoints), 'HoldOut', parameters.holdout);
    end

    score       = [];
    out         = {};
    for el = 1:numel(beh_types)
        type                = strrep(beh_types{el},' ','_');
        type_corrected      = strrep(type,'\_','_');
        current_var         = observation_data(el,:);
        current_var         = normalize(current_var);
        assignin('base',    ['var',type_corrected], observation_data(el,:))    % why 3 rows of the same?
        
        
        assignin('base',    'Predictors',           predictor_data(roi_subset,:))

        %% Compute cost matrix
        if islogical(current_var)
            ratio   = sum(current_var)/numel(current_var);
            cost    = [0,ratio;1-ratio,0];
            assignin('base',    ['cost_',type_corrected],   cost)
        else
            cost = [];
        end   
        
        if all(isnan(current_var))
            warning([type_corrected, ' has only NaNs']);
            out.calcium         = calcium_ref;
            out.bin_beh{el}     = current_var;
            out.peak_tp{el}     = timepoints;
            out.train_range{el} = [];  
            out.prediction{el}  = [];
            out.full_beh{el}    = raw_behaviour(el,:);
            out.beh_type{el}    = type_corrected;
            out.score{el}       = NaN(1,4);
            out.model{el}       = {};
            continue
        end
        
        %% Predict behaviour
        [y_predict, y_test, cross_val_score, x_test, x_train, y_train, model] = prediction(predictor_data(roi_subset,:), current_var, partition, method, cost, merge_params_obj(parameters, struct('behaviour',type_corrected)));

        %% Get accuracy score
        if isempty(y_test) && parameters.kFold > 1
            %  score(el,:) = repmat(kfoldLoss(model)*100,1,4); % reveals the fraction of predictions that were incorrect, i.e. (1 - accuracy)
            parameters.rendering = 1;
            temp = kfoldLoss(model, 'Mode','individual', 'LossFun', @pearson_correlation_coefficient)*100; 
            score(el,:) = repmat(nanmean(temp), 1, 4);
        else
            [score(el,1), score(el,2), score(el,3), score(el,4)] = get_classifier_score(y_test, y_predict);
        end
        
        %% Plot training result
        if ~isempty(raw_behaviour)
            current_raw_behaviour = raw_behaviour(el,:);
        end
        
        if parameters.rendering >= 2    
            plot_prediction(calcium_ref, current_var, timepoints, partition, y_predict, current_raw_behaviour , type, score(el,:), el);
        end

        %% Store result
        out.calcium         = calcium_ref;
        out.bin_beh{el}     = current_var;
        out.peak_tp{el}     = timepoints;
        if ~isempty(y_test) && ~parameters.block_shuffling
            out.train_range{el} = training(partition);
        elseif parameters.block_shuffling
            out.train_range{el} = partition.x_test;
        else
            out.train_range{el} = [];
        end
        out.prediction{el}  = y_predict;
        out.full_beh{el}    = raw_behaviour(el,:);
        out.beh_type{el}    = type_corrected;
        out.score{el}       = score(el,:);
        out.model{el}       = model;
    end
    if numel(beh_types) > 1 && parameters.rendering >= 1
        labels = reordercats(categorical(beh_types),beh_types);
        %figure(1000);clf();bar(labels, score(:,4));
    else
       % fprintf(['Prediction score is ',num2str(out.score{el}(end)),'\n'])
    end
    arrangefigures; drawnow
end

function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
    score = corr(y_true, y_pred); % ,'Type','Spearman'
end


function partition = block_shuffle(timepoints, shuffle_window_pts, total_tp, holdout)

    shuffle_window_pts = abs(shuffle_window_pts);

	%% Cut data into windows
    win_start = 1:shuffle_window_pts:total_tp;
    windows  = {};
    for t = win_start
         current_win = t:t+shuffle_window_pts-1;
         windows{end+1} = timepoints(ismember(timepoints, current_win));      
    end
    windows = windows(~cellfun(@(x) isempty(x), windows));
    
    %% Shuffle windows
    windows = windows(randperm(numel(windows)));
    
    %% Extract held out blocks
    n_pts_target = ceil(numel(timepoints) * holdout);
    n_cells_to_use = find(cumsum(cellfun(@(x) numel(x), windows)) < n_pts_target, 1, 'last');
    heldout_batch = windows(1:n_cells_to_use);
    training_batch = windows(n_cells_to_use+1:end);

    partition = {};
    partition.x_train = [training_batch{:}];
    partition.x_test  = [heldout_batch{:}];
    
    [~, partition.x_train] = ismember(partition.x_train, timepoints);
    [~, partition.x_test] = ismember(partition.x_test, timepoints);
%     %% Rearrange order
%     timepoints = [windows{:}]; 
%     data = data(:, timepoints);
%     processed_behaviours = processed_behaviours(:,timepoints);
%     

end