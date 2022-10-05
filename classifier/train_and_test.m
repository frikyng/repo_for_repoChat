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

    %% Get a random set of timepoints for traing vs testing
    partition     = cvpartition(numel(timepoints), 'HoldOut', parameters.holdout); % randperm(numel(timepoints));       
    score       = [];
    out         = {};
    for el = 1:numel(beh_types)
        type                = strrep(beh_types{el},' ','_');
        type_corrected      = strrep(type,'\_','_');
        current_var         = observation_data(el,:);
        assignin('base',    ['var',type_corrected], observation_data(el,:))
        assignin('base',    ['var',type_corrected], observation_data(el,:))
        assignin('base',    ['var',type_corrected], observation_data(el,:))
        assignin('base',    'Predictors',           predictor_data(roi_subset,:))

        %% Compute cost matrix
        if islogical(current_var)
            ratio   = sum(current_var)/numel(current_var);
            cost    = [0,ratio;1-ratio,0];
            assignin('base',    ['cost_',type_corrected],   cost)
        else
            cost = [];
        end   
        
        %% Predict behaviour
        [y_predict, y_test, cross_val_score, x_test, x_train, y_train] = prediction(predictor_data(roi_subset,:), current_var, partition, method, cost, parameters);

        %% Get accuracy score
        [score(el,1), score(el,2), score(el,3), score(el,4)] = get_classifier_score(y_test, y_predict);
        
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
        out.train_range{el} = training(partition);
        out.prediction{el}  = y_predict;
        out.full_beh{el}    = raw_behaviour(el,:);
        out.beh_type{el}    = type_corrected;
        out.score{el}       = score(el,:);
    end
    if numel(beh_types) > 1 && parameters.rendering >= 1
        labels = reordercats(categorical(beh_types),beh_types);
        figure(1000);clf();bar(labels, score(:,4));
    else
       % fprintf(['Prediction score is ',num2str(out.score{el}(end)),'\n'])
    end
    arrangefigures; drawnow