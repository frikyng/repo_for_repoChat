%% Given some predictors and observation classes, train a classfier with a 
%% part of the points, and predict the remaining data

function out = train_and_test(predictor_data, observation_data, timepoints, roi_subset, method, beh_types, raw_behaviour, calcium_ref)
    if nargin < 3 || isempty(timepoints)
        timepoints      = 1:size(observation_data, 2);
        raw_behaviour   = [];
        calcium_ref     = [];
    end
    if nargin < 4 || isempty(roi_subset) %% Define the ROIs we want to use for the training   
        roi_subset    = 1:size(predictor_data, 1);   % cl_4; %45:5:70 in 2019-09-26_exp_1
    end
    if nargin < 5 || isempty(method) %% defines clasfier or learner method 
        method    = 'svm';  
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
    
    %% Get a random set of timepoints for traing vs testing
    rand_tp     = randperm(numel(timepoints));       
    score       = [];
    out         = {};
    for el = 1:numel(beh_types)
        type = strrep(beh_types{el},' ','_');
        type_corrected = strrep(type,'\_','_');
        current_var = observation_data(el,:);
        assignin('base',['var',type_corrected],observation_data(el,:))
        assignin('base','Predictors',predictor_data(roi_subset,:))
        %eval(['var',type_corrected,'=all_beh(',num2str(el),',:);'])

        train_R     = rand_tp(1:floor(numel(timepoints)/2)); 

        if contains(type, 'trigger')
            cost = [0,1;10,0];
        else
            cost = [0,10;1,0];
        end   
        
        %phate(double(Ca(roi_subset,:)'),'ndim',10)'
        if islogical(current_var)
            [y_predict, y_test, cross_val_score, x_test, x_train, y_train] = classifier_prediction(predictor_data(roi_subset,:), current_var, train_R, method, cost);
        else
            [y_predict, y_test, cross_val_score, x_test, x_train, y_train] = regression_prediction(predictor_data(roi_subset,:), current_var, train_R, method, cost);
        end
        [score(el), sensitivity, specifity] = get_classifier_score(y_test, y_predict);
        
        
        %     all_scores = [];
        %     for roi = 1:size(Ca, 1)
        %         roi
        %         temp_Ca = Ca;
        %         temp_Ca(roi,:) = [];
        %         [~, ~, all_scores(roi), ~, ~, ~] = csvm_prediction(temp_Ca, current_var, train_R, cost);        
        %     end
        
        if ~isempty(raw_behaviour)
            current_raw_behaviour = raw_behaviour(el,:);
        end
        plot_prediction(calcium_ref, current_var, timepoints, train_R, y_predict, current_raw_behaviour , type, score(el), el);

        out.calcium         = calcium_ref;
        out.bin_beh{el}     = current_var;
        out.peak_tp{el}     = timepoints;
        out.train_range{el} = train_R;
        out.prediction{el}  = y_predict;
        out.full_beh{el}    = raw_behaviour(el,:);
        out.beh_type{el}    = type_corrected;
        out.score{el}       = score(el);
    end
    figure(el+1);clf();bar(categorical(beh_types), score);
    arrangefigures