        %% Need to discuss 
        %% 'Standardize' — Flag to standardize predictor data --> false (default) | true
        %% 'Solver' — Optimization routine --> 'ISDA' | 'L1QP' | 'SMO'
        %% 'Nu' — ν parameter for one-class learning --> 0.5 (default) | positive scalar
        %% 'Weights' — Observation weights --> numeric vector | name of variable in Tbl % should we use SNR?

function [y_predict, y_test, score, x_test, x_train, y_train] = classifier_prediction(XData, YData, train_tp, method, cost)
    %% Split data between train and test
    train_tp    = sort(train_tp);
    all_tp      = 1:size(XData, 2);
    test_tp     = all_tp(~ismember(all_tp, train_tp));

    x_train     = double(XData(:,train_tp))';   % Predictors for training
    x_test      = double(XData(:,test_tp))';    % Predictors for testing
    y_train     = YData(:, train_tp);              % Variable to predict used for training  
    y_test      = YData(:, test_tp);               % Variable to predict used for testing (ground truth)     

    if strcmp(method, 'forest')
        classificationTree = TreeBagger(500, x_train, y_train, 'OOBPrediction', 'on','Method', 'classification','Cost', cost); 
        y_predict = str2double(classificationTree.predict(x_test));
    elseif strcmp(method, 'svm')        
        %% Generate classifier % fitcsvm
        classificationSVM = fitcsvm(...
                                    x_train, ...
                                    y_train, ...
                                    'KernelFunction', 'gaussian', ...
                                    'Nu', 0.5,...
                                    'KernelScale', 'auto', ...
                                    'BoxConstraint', 1, ...
                                    'Standardize', true, ...
                                    'ClassNames', [false; true],...
                                    'CrossVal','on',...
                                    'KFold',10,...
                                    'Cost', cost); % [0,10;1,0]
                                
%         classificationSVM = fitclinear(...
%                                     x_train, ...
%                                     y_train,...
%                                     'KFold',10,...
%                                     'Regularization','lasso'); % [0,10;1,0]

        y_predict = get_consensus_prediction(classificationSVM, x_test);

        %% Estimate accuracy of the corr validated model, and of the predicted values vs ground thruth
        score = kfoldLoss(classificationSVM)*100;
        fprintf(['The out-of-sample misclassification rate is ',num2str(score,3),'%%\n']);  
    else
        error('classification method not implemented')
    end
end

function y_predict = get_consensus_prediction(Classifier, x_test)
    %% Buld predicted variable using all different cross valdiated models
    y_predict = [];
    for trained_idx = 1:numel(Classifier.Trained)
        mdl = Classifier.Trained{trained_idx};                        
        y_predict(trained_idx, :) = predict(mdl, x_test);
    end
    y_predict = logical(round(nanmean(y_predict)))';
end