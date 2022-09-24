function [y_predict, y_test, score, x_test, x_train, y_train] = csvm_prediction(XData, YData, train_tp, cost)
    %% Split data between train and test
    train_tp    = sort(train_tp);
    all_tp      = 1:size(XData, 2);
    test_tp     = all_tp(~ismember(all_tp, train_tp));

    x_train     = double(XData(:,train_tp))';   % Predictors for training
    x_test      = double(XData(:,test_tp))';    % Predictors for testing
    y_train     = YData(:, train_tp);              % Variable to predict used for training  
    y_test      = YData(:, test_tp);               % Variable to predict used for testing (ground truth)     

    type = 'svm'
    if strcmp(type, 'forest')
        classificationTree = TreeBagger(500, x_train, y_train, 'OOBPrediction', 'on','Method', 'classification','Cost', cost); 
        y_predict = str2double(classificationTree.predict(x_test));
    elseif strcmp(type, 'svm')
        %% Generate classifier % fitcsvm
        classificationSVM = fitcsvm(...
                                    x_train, ...
                                    y_train, ...
                                    'KernelFunction', 'gaussian', ...
                                    'PolynomialOrder', [], ...
                                    'KernelScale', 13, ...
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


        %% Buld predicted variable using all different cross valdiated models
        y_predict = [];
        for trained_idx = 1:numel(classificationSVM.Trained)
            mdl = classificationSVM.Trained{trained_idx};                        
            y_predict(trained_idx, :) = predict(mdl, x_test);
        end
        y_predict = logical(round(nanmean(y_predict)))';

        %% Estimate accuracy of the corr validated model, and of the predicted values vs ground thruth
        score = kfoldLoss(classificationSVM)*100;
        fprintf(['The out-of-sample misclassification rate is ',num2str(score,3),'%%\n'])
      
    end
    score = 100 * sum((y_predict - y_test') == 0) / numel(y_predict);
end