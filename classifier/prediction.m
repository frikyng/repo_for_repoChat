function [y_predict, y_test, score, x_test, x_train, y_train, model] = prediction(XData, YData, partition, method,cost, parameters)
    if nargin < 6 || isempty(parameters)
        parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        parameters = DEFAULT_CLASSIFIER_OPTION(parameters);
    end

    %% Define timepoints for train and test
    if isempty(partition) % if we just do many kfolds
        train_tp    = 1:numel(YData);
        test_tp     = [];
    else % if we hold out data for final testing
        train_tp    = find(training(partition));
        test_tp     = find(test(partition));
    end

    %% Assign predictor and observation data to test and train dataset
    x_train     = double(XData(:,train_tp))';       % Predictors for training
    x_test      = double(XData(:,test_tp))';        % Predictors for testing
    y_train     = YData(:, train_tp);               % Variable to predict used for training  
    y_test      = YData(:, test_tp);                % Variable to predict used for testing (ground truth)     

    %% Run the correct classifier or regression
    if strcmpi(method, 'svm') 
        %% https://fr.mathworks.com/help/stats/support-vector-machine-classification.html
        base_varargin = {   x_train         , ...
                            y_train         , ...
                            'KernelFunction', parameters.svm_kernel, ...
                            'Standardize'   , true, ...
                            'CrossVal'      , 'on',...
                            'KFold'         , parameters.kFold};
        
        %% Set training model                
        if islogical(y_train)
            func = @fitcsvm;
            base_varargin = [   base_varargin,...
                                {'Cost', cost,...
                                'ClassNames', [false; true]}];
        else
            func = @fitrsvm;
        end

        %% Define kernel and box size, and update model settings
        [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters);   
        base_varargin = [base_varargin, 'KernelScale', kmax, 'BoxConstraint', bmax];

        %% Train model
        model = func(base_varargin{:});

        %% Test model
        y_predict = get_consensus_prediction(model, x_test); 

        %% Estimate accuracy of the cross-validated model, and of the predicted values vs ground thruth
        score = kfoldLoss(model)*100; % reveals the fraction of predictions that were incorrect, i.e. (1 - accuracy)
        %fprintf(['The model out-of-sample misclassification rate is ',num2str(score,3),'%%\n']);  
    elseif strcmpi(method, 'linear')        
        if isempty(parameters.alpha) % this should be no regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , true,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares'}; % default is no regularization, no standardization (e.g. data is not mean centered)
        elseif parameters.alpha == 0 % ridge regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , true,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares',...
                            'Regularization', 'ridge'}; 
        elseif parameters.alpha == 1 % lasso regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , true,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares',...
                            'Regularization', 'lasso'};
        else
            error('regularization method not implemented')
        end
                              
        %% If you manualy set the solver, adjust it here
        if ~isempty(parameters.solver)                
            base_varargin = [base_varargin, 'Solver', parameters.solver];   
        end
        
        %% IF eleastic Net, adjust alpha
        %         if ~isempty(parameters.alpha) && parameters.alpha ~= 1 && parameters.alpha ~= 0 % ~isnan(parameters.alpha) %
        %             base_varargin = [base_varargin, {'Alpha', parameters.alpha}];
        %         end

        %% Set training model
        if islogical(y_train)
            func = @fitclinear;
            base_varargin = [base_varargin, {'Cost', cost, 'ClassNames', [false; true]}];
        else
            func = @fitrlinear;
        end

        Lambda        = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters);   
        base_varargin = [base_varargin, {'Lambda', Lambda}];
        base_varargin = [base_varargin, {'KFold', 5}]; 
           
        %% Train models
        model         = func(base_varargin{:});

        y_predict     = get_consensus_prediction(model, x_test, x_train); %x_train only when not using held out data 
        if islogical(y_train)
            y_predict = logical(round(y_predict));
        end
    
        score = [];
        
    %% still linear regression but we need to use lasso function and fitrlinear doesn't have alpha as parameter
    elseif strcmpi(method, 'elastic')
       error('to review')
       %convert roi_subset into array of strings so lasso is happy
       %for alpha = [0.1,0.9]
       roi_subset = [1:111];
       roi_subset_cell = num2cell(roi_subset);
       list_ROIs_lasso =cellfun(@num2str,roi_subset_cell,'un',0);
       
       % run lasso, each value in the columns of the first output corresponds to a particular regularization
       % coefficient in Lambda, using a geometric sequence of Lambdas 
       [rglrzd_Lambda_coefs, FitInfo] = lasso( x_train         , ...
                                          y_train         , ... 
                                          'CV', 10        , ...
                                          'PredictorNames', list_ROIs_lasso,...
                                          'Alpha', 0.5);   %...       % params.alpha or use 'Alpha', alpha, where alpha can be between 0 and 1
                                      

        % Lasso plot with x-validated fits
        % look at this plot across a range of alphas
        %lassoPlot(rglrzd_Lambda_coefs,FitInfo,'PlotType','CV'); legend('show');
        
        % Display variable names (or ROI list) in a model that corresponds to minimum x-validated mean squared error
        idxLambdaMinMSE = FitInfo.IndexMinMSE;
        minMSEModelPredictors = FitInfo.PredictorNames(rglrzd_Lambda_coefs(:,idxLambdaMinMSE)~=0);
        
        % get betas that correspond to optimal model (from MinMSE of whatever lambda)
        coef_minMSE = rglrzd_Lambda_coefs(:,idxLambdaMinMSE);
        coef0_minMSE = FitInfo.Intercept(idxLambdaMinMSE);
        
        y_predict = x_test*coef_minMSE + coef0_minMSE;
        score = NaN;
        model = NaN;
        
%         hold on
%         scatter(y_test,y_predict,'red');scatter(x_test,y_test,'blue')
%         plot(y_test,y_test)
%         plot((x_test*coef_minMSE + coef0_minMSE), y_predict)
%         xlabel('Test')
%         ylabel('Predict')
%         hold off
                    
        % Display variable names (or ROI list) in sparsest model within one SE of minimum MSE
        % (this 1SE from the minMSE is a heuristic that is widely accepted by the lasso community and accords with parsimony)
        % but this heuristic doesn't work in this cell.. just gives empty array
%         idxLambda1SE = FitInfo.Index1SE;
%         sparseModelPredictors = FitInfo.PredictorNames(rglrzd_Lambda_coefs(:,idxLambda1SE)~=0);
        
         % get betas that correspond to optimal model (from 1MSE of whatever lambda)
%         coef_1SE = rglrzd_Lambda_coefs(:,idxLambda1SE);
%         coef0_1SE = FitInfo.Intercept(idxLambda1SE);
        
%         y_predict = x_test*coef_1SE + coef0_1SE;
%         hold on
%         scatter(y_test,y_predict)
%         plot(y_test,y_test)
%         xlabel('Test')
%         ylabel('Predict')
%         hold off
        
     %%    
    elseif strcmpi(method, 'forest')
        classificationTree = TreeBagger(500, x_train, y_train, 'OOBPrediction', 'on','Method', 'classification','Cost', cost); 
        y_predict = str2double(classificationTree.predict(x_test));  
    elseif strcmpi(method, 'glm')  
            model =  fitglm(x_train,y_train,'linear','Distribution','inverse gaussian'); %  'binomial' | 'poisson' | 'gamma' | 'inverse gaussian'
            [y_predict,~] = predict(model,x_test);
            score = NaN;
            % perform cross-validation, and return average MSE across folds
            %mse = crossval('mse', x_train,y_train', 'Predfun',fcn, 'kfold',10);
            % compute root mean squared error
            %avrg_rmse = sqrt(mse)

%             mask = ones(size(mdl.Coefficients.SE)-1);
%             mask = abs(mdl.Coefficients.SE(2:end)) < 0.001;
%             mdl2 = removeTerms(mdl, mask')
%             score  = NaN
%             [y_predict2,~] = predict(mdl2,x_test)
%             assignin('base','p',mdl.Coefficients.pValue)
%             assignin('base','SE',mdl.Coefficients.SE)
%             assignin('base','tStat',mdl.Coefficients.tStat)
%             assignin('base','Estimate',mdl.Coefficients.Estimate)
            %obj.ref.plot_value_tree(SE, find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)))
            
%             [y_predict,~] = predict(mdl,x_test)
%             score  = NaN 
%             
%             

    else
        error('classification method not implemented')
    end
end

function y_predict = get_consensus_prediction(Model, x_test, x_train)
    %% Build predicted variable using all different cross-validated models
    y_predict = [];
    for trained_idx = 1:numel(Model.Trained)
        mdl = Model.Trained{trained_idx};  
        if isempty(x_test) % when using only KFOlds
            temp = predict(mdl, x_train(Model.Partition.test(trained_idx), :));
            y_predict(trained_idx, :) = temp(1:min(Model.Partition.TestSize));
        else % when using held out data
            y_predict(trained_idx, :) = predict(mdl, x_test);
        end
    end
    %figure(11111);clf();plot(y_predict'); hold on;plot(nanmean(y_predict), 'k', 'Linewidth', 2)
    y_predict = nanmean(y_predict)';
    if islogical(Model.Y)
        y_predict = logical(round(y_predict));
    end
end

function [Lmax] = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
    %% Set default Hyperparameter optimization options
    HyperparameterOptimizationOptions = DEFAULT_HYPERPARAMS;    

    if ~parameters.optimize_hyper
    	Lmax = 1e-4;    
    elseif strcmpi(parameters.optimization_method, 'native')
        optimize_varargin   = [base_varargin, {'OptimizeHyperparameters',"Lambda", 'HyperparameterOptimizationOptions', HyperparameterOptimizationOptions}];    
        model               = func(optimize_varargin{:}); 
        Lmax                = model.Lambda;
    elseif strcmpi(parameters.optimization_method, 'manual')
        lrange              = logspace(-4,4,40);
        [score, TPR, TNR, MCC] = deal(NaN(numel(lrange),1));
        base_varargin([find(strcmp(base_varargin, 'Lambda')), find(strcmp(base_varargin, 'Lambda'))+1]) = [];
        
        if isempty(x_test)
            base_varargin = [base_varargin, {'KFold', HyperparameterOptimizationOptions.KFold}];
            Timeout = 3;
        else
            Timeout = 0.5;
        end
        for k = 1:numel(lrange)            
            temp_varargin           = [base_varargin, 'Lambda', lrange(k)];
            fut(k)                  = parfeval(func, 1, temp_varargin{:});
        end 
        for k = 1:numel(lrange)
            [idx,mdl]               = fetchNext(fut, Timeout);
            if ~isempty(mdl)
                if isempty(x_test)
                    temp            = kfoldLoss(mdl, 'Mode','individual', 'LossFun', @pearson_correlation_coefficient)*100; 
                    [score(idx), TPR(idx), TNR(idx), MCC(idx)] = deal(nanmean(temp));
                else
                    y_predict            = mdl.predict(x_test);                
                    [score(idx), TPR(idx), TNR(idx), MCC(idx)] = get_classifier_score(y_test, y_predict);
                end
            end
        end
        
        MCC     = smoothdata(MCC, 'gaussian', 5);
        score   = smoothdata(score, 'gaussian', 5);
        TPR     = smoothdata(TPR, 'gaussian', 5);
        TNR     = smoothdata(TNR, 'gaussian', 5);

        if parameters.rendering >= 3
            figure(1002);clf();subplot(2,2,1);im = plot(MCC, 'XData', lrange);im.Parent.XScale = 'log';title('MCC');xlabel('Lambda');
            subplot(2,2,2);im = plot(score, 'XData', lrange);im.Parent.XScale = 'log';title('Accuracy');xlabel('Lambda');
            subplot(2,2,3);im = plot(TPR, 'XData', lrange);im.Parent.XScale = 'log';title('TPR');xlabel('Lambda');
            subplot(2,2,4);im = plot(TNR, 'XData', lrange);im.Parent.XScale = 'log';title('TNR');xlabel('Lambda');
            drawnow();
        end
        [~, loc] = max(MCC); Lmax = lrange(loc);       
    end 
    
    function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
        score = corr(y_true, y_pred);
    end
    
end



function [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
    if ~parameters.optimize_hyper
    	bmax = 1000; kmax = 20;                 
    elseif strcmpi(parameters.optimization_method, 'manual')
        krange = logspace(-1,4,40);
        brange = logspace(0,6,40);
        [score, TPR, TNR, MCC] = deal(NaN(numel(krange),numel(brange)));
        base_varargin{find(strcmp(base_varargin, 'CrossVal'))+1} = 'off';
        base_varargin([find(strcmp(base_varargin, 'KFold')), find(strcmp(base_varargin, 'KFold'))+1]) = [];
        for k = 1:numel(krange)
            fprintf(['optimization at : ',num2str(100*k/numel(krange)),' %%\n'])
            clear fut
            for b = 1:numel(brange)
                temp_varargin           = [base_varargin, 'KernelScale', krange(k), 'BoxConstraint', brange(b)];
                fut(b)                  = parfeval(func, 1, temp_varargin{:}); 
            end            
            for b = 1:numel(brange)
                [idx,mdl]            = fetchNext(fut, 0.5);
                if ~isempty(mdl)
                    y_predict            = mdl.predict(x_test);                
                    [score(k,idx), TPR(k,idx), TNR(k,idx), MCC(k,idx)] = get_classifier_score(y_test, y_predict);
                end
            end
        end  
        MCC     = imfilter(MCC, fspecial('gaussian', 5, 5), 'replicate');
        score   = imfilter(score, fspecial('gaussian', 5, 5), 'replicate');
        TPR     = imfilter(TPR, fspecial('gaussian', 5, 5), 'replicate');
        TNR     = imfilter(TNR, fspecial('gaussian', 5, 5), 'replicate');

        if parameters.rendering >= 3
            figure(1002);clf();subplot(2,2,1);im = imagesc(MCC, 'XData', brange, 'YData', krange);im.Parent.XScale = 'log';im.Parent.YScale = 'log';im.Parent.YLim = [krange(1),krange(end)];colorbar;title('MCC');xlabel('box');ylabel('kernel')
            subplot(2,2,2);im = imagesc(score, 'XData', brange, 'YData', krange);im.Parent.XScale = 'log';im.Parent.YScale = 'log';im.Parent.YLim = [krange(1),krange(end)];colorbar;title('Accuracy');xlabel('box');ylabel('kernel')
            subplot(2,2,3);im = imagesc(TPR, 'XData', brange, 'YData', krange);im.Parent.XScale = 'log';im.Parent.YScale = 'log';im.Parent.YLim = [krange(1),krange(end)];colorbar;title('TPR');xlabel('box');ylabel('kernel')
            subplot(2,2,4);im = imagesc(TNR, 'XData', brange, 'YData', krange);im.Parent.XScale = 'log';im.Parent.YScale = 'log';im.Parent.YLim = [krange(1),krange(end)];colorbar;title('TNR');xlabel('box');ylabel('kernel')
            drawnow();
        end
        [~, loc] = max(max(MCC));  bmax = brange(loc);
        [~, loc] = max(max(MCC')); kmax = krange(loc);        
    else
        try
            Mdl = fitcsvm(  x_train ,...
                            y_train ,...
                            'KernelFunction'    , 'gaussian', ...
                            'Standardize'       , true, ...
                            'ClassNames'        , [false; true],...
                            'Cost'              , cost,...
                            'OptimizeHyperparameters','auto' ,... 
                            'HyperparameterOptimizationOptions',struct(...
                                                                    'MaxObjectiveEvaluations',150,...
                                                                    'ShowPlots',parameters.rendering >= 3,...
                                                                    'MaxTime',Inf,...
                                                                    'UseParallel',true,...
                                                                    'Repartition',true)); 
            kmax =  Mdl.HyperparameterOptimizationResults.XAtMinObjective.KernelScale;
            bmax =  Mdl.HyperparameterOptimizationResults.XAtMinObjective.BoxConstraint;                                         
        catch
            parameters.optimization_method = 'manual';
            [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters);   
        end
    end              
end