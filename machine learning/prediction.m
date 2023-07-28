function [y_predict, y_test, score, x_test, x_train, y_train, model] = prediction(XData, YData, partition, cost, ml_parameters)
    if nargin < 5 || isempty(ml_parameters)
        ml_parameters = machine_learning_params;
    else
        ml_parameters = machine_learning_params(ml_parameters);
    end

    %% Define timepoints for train and test
    if isempty(partition) % if we just do many kfolds
        train_tp    = 1:numel(YData);
        test_tp     = [];
    elseif isstruct(partition)
    	train_tp    = partition.x_train;
        test_tp     = partition.x_test;
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
    if strcmpi(ml_parameters.method, 'svm') 
        %% https://fr.mathworks.com/help/stats/support-vector-machine-classification.html
        base_varargin = {   x_train         , ...
                            y_train         , ...
                            'KernelFunction', ml_parameters.svm_kernel, ...
                            'Standardize'   , true, ...
                            'CrossVal'      , 'on'};
        
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
        [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, ml_parameters);   
        base_varargin = [base_varargin, 'KernelScale', kmax, 'BoxConstraint', bmax];
        if ml_parameters.kFold > 1
            base_varargin = [base_varargin, {'KFold', ml_parameters.kFold}]; 
        end

        %% Train models
        model         = func(base_varargin{:});
        y_predict     = get_consensus_prediction(model, x_test, x_train); %x_train only when not using held out data 

        if islogical(y_train)
            y_predict = logical(round(y_predict));
        end
    
        score = [];
    elseif strcmpi(ml_parameters.method, 'linear')      
        if ml_parameters.use_classifier
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , ml_parameters.postfit_bias,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'logistic'}; 
        elseif isempty(ml_parameters.alpha) % this should be no regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , ml_parameters.postfit_bias,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares'}; % default is no regularization, no standardization (e.g. data is not mean centered)
        elseif ml_parameters.alpha == 0 % ridge regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , ml_parameters.postfit_bias,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares',...
                            'Regularization', 'ridge'}; 
        elseif ml_parameters.alpha == 1 % lasso regularization
            base_varargin = {   x_train         , ...
                            y_train         , ...                              
                            'PostFitBias'   , ml_parameters.postfit_bias,...
                            'PassLimit'     , 10,...
                            'Learner'       , 'leastsquares',...
                            'Regularization', 'lasso'};
        else
            error('regularization method not implemented')
        end
                              
        %% If you manualy set the solver, adjust it here
        if ~isempty(ml_parameters.solver)                
            base_varargin = [base_varargin, 'Solver', ml_parameters.solver];   
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

        Lambda        = linear_hyperparameters_optimization(base_varargin, func, x_test, y_test, ml_parameters);   
        base_varargin = [base_varargin, {'Lambda', Lambda}];
        if ml_parameters.kFold > 1
            base_varargin = [base_varargin, {'KFold', ml_parameters.kFold}]; 
        end
           
        %% Train models
        model         = func(base_varargin{:});
        
        %% Get predictove score
        if ml_parameters.kFold > 1
            y_predict     = get_consensus_prediction(model, x_test, x_train); %x_train only when not using held out data 
        elseif ~isempty(x_test)
            y_predict     = model.predict(x_test);
        else
            y_predict     = model.predict(x_train);
        end
        
        if islogical(y_train)
            y_predict = logical(round(y_predict));
        end
    
        score = [];
    elseif strcmpi(ml_parameters.method, 'glm')  
        % Set the base hyperparameters
        base_varargin = {x_train, ...
                         y_train, ...
                         'Distribution', 'normal', ... % assuming a normal distribution as in linear regression
                         'Link', 'identity'}; % 'identity' link function for linear regression


                     
        % If a solver has been manually set, adjust it here
        if ~isempty(ml_parameters.solver)                
            base_varargin = [base_varargin, 'Solver', ml_parameters.solver];   
        end

        % Set the training model function
        func = @fitglm;

        % Lambda is the equivalent to the regularization strength in ridge regression. For 'fitglm', 
        % it doesn't directly take Lambda. However, you can include binomial size, weights, offset, 
        % etc. based on your data and problem characteristics. The line below is a placeholder. 
        % Please adjust accordingly.
        
        if ml_parameters.optimize_hyper || strcmp(ml_parameters.optimization_method, 'manual') || strcmp(ml_parameters.optimization_method, 'native')
            Lambda = logspace(-3, 0, 10);%logspace(-5, 5, 100);
            % Fit a generalized linear model with Lasso regularization
            [B,FitInfo] = lassoglm(x_train,y_train,'normal','Link','identity','Lambda',Lambda,'CV',5);
            lassoPlot(B,FitInfo,'plottype','CV','XScale','log'); 
            legend('show'); drawnow % Show legend
            optimalLambda  = FitInfo.Lambda(FitInfo.IndexMinDeviance);
            
            % Fit the final model using the optimal Lambda
            [B_final, FitInfo_final] = lassoglm(x_train, y_train, 'normal', 'Link', 'identity', 'Lambda', optimalLambda);
            
            % Define a custom predict function for the lassoglm model
            model = @(X) [ones(size(X, 1), 1), X] * [FitInfo_final.Intercept; B_final];

            % Get the predictive score
            if ml_parameters.kFold > 1
                y_predict = get_consensus_prediction(model, x_test, x_train); % Pass the custom predict function to get_consensus_prediction
            elseif ~isempty(x_test)
                y_predict = model(x_test);
            else
                y_predict = model(x_train);
            end 
        else
            % If you're using k-fold cross-validation, include the 'KFold' parameter.
            if ml_parameters.kFold > 1
                base_varargin = [base_varargin, {'KFold', ml_parameters.kFold}];
            end

            % Train the model
            model = func(base_varargin{:});
            
            % Get the predictive score
            if ml_parameters.kFold > 1
                y_predict = get_consensus_prediction(model, x_test, x_train); % x_train only when not using held-out data 
            elseif ~isempty(x_test)
                y_predict = model.predict(x_test);
            else
                y_predict = model.predict(x_train);
            end
        end

        % Initialize score variable
        score = [];

    else
        error('classification or regression method not implemented')
    end
    
    
    
    if ~iscolumn(y_test)
        y_test = y_test';
    end      
    
    
    
    
end

function y_predict = get_consensus_prediction(ModelOrPredictFunc, x_test, x_train)
    %% Build predicted variable using all different cross-validated models
    y_predict = [];
    
    if isa(ModelOrPredictFunc, 'function_handle') % Custom predict function for lassoglm
        if isempty(x_test) % when using only KFOlds
            y_predict = ModelOrPredictFunc(x_train);
        else % when using held out data
            y_predict = ModelOrPredictFunc(x_test);
        end
    else % Model object
        for trained_idx = 1:numel(ModelOrPredictFunc.Trained)
            mdl = ModelOrPredictFunc.Trained{trained_idx};  
            if isempty(x_test) % when using only KFOlds
                temp = predict(mdl, x_train(ModelOrPredictFunc.Partition.test(trained_idx), :));
                y_predict(trained_idx, :) = temp(1:min(ModelOrPredictFunc.Partition.TestSize));
            else % when using held out data
                y_predict(trained_idx, :) = predict(mdl, x_test);
            end
        end
    end
    %figure(11111);clf();plot(y_predict'); hold on;plot(nanmean(y_predict), 'k', 'Linewidth', 2)
    
    if islogical(ModelOrPredictFunc.Y)
        y_predict = nanmedian(y_predict,1);
    else
        y_predict = nanmean(y_predict)';
    end
    
    if ~iscolumn(y_predict)
        y_predict = y_predict';
    end    
end





function [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
    error('to fix with the new flexible scioring method')
    if ~parameters.optimize_hyper
    	bmax = 1000; kmax = 20;                 
    elseif strcmpi(parameters.optimization_method, 'manual')
        krange = logspace(-1,4,10);
        brange = logspace(0,6,10);
        [score, TPR, TNR, MCC] = deal(repmat({},numel(krange),numel(brange)));
        base_varargin{find(strcmp(base_varargin, 'CrossVal'))+1} = 'off';
        base_varargin([find(strcmp(base_varargin, 'KFold')), find(strcmp(base_varargin, 'KFold'))+1]) = [];
        for k = 1:numel(krange)
            %fprintf(['optimization at : ',num2str(100*k/numel(krange)),' %%\n'])
            clear fut
            for b = 1:numel(brange)
                temp_varargin           = [base_varargin, 'KernelScale', krange(k), 'BoxConstraint', brange(b)];
                fut(b)                  = parfeval(func, 1, temp_varargin{:}); 
            end            
            for b = 1:numel(brange)
                [idx,mdl]            = fetchNext(fut, 0.5);
                if ~isempty(mdl)
                    y_predict            = mdl.predict(x_test);                
                    score(k,idx) = get_ml_score(y_test, y_predict);
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
        [~, loc]    = max(MCC(:));
        [kmax,bmax] = ind2sub(size(MCC),loc);
        bmax        = brange(bmax);
        kmax        = krange(kmax);        
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