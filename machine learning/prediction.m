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

        Lambda        = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, ml_parameters);   
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
        
%     %% still linear regression but we need to use lasso function and fitrlinear doesn't have alpha as parameter
%     elseif strcmpi(ml_parameters.method, 'elastic')
%        error('to review')
%        %convert roi_subset into array of strings so lasso is happy
%        %for alpha = [0.1,0.9]
%        roi_subset = [1:111];
%        roi_subset_cell = num2cell(roi_subset);
%        list_ROIs_lasso =cellfun(@num2str,roi_subset_cell,'un',0);
%        
%        % run lasso, each value in the columns of the first output corresponds to a particular regularization
%        % coefficient in Lambda, using a geometric sequence of Lambdas 
%        [rglrzd_Lambda_coefs, FitInfo] = lasso( x_train         , ...
%                                           y_train         , ... 
%                                           'CV', 10        , ...
%                                           'PredictorNames', list_ROIs_lasso,...
%                                           'Alpha', 0.5);   %...       % params.alpha or use 'Alpha', alpha, where alpha can be between 0 and 1
%                                       
% 
%         % Lasso plot with x-validated fits
%         % look at this plot across a range of alphas
%         %lassoPlot(rglrzd_Lambda_coefs,FitInfo,'PlotType','CV'); legend('show');
%         
%         % Display variable names (or ROI list) in a model that corresponds to minimum x-validated mean squared error
%         idxLambdaMinMSE = FitInfo.IndexMinMSE;
%         minMSEModelPredictors = FitInfo.PredictorNames(rglrzd_Lambda_coefs(:,idxLambdaMinMSE)~=0);
%         
%         % get betas that correspond to optimal model (from MinMSE of whatever lambda)
%         coef_minMSE = rglrzd_Lambda_coefs(:,idxLambdaMinMSE);
%         coef0_minMSE = FitInfo.Intercept(idxLambdaMinMSE);
%         
%         y_predict = x_test*coef_minMSE + coef0_minMSE;
%         score = NaN;
%         model = NaN;
%         
% %         hold on
% %         scatter(y_test,y_predict,'red');scatter(x_test,y_test,'blue')
% %         plot(y_test,y_test)
% %         plot((x_test*coef_minMSE + coef0_minMSE), y_predict)
% %         xlabel('Test')
% %         ylabel('Predict')
% %         hold off
%                     
%         % Display variable names (or ROI list) in sparsest model within one SE of minimum MSE
%         % (this 1SE from the minMSE is a heuristic that is widely accepted by the lasso community and accords with parsimony)
%         % but this heuristic doesn't work in this cell.. just gives empty array
% %         idxLambda1SE = FitInfo.Index1SE;
% %         sparseModelPredictors = FitInfo.PredictorNames(rglrzd_Lambda_coefs(:,idxLambda1SE)~=0);
%         
%          % get betas that correspond to optimal model (from 1MSE of whatever lambda)
% %         coef_1SE = rglrzd_Lambda_coefs(:,idxLambda1SE);
% %         coef0_1SE = FitInfo.Intercept(idxLambda1SE);
%         
% %         y_predict = x_test*coef_1SE + coef0_1SE;
% %         hold on
% %         scatter(y_test,y_predict)
% %         plot(y_test,y_test)
% %         xlabel('Test')
% %         ylabel('Predict')
% %         hold off
%         
%      %%    
%     elseif strcmpi(ml_parameters.method, 'forest')
%         classificationTree = TreeBagger(500, x_train, y_train, 'OOBPrediction', 'on','parameters.method', 'classification','Cost', cost); 
%         y_predict = str2double(classificationTree.predict(x_test));  
%     elseif strcmpi(ml_parameters.method, 'glm')  
%             model =  fitglm(x_train,y_train,'linear','Distribution','inverse gaussian'); %  'binomial' | 'poisson' | 'gamma' | 'inverse gaussian'
%             [y_predict,~] = predict(model,x_test);
%             score = NaN;
%             % perform cross-validation, and return average MSE across folds
%             %mse = crossval('mse', x_train,y_train', 'Predfun',fcn, 'kfold',10);
%             % compute root mean squared error
%             %avrg_rmse = sqrt(mse)
% 
% %             mask = ones(size(mdl.Coefficients.SE)-1);
% %             mask = abs(mdl.Coefficients.SE(2:end)) < 0.001;
% %             mdl2 = removeTerms(mdl, mask')
% %             score  = NaN
% %             [y_predict2,~] = predict(mdl2,x_test)
% %             assignin('base','p',mdl.Coefficients.pValue)
% %             assignin('base','SE',mdl.Coefficients.SE)
% %             assignin('base','tStat',mdl.Coefficients.tStat)
% %             assignin('base','Estimate',mdl.Coefficients.Estimate)
%             %obj.ref.plot_value_tree(SE, find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)))
%             
% %             [y_predict,~] = predict(mdl,x_test)
% %             score  = NaN 
% %             
% %             

    else
        error('classification or regression method not implemented')
    end
    
    
    
    if ~iscolumn(y_test)
        y_test = y_test';
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
    
    if islogical(Model.Y)
        y_predict = nanmedian(y_predict,1);
    else
        y_predict = nanmean(y_predict)';
    end
    
    if ~iscolumn(y_predict)
        y_predict = y_predict';
    end    
end

function Lmax = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
    %% Set default Hyperparameter optimization options
    HyperparameterOptimizationOptions = machine_learning_hyper_params; 

    if strcmpi(parameters.score_metrics, 'pearson')
        LossFun = @pearson_correlation_coefficient;
    elseif strcmpi(parameters.score_metrics, 'rmse')
        LossFun = @rmse_score;
    elseif strcmpi(parameters.score_metrics, 'mse')
        LossFun = @mse_score;
    elseif strcmpi(parameters.score_metrics, 'explained_variance')
        LossFun = @explained_variance_score;
    elseif ishandle(parameters.score_metrics)
        LossFun = ml_parameters.score_metrics;        
    elseif strcmpi(parameters.score_metrics, 'all')
        LossFun = {@pearson_correlation_coefficient, @mse_score, @rmse_score, @explained_variance_score};
    else
        error('score function not recognized. Use "pearson", "mse", "rmse" or a valid function handle (see pearson_correlation_coefficient.m as an example)')
    end

    if ~parameters.optimize_hyper
    	Lmax = 1e-4;    
    elseif strcmpi(parameters.optimization_method, 'native')
        optimize_varargin   = [base_varargin, {'OptimizeHyperparameters',"Lambda", 'HyperparameterOptimizationOptions', HyperparameterOptimizationOptions}];    
        model               = func(optimize_varargin{:}); 
        Lmax                = model.Lambda;
    elseif strcmpi(parameters.optimization_method, 'manual')
        lrange              = logspace(-4,4,40);
        score               = num2cell(NaN(numel(lrange),1));
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
                    score{idx}      = kfoldLoss(mdl, 'Mode','individual', 'LossFun', LossFun)*100; 
                else
                    y_predict       = mdl.predict(x_test);                
                    score{idx}      = get_ml_score(y_test, y_predict, '' ,LossFun);
                end
            end
        end
        
        names = fieldnames(score{1})';
        for name = names
            name = name{1};
            for iter = 1:numel(score)
                score{iter}.(name)   = smoothdata(score{iter}.(name), 'gaussian', 5);
            end
        end

        if parameters.rendering >= 3
            error('to update with nexttile and score(idx).(name) to handle the new flexible sciring method')
            %             figure(1002);clf();
            %             
            %             subplot(2,2,1);im = plot(MCC, 'XData', lrange);im.Parent.XScale = 'log';title('MCC');xlabel('Lambda');
            %             subplot(2,2,2);im = plot(score, 'XData', lrange);im.Parent.XScale = 'log';title('Accuracy');xlabel('Lambda');
            %             subplot(2,2,3);im = plot(TPR, 'XData', lrange);im.Parent.XScale = 'log';title('TPR');xlabel('Lambda');
            %             subplot(2,2,4);im = plot(TNR, 'XData', lrange);im.Parent.XScale = 'log';title('TNR');xlabel('Lambda');
            drawnow();
        end
        score_used = cellfun(@(x) x.(names{1}), score);
        [~, loc] = max(score_used); Lmax = lrange(loc);       
    end     
end



function [bmax, kmax] = svm_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
    error('to fix with the new flexible scioring method')
    if ~parameters.optimize_hyper
    	bmax = 1000; kmax = 20;                 
    elseif strcmpi(parameters.optimization_method, 'manual')
        krange = logspace(-1,4,10);
        brange = logspace(0,6,10);
        [score, TPR, TNR, MCC] = deal(NaN(numel(krange),numel(brange)));
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