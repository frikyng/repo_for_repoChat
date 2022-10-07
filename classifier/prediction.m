%% Need to discuss 
%% 'Solver' — Optimization routine --> 'ISDA' | 'L1QP' | 'SMO'
%% 'Weights' — Observation weights --> numeric vector | name of variable in Tbl % should we use SNR?

function [y_predict, y_test, score, x_test, x_train, y_train, model] = prediction(XData, YData, partition, method, cost, parameters)
    if nargin < 6 || isempty(parameters)
        parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        parameters = DEFAULT_CLASSIFIER_OPTION(parameters);
    end

    %% Define timepoints for train and test
    train_tp    = find(training(partition));
    test_tp     = find(test(partition));

    %% Assign predictor and observation data to test and train dataset
    x_train     = double(XData(:,train_tp))';       % Predictors for training
    x_test      = double(XData(:,test_tp))';        % Predictors for testing
    y_train     = YData(:, train_tp);               % Variable to predict used for training  
    y_test      = YData(:, test_tp);                % Variable to predict used for testing (ground truth)     

    %% Set default Hyperparameter optimization options
    HyperparameterOptimizationOptions = DEFAULT_HYPERPARAMS;
    
    %% Run the correst classifier or regression
    if strcmp(method, 'svm') 
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
        fprintf(['The model out-of-sample misclassification rate is ',num2str(score,3),'%%\n']);  
    elseif strcmp(method, 'linear')        
        base_varargin = {   x_train         , ...
                            y_train         , ...  
                            'Learner'       , 'leastsquares',...
                            'Regularization','ridge',...
                            'PostFitBias'   , true,...
                            'PassLimit'     , 10};

        %% Set training model
        if islogical(y_train)
            func = @fitclinear;
            base_varargin = [base_varargin, {'Cost', cost, 'ClassNames', [false; true]}];
        else
            updated_varargin = base_varargin;
            func = @fitrlinear;
        end

        Lambda        = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters);   
        base_varargin = [base_varargin, {'Lambda', Lambda}];
        base_varargin = [base_varargin, {'KFold', 5}]; 
        
        %% Train models
        model         = func(base_varargin{:});

        y_predict = get_consensus_prediction(model, x_test); 
        if islogical(y_train)
            y_predict = logical(round(y_predict));
        end

        %% Estimate accuracy of the cross-validated model, and of the predicted values vs ground thruth
        score = kfoldLoss(model)*100; % reveals the fraction of predictions that were incorrect, i.e. (1 - accuracy)
        fprintf(['The out-of-sample misclassification rate is ',num2str(score,3),'%%\n']);  

    
    elseif strcmp(method, 'forest')
        classificationTree = TreeBagger(500, x_train, y_train, 'OOBPrediction', 'on','Method', 'classification','Cost', cost); 
        y_predict = str2double(classificationTree.predict(x_test));  
    elseif strcmp(method, 'glm')  
            mdl =  fitglm(x_train,y_train,'linear','Distribution','normal');
            [y_predict,~] = predict(mdl,x_test)
            mask = ones(size(mdl.Coefficients.SE)-1);
            mask = abs(mdl.Coefficients.SE(2:end)) < 0.001;
            mdl2 = removeTerms(mdl, mask')
            score  = NaN
            [y_predict2,~] = predict(mdl2,x_test)
            assignin('base','p',mdl.Coefficients.pValue)
            assignin('base','SE',mdl.Coefficients.SE)
            assignin('base','tStat',mdl.Coefficients.tStat)
            assignin('base','Estimate',mdl.Coefficients.Estimate)
            %obj.ref.plot_value_tree(SE, find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)))
            
%             [y_predict,~] = predict(mdl,x_test)
%             score  = NaN 
%             
%             

    else
        error('classification method not implemented')
    end
end

function y_predict = get_consensus_prediction(Model, x_test)
    %% Build predicted variable using all different cross-validated models
    y_predict = [];
    for trained_idx = 1:numel(Model.Trained)
        mdl = Model.Trained{trained_idx};                        
        y_predict(trained_idx, :) = predict(mdl, x_test);
    end
    y_predict = nanmean(y_predict)';
    if islogical(Model.Y)
        y_predict = logical(round(y_predict));
    end
end

function [Lmax] = linear_hyperparameters_optimization(x_train, y_train, x_test, y_test, cost, func, base_varargin, parameters)
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
        for k = 1:numel(lrange)            
            temp_varargin           = [base_varargin, 'Lambda', lrange(k)];
            fut(k)                  = parfeval(func, 1, temp_varargin{:});
        end  
        for k = 1:numel(lrange)
            [idx,mdl]            = fetchNext(fut, 0.5);
            if ~isempty(mdl)
                y_predict            = mdl.predict(x_test);                
                [score(idx), TPR(idx), TNR(idx), MCC(idx)] = get_classifier_score(y_test, y_predict);
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