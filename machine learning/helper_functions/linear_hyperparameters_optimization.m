
% xtrain and ytrain are included in base_varargin
function Lmax = linear_hyperparameters_optimization(base_varargin, func, x_test, y_test, parameters)
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
        score               = repmat({},1, numel(lrange));
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
                    if ishandle(LossFun)
                        LossFun = {LossFun};            
                    end
                    score{idx} = {};
                    for el = 1:numel(LossFun)
                        name = erase(func2str(LossFun{el}),'_score');
                        score{idx}.(name) = kfoldLoss(mdl, 'Mode','individual', 'LossFun', LossFun{el});
                        score{idx}.(name) = nanmean(score{idx}.(name));
                    end
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
        score_used = cellfun(@(x) x.(names{end}), score);
        [~, loc] = max(score_used); Lmax = lrange(loc);       
    end     
end
