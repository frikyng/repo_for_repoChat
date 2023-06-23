%% Optimizes linear hyperparameters.
%  Performs linear hyperparameters optimization based on given parameters and
%  score metrics. The optimization can be either manual or native. It uses
%  machine learning to find optimal values.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	Lmax = linear_hyperparameters_optimization(base_varargin, func, x_test, 
%                                              y_test, parameters)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	base_varargin({Cell} of various object types) - Optional:
%                                   Cell array with all base variable 
%                                   arguments needed for the function.
%
% 	func(Function Handle):
%                                   Function handle for the function to be 
%                                   optimized.
%
% 	x_test(Matrix):
%                                   Test data set used for optimization.
%
% 	y_test(Matrix):
%                                   Test labels corresponding to the test
%                                   data set.
%
% 	parameters(Struct):
%                                   Structure with fields controlling the 
%                                   optimization, including 'score_metrics', 
%                                   'optimize_hyper' and 'optimization_method'.
%                                   Each field can take on various values and 
%                                   control different aspects of the 
%                                   optimization process.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	Lmax(Numeric):
%                                   The optimal value of the hyperparameter 
%                                   Lambda.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * In the case of manual optimization, this function will iterate over a
%   range of Lambda values, storing the corresponding scores for each. The
%   function will then select the Lambda value that yields the highest score.
%
% -------------------------------------------------------------------------
%% Examples:
% * Optimizing linear hyperparameters:
% 	Lmax = linear_hyperparameters_optimization(base_varargin, func, x_test, 
%                                              y_test, parameters);
%
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% This function was initially released as part of The SilverLab MatLab
% Imaging Software, an open-source application for controlling an
% Acousto-Optic Lens laser scanning microscope. The software was 
% developed in the laboratory of Prof Robin Angus Silver at University
% College London with funds from the NIH, ERC and Wellcome Trust.
%
% Copyright © 2015-2020 University College London
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License. 
% -------------------------------------------------------------------------
% Revision Date:
% 	09-06-2023
% -------------------------------------------------------------------------
% See also: 
%	machine_learning_hyper_params, pearson_correlation_coefficient,
%   rmse_score, mse_score, explained_variance_score

% TODO : Improve the rendering for the new flexible scoring method
% 	      Handle errors related to unrecognized score function

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
