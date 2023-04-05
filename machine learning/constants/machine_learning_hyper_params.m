function params = machine_learning_hyper_params(varargin)    
    params = struct(...
                    'MaxObjectiveEvaluations'   ,150    ,...
                    'ShowPlots'                 ,false  ,...
                    'MaxTime'                   ,Inf    ,...
                    'KFold'                     ,5      , ...
                    'UseParallel'               ,true   ,...
                    'Verbose'                   ,0      ,...
                    'Repartition'               ,true	,...
                    'Optimizer'                 ,'gridsearch'); % bayesopt
end