function params = DEFAULT_CLASSIFIER_OPTION(varargin)    
    params                          = {};
    params.holdout                  = 0.2   ; % This is to seperate the train and test dataset
    params.kFold                    = 5     ; % This is defining the crossvalidation settings for the training dataset
    params.optimize_hyper           = false ; % if true, hyperparameters are optimized during training
    params.optimization_method      = 'native' ; % 'manual' or 'native'
    params.rendering                = 2     ; % (0) no plot, (1) summary plot, (2) 1 + predictions, (3) 2 + hyperparameter tuning 
    params.svm_kernel               = 'gaussian' ; % gaussian or linear
    
    %% Unwrap varargin
    while nargin > 0 && iscell(varargin) && iscell(varargin{1})
        varargin = varargin{1};
    end
    
    %% Update aparameters if required
    if nargin > 0 && iscell(varargin) && ~isstruct(varargin{1}) && rem(numel(varargin), 2)
        error('input option must be an even number of cells')
    elseif nargin > 0 && isstruct(varargin{1})
        params = varargin{1};
    elseif nargin > 0        
        %% Update default if any
        for i = 1:length(varargin)
            if(strcmpi(varargin{i},'holdout'))
               params.holdout = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'kFold'))
               params.kFold = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'optimize_hyper'))
               params.optimize_hyper = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'rendering'))
               params.rendering = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'svm_kernel'))
               params.svm_kernel = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'optimization_method'))
               params.optimization_method = lower(varargin{i+1});
            end
        end
    end
end
