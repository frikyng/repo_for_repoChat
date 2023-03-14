function params = DEFAULT_CLASSIFIER_OPTION(varargin)    
    params                          = {};
    params.holdout                  = 0.5   ; % This is to seperate the train and test dataset
    params.kFold                    = 5     ; % This is defining the kFold crossvalidation settings for the training dataset
    params.optimize_hyper           = false ; % if true, hyperparameters are optimized during training
    params.optimization_method      = 'native' ; % 'manual' or 'native'. manual is a custom grid search. native is the matlab baysian optimized approach
    params.rendering                = 2     ; % (0) no plot, (1) summary plot, (2) 1 + predictions, (3) 2 + hyperparameter tuning 
    params.svm_kernel               = 'gaussian' ; % gaussian or linear
    params.solver                   = ''    ; % see https://fr.mathworks.com/help/stats/fitrlinear.html, "Solver" section
    params.shuffling                = ''    ; % set events to shuffle timpoints, and ROIs to shuffle the spatial structure, and both to do all
    params.block_shuffling          = 0     ; % If non 0, then we do block shuffling using this window size. holdout must be > 0. kFold must be 1  
    params.title                    = ''    ; % set final bar chart title
    params.alpha                    = []    ; % Set a value between 0 and 1 for elastic Net (1 is lasso and 0 is ridge)
    params.save                     = false ; % If true, the result is saved
    params.savefig                  = false ; % If true, the output figure is saved. If rendering is 0, this set rendering to 1
    
    %% Unwrap varargin
    while nargin > 0 && iscell(varargin) && iscell(varargin{1})
        varargin = varargin{1};
    end
    
    %% Update parameters if required
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
            if(strcmpi(varargin{i},'solver'))
                if ischar(varargin{i+1})
                    params.solver = lower(varargin{i+1});
                else
                    params.solver = cellfun(@lower, varargin{i+1}, 'UniformOutput', false); % for multiple solver sequences
                end
            end
            if(strcmpi(varargin{i},'shuffling'))
                params.shuffling = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'title'))
                params.title = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'alpha'))
                params.alpha = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'save'))
                params.save = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'savefig'))
                params.savefig = lower(varargin{i+1});
            end
            if(strcmpi(varargin{i},'block_shuffling'))
                params.block_shuffling = lower(varargin{i+1});
            end
        end
    end
    
    if any(params.savefig) && ~params.rendering
        params.rendering = 1;
    end
    if params.block_shuffling && params.kFold > 1
        warning('params.kFold > 1 not supported for block shuffling. Kfold set to 1')
        params.kFold = 1;
    end  
    if (params.block_shuffling || params.kFold == 1) && params.holdout == 0
        error('When using block_shuffling or KFold == 1, you must specify holdout > 0')
    end    
end
