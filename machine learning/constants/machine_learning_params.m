function params = machine_learning_params(varargin)    
    params                          = {};
    params.use_classifier           = false ; % If true, a classifier is used instead of a regression model
    params.method                   = 'linear'; % The type of model used for regression / classification
    params.holdout                  = 0.5   ; % This is to seperate the train and test dataset
    params.kFold                    = 5     ; % This is defining the kFold crossvalidation settings for the training dataset
    params.optimize_hyper           = false ; % if true, hyperparameters are optimized during training
    params.optimization_method      = 'native' ; % 'manual' or 'native'. manual is a custom grid search. native is the matlab baysian optimized approach
    params.rendering                = 2     ; % (0) no plot, (1) summary plot, (2) 1 + predictions, (3) 2 + hyperparameter tuning 
    params.svm_kernel               = 'gaussian' ; % gaussian or linear
    params.solver                   = ''    ; % see https://fr.mathworks.com/help/stats/fitrlinear.html, "Solver" section
    params.shuffling                = ''    ; % 'events' to shuffle timepoints, 'ROIs' to shuffle the spatial structure, 'behaviours' to shuffle observations. Use a cell array to cet more than one.
    params.block_shuffling          = 0     ; % If non 0, then we do block shuffling using this window size. holdout must be > 0. kFold must be 1  
    params.title                    = ''    ; % set final bar chart title
    params.alpha                    = []    ; % Set a value between 0 and 1 for elastic Net (1 is lasso and 0 is ridge)
    params.save                     = false ; % If true, the result is saved
    params.saving_level             = 1     ; % level 1 saves the model and its results, level 2 also saves the full predictor and observation traces
    params.savefig                  = false ; % If true, the output figure is saved. If rendering is 0, this set rendering to 1
    params.N_iter                   = 1     ; % If true, the output figure is saved. If rendering is 0, this set rendering to 1
    params.randomize_ROIs           = 0     ; % If 1, ROis are randomized. If you passed groups, random groups will be sized matched. If -1, only ROIs NOT listed in the predictor list will be used (if possible)
    params.add_shuffle              = false ; % If true, behaviour temporal shuffling is computed too 
    params.obs_shuf_block_sz        = 1     ; % when params.shuffling is 'behaviours', this defines the size of the block to shuffle. 1 will shuffle every timepoint, and values > 1 will shuffle blocks
    params.score_metrics            = 'all';  % pearson, mse, rmse, explained_variance, all (for all 4) or a function handle
    params.show_dots                = true  ; % true or false. If true, individual model training are display in the bar chart
    params.split_shuffle            = true  ; % true or false. Set to false to make the charts a bit lighter
    params.use_shuffle              = true  ; % true or false. Set to false to ignore shuffling results in the plots and the stats
    params.show_iterations          = true  ; % true or false. Set to false to hide the dots
    
    %% Unwrap varargin
    while nargin > 0 && iscell(varargin) && iscell(varargin{1})
        varargin = varargin{1};
    end
    
    %% Extract existing params object
    if nargin > 0 && iscell(varargin) && isstruct(varargin{1})
        params = varargin{1};
        varargin(1) = [];
    end
    
    %% Update parameters if required
    if nargin > 0 && ~isempty(varargin) && iscell(varargin) && rem(numel(varargin), 2)
        error('input option must be an even number of cells')
    elseif nargin > 0 && ~isempty(varargin) && iscell(varargin)     
        %% Update default if any
        for i = 1:length(varargin)
            if (strcmpi(varargin{i},'use_classifier'))
                params.use_classifier = lower(varargin{i+1});
            end    
            if (strcmpi(varargin{i},'method'))
                params.method = varargin{i+1};
            end 
            if (strcmpi(varargin{i},'holdout'))
                params.holdout = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'kFold'))
                params.kFold = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'optimize_hyper'))
                params.optimize_hyper = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'optimization_method'))
                params.optimization_method = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'rendering'))
                params.rendering = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'svm_kernel'))
                params.svm_kernel = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'solver'))
                if ischar(varargin{i+1})
                    params.solver = lower(varargin{i+1});
                else
                    params.solver = cellfun(@lower, varargin{i+1}, 'UniformOutput', false); % for multiple solver sequences
                end
            end
            if (strcmpi(varargin{i},'shuffling'))
                params.shuffling = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'block_shuffling'))
                params.block_shuffling = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'title'))
                params.title = varargin{i+1};
            end
            if (strcmpi(varargin{i},'alpha'))
                params.alpha = lower(varargin{i+1});
            end
            if (strcmpi(varargin{i},'save'))
                params.save = varargin{i+1};
            end
            if (strcmpi(varargin{i},'saving_level'))
                params.saving_level = varargin{i+1};
            end
            if (strcmpi(varargin{i},'savefig'))
                params.savefig = varargin{i+1};
            end
            if (strcmpi(varargin{i},'N_iter'))
                params.N_iter = lower(varargin{i+1});
            end 
            if (strcmpi(varargin{i},'randomize_ROIs'))
                params.randomize_ROIs = lower(varargin{i+1});
            end  
            if (strcmpi(varargin{i},'add_shuffle'))
                params.add_shuffle = varargin{i+1};
            end
            if (strcmpi(varargin{i},'obs_shuf_block_sz'))
                params.obs_shuf_block_sz = varargin{i+1};
            end     
            if (strcmpi(varargin{i},'score_metrics'))
                params.score_metrics = varargin{i+1};
            end  
            if (strcmpi(varargin{i},'split_shuffle'))
                params.split_shuffle = varargin{i+1};
            end  
            if (strcmpi(varargin{i},'show_dots'))
                params.show_dots = varargin{i+1};
            end  
            if (strcmpi(varargin{i},'use_shuffle'))
                params.use_shuffle = varargin{i+1};
            end               
            if (strcmpi(varargin{i},'show_iterations'))
                params.show_iterations = varargin{i+1};
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
    if params.obs_shuf_block_sz < 1
        error('Minimal value for obs_shuf_block_sz is 0. Shuffling is only used if params.shuffling is "behaviours"');
    end
end
