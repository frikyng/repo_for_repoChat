%% STANDARD SETTINGS FOR DEMO SRIPTS AND ANALYSIS. THIS IS CALLED ACROSS SEVERAL FUNCTIONS

if ~exist('skip_path','var') || (exist('skip_path','var') && ~skip_path)
    %% Load an arboreal_scan_experiment object if required
    path_to_use = 'C:\Users\Antoine.Valera\MATLAB\newextraction_raw_zscored\2019-10-01_exp_1';  
    if ~exist('obj','var')     
        obj = arboreal_scan_experiment(path_to_use);
        path_to_use = obj.source_folder;  
    end
end

%% Arboreal Scan settings
obj.detrend                         = 0;            % -1; %'auto'
obj.rendering                       = 0;
obj.time_smoothing                  = [-1.5, 0];
obj.thr_for_detection               = 0.2;
obj.bad_ROI_thr                     = 0.6;

%% Behaviour Settings
behaviours                          = {'encoder','EyeCam_L_forelimb','EyeCam_R_forelimb','BodyCam_L_whisker','BodyCam_R_whisker','EyeCam_Perioral','trigger'};
obj.beh_smoothing                   = [-1.5, 0];
obj.beh_sm_func                     = 'gaussian';
obj.shuffling_block_size            = -3;
obj.shuffling_block_size            = -3;
obj.detrend_behaviour               = false;

%% Machine Learning settings
use_classifier  = false;                % true or false
method          = 'linear';             % 'linear' or 'svm'    
N_iter          = 20;                   % N iterations
N_iter_for_lag  = 20;                   % N iterations when testing lag
KFold           = 1;                    % N-Kfold. Set to 1 if using  
holdout         = 0.2;                  % percentage of data held out for final testing. If 0, we use KFoldLoss with custom Pearson coeff to estimate the score
rendering       = 1;                    % If 0 no figure, if 1 summary plot is display per iteraton + averages, if 2, individual iteration data is displayed , if 3 hyperparameter optimization plot are displayed.
block_shuffling = -3;                   % size of the blocks to shuffle for ML. values < 0 for block size in s (otherwise in points)
obs_shuf_block_sz= obj.shuffling_block_size;

% Build ml_object
ml_parameters = machine_learning_params(   'optimize_hyper'   , true          ,... % ML varargin ; if true, Lambda is optimized for each iteration
                                           'rendering'        , rendering     ,... % If 0 no figure, if 1 summary plot is display per iteraton + averages, if 2, individual iteration data is displayed , if 3 hyperparameter optimization plot are displayed.
                                           'optimization_method', 'manual'    ,... % "manual" or "native". Native uses MATLAB baysian search parameter optimization for Lambda. Native generate a search grid and find Lambda with the best prediction score.. 
                                           'holdout'          , holdout       ,... % percentage of data held out for final testing. If 0, we use KFoldLoss with custom Pearson coeff to estimate the score
                                           'kfold'            , KFold         ,... %
                                           'block_shuffling'  , block_shuffling,... % in seconds, the time of the block to shuffle (for cross validation));
                                           'savefig'          , path_to_use   ,... % Save figure. If you set a path, the figure is saved there
                                           'save'             , path_to_use   ,...
                                           'N_iter'           , N_iter        ,...
                                           'obs_shuf_block_sz', obs_shuf_block_sz,...
                                           'shuffling'        , 'behaviours')  ;  
                                       
% For lag analysis, define the range and step size
time_step       = 0.25;
time_ext        = 5;
sr              = 1/nanmedian(obj.timescale.sr);
lag_list        = [(0:time_step:time_ext)*-1, (0:time_step:time_ext)];
lag_list        = unique(round(sort(lag_list) * sr));
lag_smoothing   = 0;