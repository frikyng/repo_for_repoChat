%% Extract and display results from a machine elarning ouput
% -------------------------------------------------------------------------
%% Syntax:
%   [all_scores, valid_trace_idx]  = extract_ML_coefs(obj,
%           training_output, rendering, metric, cmap)
%
% -------------------------------------------------------------------------
%% Inputs:
%   obj(arboreal_scan_experiment object, or CHAR) - Optional - Default is
%           empty
%       obj can be :
%       * An arboreal_scan_experiment object (Recommended)
%       * A STR PATH to an extracted arboreal_scan_experiment
%       * an empty value, in whic case it will look for the current folder
%
%   training_output(1xN CELL ARRAY of MACHINE LEARNING OUPTUT STRUCTURE) :
%       * Structure containing one cell per model training, as generated by
%       the predict_behaviour function
%
%   rendering(char) - Optional - Default is 0
%       * If 1, a tree with the selected metrics is display.
%       * If 2, a plot of the model score vs number of ROIs is generated
%
%   metric(CHAR) - Optional - Default is Beta 
%       * One of the valid score metrics.
%           - beta will plot the beta coef from the model
%           - 'explained variance' will display the individual explained
%               variance for each ROI
%           - 'order' will display the ranking of each ROi based on they
%               contribution to the model using the beta coef
%
%   cmap(char) - Optional - Default is 'redblue'
%       The colormap applied to the tree
% -------------------------------------------------------------------------
%% Outputs:
%
%	all_scores (1xN CELL ARRAY of 1xM VALUES) :
%       * For each N behaviour, a 1xM score with M being the number of
%       valid ROIs (see valid_trace_idx to know which ones)
%
%   valid_trace_idx(1 x N INT)
%       The ROIs used for precition
%
% -------------------------------------------------------------------------
%% Extra Notes:
% -------------------------------------------------------------------------
%% Examples:
%
% * Plot explained_variance chart from a previously computed model
%   expe_path = 'C:\Users\Antoine.Valera\MATLAB\newextraction_raw_zscored\2019-10-01_exp_1'
%   % load previous model
%   model_fname = 'Prediction on all behaviours using events modulation on all ROIs.mat'
%   load([expe_path, '\', model_fname]);
%   % load arboreal_experiment_object
%   obj = arboreal_scan_experiment(expe_path);
%   % plot explained variance per ROI
%   extract_ML_coefs(obj, results, true, 'explained_variance');
%   
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera  
% -------------------------------------------------------------------------
%% Revision Date:
%   17-04-2023
% -------------------------------------------------------------------------
%% See also:
%    


function [all_scores, valid_trace_idx]  = extract_ML_coefs(obj, training_output, rendering, metric, cmap)
    if ~isa(obj, 'arboreal_scan_experiment')
        error('Input must be an extarcted arboreal_Scan_experiment object')
    end
    if nargin < 3 || isempty(rendering)
        rendering = true;
    end
    if nargin < 4 || isempty(metric)
        metric = 'beta';
    end
    if nargin < 5 || isempty(cmap)
        cmap = 'redblue';
    end

    %ROIs            = training_output{1}.used_ROIs; % Note that in the case of randomized shuffle, this actually changes from one training to the next
    labels          = fix_labels(training_output{1}.beh_type);
    valid_trace_idx = training_output{1}.used_ROIs;    
    behaviours      = training_output{1}.beh_type;
    all_scores      = {};
    not_shuffled_idx= find(~contains(training_output{1}.beh_type, '_shuffled'));
    
    %% If you use explained variance, we need to rextract the predictors and observations if they were not saved.
    ROI_list                        = cell2mat(training_output{1}.used_ROIs); % List of ROis used during training
    
    %% Now, get value for each type of observation
    for beh_idx = not_shuffled_idx%(4)
        all_scores{beh_idx}   = [];   
        ROI_list_sorted = [];
        for iter = 1:numel(training_output)            
            model           = training_output{iter}.model{beh_idx};
            if isprop(model, 'Trained')
                for trained_idx = 1:numel(model.Trained)
                    mdl = model.Trained{trained_idx}; 
                    if strcmpi(metric, 'beta')
                        all_scores{beh_idx}(trained_idx, :) = mdl.Beta;
                    elseif strcmpi(metric, 'order')
                        [~, all_scores{beh_idx}(trained_idx, :)] = sort(abs(model.Beta));  
                    end
                end
            elseif isempty(model) % empty if the behaviour was missing
                continue
            else
                if strcmpi(metric, 'beta')
                    all_scores{beh_idx}(iter, :)            = model.Beta;     
                elseif strcmpi(metric, 'order')
                    [~, all_scores{beh_idx}(iter, :)]       = sort(abs(model.Beta));    
                elseif contains(metric, 'explained_variance')
                    true_behaviours                         = training_output{iter}.observation{beh_idx};
                    predictors                              = training_output{iter}.predictors;
                    predictions                             = model.Bias + model.Beta .* predictors;
                    for roi = 1:numel(ROI_list)
                        all_scores{beh_idx}(iter, roi)      = explained_variance_score(predictions(roi, :), true_behaviours, '');
                    end                  
                    if strcmpi(metric, 'cum_explained_variance')    
                        [~, ordered_beta]                       = sort(abs(model.Beta), 'descend');
                        beta_sorted                             = model.Beta(ordered_beta);
                        ROI_list_sorted(iter, :)                = ordered_beta;
                        %ROI_list_sorted                        = ROI_list(order);
                        n_coefs                                 = 1:numel(model.Beta);
                        for coef = n_coefs
                            current_pred = model.Bias;
                            for i = 1:coef
                                current_pred = current_pred + beta_sorted(i) * predictors(ROI_list_sorted(iter, i),:);
                            end  
                            expl_var(iter, coef) = explained_variance_score(current_pred, true_behaviours, '');                        
                        end
                    end
                elseif strcmpi(metric, 'partialcorr')
                    true_behaviours                         = training_output{iter}.observation{beh_idx};
                    predictors                              = training_output{iter}.predictors;
                    all_scores{beh_idx}(iter, :)            = semi_partial_correlation_coefficients(true_behaviours, predictors, training_output{iter}.ml_parameters);
                end
            end
        end
        
        %% Plot learning performance
        if rendering && strcmpi(metric, 'cum_explained_variance')
            added = NaN(size(all_scores{1}));
            for iter = 1:numel(training_output) 
                % get current variance values
                current = expl_var(iter, :);
                % get the added variance when each ROi is added to the model
                current_diff = diff([0, current]);
                % set the added variance to the correct ROI # (because they
                % were sorted by beta coef until now)
                added(iter, ROI_list_sorted(iter, :)) = current_diff;
            end
            
            
            figure();plot(n_coefs, expl_var, 'Color', [0.8,0.8,0.8]); hold on
            plot(n_coefs, nanmean(expl_var, 1),'ko-','MarkerFaceColor', 'k'); hold on;ylim([0,75])   
            set(gca,'box','off');hold on;xlabel('number of ROI');ylabel('Explained Variance');set(gcf, 'Color', 'w');title(['Explained variance for ',labels{beh_idx},' with increasing number of ROIs']);
            all_scores{beh_idx} = nanmedian(added,1);
        
        end
        
%             %all_scores{beh_idx} = expl_var;            
%             next_fig
%             [v, order] = sort(mean(expl_var,1));
%             v = diff([0, v]);
%             obj.ref.plot_value_tree(v, ROI_list(ordered_beta(order)),'',labels{beh_idx},'',beh_idx,'','viridis','invalid_color');
%             if strcmpi(metric, 'explained_variance')
%                 caxis([0,max(v)]);
%             end
            
        if rendering && ~isempty(model)
            next_fig
            if strcmpi(metric, 'cum_explained_variance')
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'',[cmap, 'symmetrical'],'invalid_color');
                set(gca,'Color',[0.9,0.9,0.9])
            elseif strcmpi(metric, 'beta')
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'','parula','invalid_color');
                set(gca,'Color',[0.9,0.9,0.9])
            elseif strcmpi(metric, 'order')
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'','parula','invalid_color');
                set(gca,'Color',[0.9,0.9,0.9])
            elseif strcmpi(metric, 'explained_variance')
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'','parula','invalid_color');
                caxis([0,max(nanmean(all_scores{beh_idx},1))]);
            elseif strcmpi(metric, 'partialcorr')  
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'','hot','invalid_color');
                % Plotting the results
                figure('Color', 'w'); % Create a new figure with white background
                bar(median(all_scores{beh_idx},1)); % Create a bar chart
                xlabel('Predictor'); % Label x-axis
                ylabel('Squared Semi-partial Correlation'); % Label y-axis
                title('Squared Semi-partial Correlations for each Predictor'); % Title for the graph
                
            else
                obj.ref.plot_value_tree(median(all_scores{beh_idx},1), ROI_list,'',labels{beh_idx},'',beh_idx,'',cmap,'invalid_color');
            end 
        end
    end
    
    
    %% Control block
    % display bet coefs for all iteration for a specific behaviour
    %     close all
    %     beh_idx = find(strcmp(training_output{1}.beh_type, 'BodyCam_L_whisker'));
    %     figure();plot(abs(all_scores{beh_idx})');
    %     title('Beta coefs across 20 iterations')
    %     set(gcf, 'Color', 'w');
    %     xlabel('ROIs #')
    %     ylabel('beta (abs)')
    
%     close all
%     ref = nanmean(cat(1, all_scores{:}));
%     for beh_idx = not_shuffled_idx
%         new = nanmean(all_scores{beh_idx}, 1);
%         new = new / (new / ref);
%         new = new - ref;
%         next_fig
%         obj.ref.plot_value_tree(new, ROI_list,'',labels{beh_idx},'',beh_idx,'','redblue','invalid_color');
%         if strcmpi(metric, 'explained_variance')
%             caxis([-5, 5]);
%         end
%     end
    
    
    
    
end

function semipartial_r2 = semi_partial_correlation_coefficients(true_behaviours, predictors, ml_params)
    use_kfold           = false;
    optimize_hyper  = true;
    
    SSR_full = get_prediction(true_behaviours, predictors, use_kfold, optimize_hyper, ml_params);
       
    %     optimalIndices = sequentialfs(fun, predictors', true_behaviours', 'direction', 'backward');
    %     disp(find(optimalIndices));


    % initialize an array to hold semi-partial correlation coefficients
    semipartial_corr = zeros(1, size(predictors, 1));

    % loop over predictors
    for i = 1:size(predictors, 1)
        round(i/size(predictors, 1)*100, 1)
        % reduced model: fit model without predictor i
        predictors_reduced = predictors;
        predictors_reduced(i, :) = [];
        
        SSR_reduced = get_prediction(true_behaviours, predictors_reduced, use_kfold, optimize_hyper, ml_params);

        % semi-partial correlation coefficient
        semipartial_corr(i) = (SSR_reduced - SSR_full) / SSR_reduced;
        
        if abs(semipartial_corr(i)) > 1 % that should happen because the SSR should be lower than in the full model, however, if there is an overfitting issue, it's possible
            1
        end
    end

    % Squaring the semi-partial correlation coefficients to get squared semi-partial correlation coefficients
    semipartial_r2 = semipartial_corr.^2;

    % Plotting the results
    figure('Color', 'w'); % Create a new figure with white background
    bar(semipartial_r2); % Create a bar chart
    xlabel('Predictor'); % Label x-axis
    ylabel('Squared Semi-partial Correlation'); % Label y-axis
    title('Squared Semi-partial Correlations for each Predictor'); % Title for the graph
end

function SS_res = get_prediction(observation, predictors, use_kfold, optimize_hyper, ml_params)
    base_varargin = {predictors'     , ...
                     observation , ...                              
                     'PostFitBias'   , ml_params.postfit_bias,...
                     'PassLimit'     , 10,...
                     'Learner'       , 'leastsquares'}; % default is no regularization, no standardization (e.g. data is not mean centered)
                        
    if optimize_hyper            
        Lambda        = linear_hyperparameters_optimization(base_varargin, @fitrlinear, '', '', ml_params);   
        base_varargin = [base_varargin, {'Lambda', Lambda}];
    end
    
    if use_kfold > 1
        base_varargin = [base_varargin, {'KFold', 5}]; 
    end

    % full model
    mdl_full = fitrlinear(base_varargin{:});
    % SSR for the full model
    if ~use_kfold        
    	yhat_full = predict(mdl_full, predictors');
    else
        yhat_full = kfoldPredict(mdl_full)';
    end
    
    SS_res = sum((observation - yhat_full).^2);
end
