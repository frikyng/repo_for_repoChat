function [all_scores, valid_trace_idx]  = extract_ML_coefs(obj, training_output, rendering, metric, use_absolute)
    if nargin < 3 || isempty(rendering)
        rendering = true;
    end
    if nargin < 4 || isempty(metric)
        metric = 'beta';
    end
    if nargin < 5 || isempty(use_absolute)
        use_absolute = false;
    end
    
    ROIs            = training_output{1}.used_ROIs; % Note that in the case of randomized shuffle, this actually changes from one training to the next
    
    labels          = fix_labels(training_output{1}.beh_type);
    valid_trace_idx = training_output{1}.used_ROIs;
    
    behaviours = training_output{1}.beh_type;
    all_scores = {};
    for beh_idx = 1:2:numel(behaviours)
        all_scores{beh_idx}   = [];
        for iter = 1:numel(training_output)
            %valid_trace_idx = training_output{iter}.used_ROIs;
            model           = training_output{iter}.model{beh_idx};
            if isprop(model, 'Trained')
                for trained_idx = 1:numel(model.Trained)
                    mdl = model.Trained{trained_idx}; 
                    if strcmpi(metric, 'beta')
                        all_scores{beh_idx}(trained_idx, :) = abs(mdl.Beta);
                    elseif strcmpi(metric, 'bias')
                        all_scores{beh_idx}(trained_idx, :) = mdl.Bias;
                    end
                end
            else
                all_scores{beh_idx}(iter, :) = model.Beta;
            end
        end
        if use_absolute
            all_scores{beh_idx} = all_scores{beh_idx};
        end
        
        if rendering
            next_fig
            obj.ref.plot_value_tree(median(all_scores{beh_idx},1), find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)),'',labels{beh_idx},'',beh_idx,'','redblue', 'invalid_color');
        end
    end
end

