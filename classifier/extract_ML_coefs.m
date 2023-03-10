function all_scores = extract_ML_coefs(obj, training_output, rendering, metric)
    if nargin < 3 || isempty(rendering)
        rendering = true
    end
    if nargin < 4 || isempty(metric)
        metric = 'beta'
    end
    
    ROIs = training_output{1}.used_ROIs; % Note that in the case of randomized shuffle, this actually changes from one traiing to the next
    
    behaviours = training_output{1}.beh_type;
    all_scores = {};
    for beh_idx = 1:numel(behaviours)
        all_scores{beh_idx}   = [];
        for iter = 1:numel(training_output)
            model = training_output{iter}.model{beh_idx};
            for trained_idx = 1:numel(model.Trained)
                mdl = model.Trained{trained_idx}; 
                if strcmpi(metric, 'beta')
                    all_scores{beh_idx}(trained_idx, :) = abs(mdl.Beta);
                elseif strcmpi(metric, 'bias')
                    all_scores{beh_idx}(trained_idx, :) = mdl.Bias;
                end
            end
        end
        if rendering
            next_fig
            obj.ref.plot_value_tree(mean(all_scores{beh_idx},1), find(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)),'',strrep(training_output{1}.beh_type{beh_idx},'_',' '),'',beh_idx,'','hot', 'invalid_color');
        end
    end
end

