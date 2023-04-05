%% example bar_chart(result, multirun{1}.beh_type)
%% example bar_chart(result, 'beh_type')
%% example bar_chart(result, [strcat('group ', strsplit(num2str(1:numel(groups)),' '))])



function [meanvalue, sem_values, fig_handle, stats_results, values] = bar_chart(result, behaviour_filters, result_fieldname, additional_handle, condition_labels, rendering, do_stats, varargin)
    if nargin < 2 || isempty(behaviour_filters)
        behaviour_filters = '';
    end
    if nargin < 3 || isempty(result_fieldname)
        result_fieldname = 'score';
    end
    if nargin < 4 || isempty(additional_handle)
        if strcmp(result_fieldname, 'score')
            additional_handle = @(x) x(end);
        else
            additional_handle = [];
        end
    end
    if nargin < 5 || isempty(condition_labels)
        condition_labels = '';
    end
    if nargin < 6 || isempty(rendering)
        rendering = true;
    end
    if nargin < 7 || isempty(do_stats)
        do_stats = true;
    end
    
    N_behaviours    = 1;
    N_iter          = 1;
    N_conditions    = 1;
    
    if isstruct(result) 
        % One training
        N_behaviours    = numel(result.score);
        result          = {{result}};
    elseif iscell(result) && isstruct(result{1})    
        % Multiple training
        N_iter          = numel(result);
        N_behaviours    = numel(result{1}.score);
        result          = {result};
    elseif iscell(result{1}) && isstruct(result{1}{1})
        % Multiple training, multiple conditions
        N_conditions    = numel(result);
        N_iter          = numel(result{1});
        N_behaviours    = max(cellfun(@(x) numel(x{1}.score), result));%numel(result{1}{1}.score);
    end
    
    %% In case some condition have a different number of model iteration (eg : shuffle), add the missing ones
    N_cond_per_file = cellfun(@(x) numel(x{1}.score), result);
    [N_max, max_loc]        = max(N_cond_per_file);
    if numel(unique(N_cond_per_file)) > 1
        for cond_idx = 1:N_conditions
            for iter_idx = 1:N_iter
                if numel(result{cond_idx}{iter_idx}.score) == max(N_cond_per_file)
                    % pass  
                else
                    result{cond_idx}{iter_idx}.(result_fieldname) = cellfun(@(x) x(1), result{cond_idx}{iter_idx}.(result_fieldname), 'UniformOutput', false);
                    result{cond_idx}{iter_idx}.(result_fieldname)(1:2:N_max) = result{cond_idx}{iter_idx}.(result_fieldname);
                    empty = NaN(size(result{cond_idx}{iter_idx}.(result_fieldname){1}));
                    result{cond_idx}{iter_idx}.(result_fieldname)(2:2:N_max) = repmat({empty},1,N_max/2);
                    result{cond_idx}{iter_idx}.beh_type = result{max_loc}{iter_idx}.beh_type;   
                end
            end        
        end
    end
    
    %% Now, extract the condition
    values = num2cell(NaN(N_behaviours, N_iter, N_conditions));
    for cond_idx = 1:N_conditions
        for iter_idx = 1:N_iter
            for beh_idx = 1:N_behaviours  
                values{beh_idx, iter_idx, cond_idx} = result{cond_idx}{iter_idx}.(result_fieldname){beh_idx};   
            end 
        end        
    end
    
    labels  = result{1}{1}.('beh_type');
    labels  = strrep(labels, '_', '\_');        
    labels  = reordercats(categorical(labels),labels);
    prefilter_labels  = categories(labels);
    
    if ischar(behaviour_filters) && isempty(behaviour_filters)
        to_keep     = true(size(prefilter_labels));
    elseif (ischar(behaviour_filters) && ~isempty(behaviour_filters)) || iscell(behaviour_filters)         
        to_keep     = cellfun(@(x) any(contains(strrep(x, '\_', '_'), behaviour_filters)), prefilter_labels);
    end

    INVALID     = ~to_keep;% | all(cellfun(@isempty, values),2);
    INVALID     = INVALID(:,1,1); % may have to tweak this for more complex scenario
    values(INVALID, : , :) = [];
    for condition = 1:numel(result)
        for iter    = 1:numel(result{condition})
            for f       = fieldnames(result{condition}{iter})'
                if numel(result{condition}{iter}.(f{1})) == numel(labels)
                    result{condition}{iter}.(f{1})(:,INVALID) = [];
                end
            end        
        end
    end
    
    labels  = result{1}{1}.('beh_type');
    labels  = strrep(labels, '_', '\_');        
    labels  = reordercats(categorical(labels),labels);
    
    is_shuffle = 0;
    N_behaviours_unique = numel(labels);
    if any(contains(categories(labels), 'shuffle'))
        is_shuffle = 1;
        N_behaviours_unique = N_behaviours_unique / 2;
    end
    


    if ~isempty(additional_handle) % qq if non scalar output, need to be adjusted       
        values = cellfun(@(y) additional_handle(y), values);
    end
    
    %value = cell2mat(cellfun(@(x) cell2mat(cellfun(@(y) y.(result_fieldname)), x)'), result, 'UniformOutput', false));

    %% Get mean value if multiple iterations of more than one variable
    meanvalue   = permute(nanmean(values,2),[1,3,2]);
    sem_values  = permute(nanstd(values,[], 2)./sqrt(size(values, 2)),[1,3,2]);
    initial_names = categories(labels);
    
    %% Set default option in case no varargin is passed
    max_scores      = [];
    score_metrics   = 'performance';
    show_dots       = true;
    split_shuffle   = true;   
    if nargin >= 8 && isstruct(varargin{1})
        ml_parameters = varargin{1};
        if isfield(ml_parameters, 'max_score') && numel(ml_parameters.max_score) == numel(prefilter_labels)
            max_scores = ml_parameters.max_score;
            if numel(prefilter_labels) ~= numel(initial_names) % if you had a filter
                max_scores(INVALID) = [];
            end
            max_scores = max_scores(1:2:end);
        elseif isfield(ml_parameters, 'max_score') && numel(ml_parameters.max_score) == numel(meanvalue)
            max_scores = ml_parameters.max_score(1:2:end);                
        end

        if isfield(ml_parameters, 'score_metrics')
            score_metrics = ml_parameters.score_metrics;
            if ishandle(score_metrics)
                score_metrics = handle2str(score_metrics);
            end
            if contains(score_metrics, 'pearson')
                score_metrics = 'Predictive Score (r)';
            elseif strcmpi(score_metrics, 'mse')
                score_metrics = '1/mse';
            elseif strcmpi(score_metrics, 'rmse')
                score_metrics = '1/rmse';
            else
                score_metrics = 'model performance';
            end
        end

        if isfield(ml_parameters, 'show_dots')
            show_dots = ml_parameters.show_dots;
        end
        if isfield(ml_parameters, 'split_shuffle')
            split_shuffle = ml_parameters.split_shuffle;
        end

        if N_conditions > numel(labels) % multi conditions
            max_scores = repmat(max_scores, 1, numel(meanvalue));
        end
    end

    %% Generate figure
    fig_handle      = [];
    shuffle_values  = [];
    shuffle_semvalues     = [];
    shuffle_meanvalue    = [];
    if rendering
        %% Fix for labels
        labels          = fix_labels(labels);
        fig_handle      = figure();clf();
        fig_handle.Position(4) = fig_handle.Position(4)*1.3;
        if is_shuffle
            idx             = reshape(reshape(1:size(meanvalue,1),2,size(meanvalue,1 )/2)',[],2);
            labels          = removecats(labels,initial_names(idx(:,2)));
            labels          = labels(~isundefined(labels));
            if size(idx, 1) > 1 % otherwise reshaping error when only one behaviour
                meanvalue       = meanvalue(idx)'; 
            end
            if split_shuffle
                try
                    hb              = bar(labels, meanvalue', 'EdgeColor', 'none', 'facecolor', 'flat', 'BarWidth',1.5);hold on;
                catch
                    hb              = bar(meanvalue', 'EdgeColor', 'none', 'facecolor', 'flat', 'BarWidth',1.5);hold on;
                end
                hb(2).FaceColor = [0.8,0.8,0.8];            
                hb(1).FaceColor = [0 ,0.4470,0.7410];
                hb(1).FaceColor = 'flat';
            else
                hb(1:N_conditions)  = bar(labels,meanvalue(1,:)', 'EdgeColor', 'none', 'facecolor', 'flat', 'FaceColor', [0 ,0.4470,0.7410]);hold on;
                hb(N_conditions+1:N_conditions*2)           = bar(labels,meanvalue(2,:)', 'EdgeColor', 'none', 'facecolor', 'flat', 'FaceColor', [0.8,0.8,0.8]);hold on;
                %hb(1:N_conditions).FaceColor = 'flat';
            end
            shuffle_values      = values(idx(:,2),:,:);
            values              = values(idx(:,1),:,:);
            shuffle_semvalues   = sem_values(idx(:,2),:,:);
            sem_values          = sem_values(idx(:,1),:,:);
            shuffle_meanvalue   = meanvalue(2,:,:)';
            meanvalue           = meanvalue(1,:,:)';
        else
            hb                  = bar(labels, meanvalue, 'EdgeColor', 'none','FaceColor',lines(1));hold on;
            colororder(viridis(N_conditions))
        end        

        %% Adjust max bar width depending on shuffle condition
        if ~isempty(max_scores)
            if is_shuffle
                w = hb(1).BarWidth/4;
            else
                w = hb(1).BarWidth/2;
            end            
            %             for el = 1:numel(max_scores)
            %                 plot([el-w,el+w] ,[max_scores(el), max_scores(el)],'r:', 'LineWidth',2)
            %             end
        end

        %% Add error bars. This needs a trick for categorical data
        if ~is_shuffle
            n_sub_bars = 1;
        else
            n_sub_bars = 1:2;
        end
        if split_shuffle
            condition_list = N_conditions;
        else
            condition_list = 1;
        end        
        for condition = 1:condition_list
            for behaviour = 1:N_behaviours_unique
                x = [];
                x_shuff = [];
                if N_conditions == 1
                    bar_idx = behaviour;
                else
                    bar_idx = condition;
                end               
                
                for sub_bar = n_sub_bars % 1 if not shuffle, 2 otherwise
                    if N_conditions > 1 && ~split_shuffle
                        sub_bar_idx = (N_conditions*(sub_bar-1)+1) : N_conditions*(sub_bar);                             
                    else
                        sub_bar_idx = sub_bar;
                    end 
                    
                    if sub_bar == 1 && split_shuffle
                        x = [x ; hb(sub_bar).XEndPoints(bar_idx)];
                    elseif sub_bar == 2 && is_shuffle && split_shuffle
                        x_shuff = [x_shuff ; hb(sub_bar).XEndPoints(bar_idx)]; 
                    elseif sub_bar == 1 && is_shuffle && ~split_shuffle
                        x = [x ; arrayfun(@(x) x.XEndPoints(bar_idx), hb(sub_bar_idx))];
                    elseif sub_bar == 2 && is_shuffle && ~split_shuffle
                        x_shuff = [x_shuff ; arrayfun(@(x) x.XEndPoints(bar_idx), hb(sub_bar_idx))]; 
                    end

                    if show_dots    
                        if sub_bar == 1 && numel(sub_bar_idx) == 1
                            scatter(hb(sub_bar_idx).XEndPoints(bar_idx), values(behaviour,:,condition), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.15); hold on;
                        elseif sub_bar == 1 && numel(sub_bar_idx) > 1                            
                            scatter(arrayfun(@(x) x.XEndPoints(bar_idx), hb(sub_bar_idx)), squeeze(values(behaviour,:,:)), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.15); hold on;
                        elseif N_conditions == 1 && sub_bar == 2 && split_shuffle
                            scatter(hb(sub_bar).XEndPoints(bar_idx), shuffle_values(behaviour,:,condition), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.15); hold on;
                        elseif N_conditions > 1 && split_shuffle && sub_bar == 2
                            scatter(hb(sub_bar).XEndPoints(bar_idx), shuffle_values(behaviour,:,condition), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.15); hold on;
                        end
                    end
                end

                if is_shuffle
                    if ~split_shuffle && N_conditions > 1
                        err_bar_idx = 1:N_conditions;
                    elseif N_conditions == 1
                        err_bar_idx = behaviour;  
                    else
                        err_bar_idx = condition;
                    end
                    errorbar(x',meanvalue(err_bar_idx),sem_values(err_bar_idx),'k','linestyle','none','LineWidth',2,'CapSize',0);hold on;
                    errorbar(x_shuff',shuffle_meanvalue(err_bar_idx),shuffle_semvalues(err_bar_idx),'linestyle','none','LineWidth',2,'CapSize',0,'Color',[0.6,0.6,0.6]);hold on;
                else
                    errorbar(x',meanvalue(bar_idx),sem_values(bar_idx),'k','linestyle','none','LineWidth',2,'CapSize',0);hold on;
                end
            end
        end
        
%         if is_shuffle && ~split_shuffle && N_conditions > 1 
%             % https://fr.mathworks.com/matlabcentral/answers/462532-multi-bar-labeling-plot
%             xCnt = cell2mat(get(hb,'XOffset')); % XOffset is undocumented!
%             xCnt = xCnt(1:end/2);        
%             xLab = categorical({'A','B','C','D','E'});
%             set(gca, 'XTick', xCnt, 'XTickLabel', xLab)
%         end

        %% Adjust plot limits
        ylim([-10,100]);hold on;

        if ~isempty(condition_labels) && numel(condition_labels) == size(meanvalue, 1)
            %legend(condition_labels);
        elseif ~isempty(condition_labels) && numel(condition_labels) ~= size(meanvalue, 1)
            error(['labels must be a string array of ', num2str(size(meanvalue, 1)), ' elements'])
        end
        
        ylabel(score_metrics);
        set(fig_handle, 'Color', 'w'); hold on;
        set(gca,'box','off')
    elseif ~rendering && is_shuffle
        idx             = reshape(reshape(1:size(meanvalue,1),2,size(meanvalue,1 )/2)',[],2);
        if size(idx, 1) > 1 % otherwise reshaping error when only one behaviour
            meanvalue       = meanvalue(idx)'; 
        end
        shuffle_values      = values(idx(:,2),:,:);
        values              = values(idx(:,1),:,:);
        shuffle_semvalues   = sem_values(idx(:,2),:,:);
        sem_values          = sem_values(idx(:,1),:,:);
        shuffle_meanvalue   = meanvalue(2,:,:)';
        meanvalue           = meanvalue(1,:,:)';
    end
    
    stats_results = {};    
    if do_stats
        if N_conditions > 1
           values           = permute(values,[3,2,1]);
           shuffle_values   = permute(shuffle_values,[3,2,1]);
           labels = categorical(condition_labels);
           xticklabels(gca, categorical(condition_labels));
        end
        stats_results.equalvariances = vartestn(values','TestType','LeveneAbsolute','Display','off') > 0.05;
        if size(values, 1) > 2
            [stats_results.p,stats_results.table,stats_results.stats] = kruskalwallis(values',[],'off');
        else
            [stats_results.p,~,stats_results.stats] = ranksum(values(1,:),shuffle_values(1,:));
        end
%         if rendering
%             stat_disp = 'on';
%         else
            stat_disp = 'off';
%         end
        
        ignore_list = [];
        if ~isempty(shuffle_values)
            for beh = 1:size(values,1)
                if ~all(isnan(shuffle_values(beh,:)))
                    stats_results.shuffle_stat(beh) = ranksum(values(beh,:),shuffle_values(beh,:));
                    if stats_results.shuffle_stat(beh) >= 0.05
                        hb(1).CData(beh,:) = [0.6,0.6,0.6];
                        ignore_list(end+1) = beh;
                    end
                end
            end            
        end
        
        if numel(ignore_list) > 1
            ignore_list = nchoosek(ignore_list,2);
        else
            ignore_list = [];
        end
        
        if stats_results.p < 0.05
            if size(values, 1) > 2
                disp('some conditions are significantly different from others')
                try
                [stats_results.p,stats_results.table,stats_results.stats]   = kruskalwallis(values', labels, stat_disp);
                catch
                    return
                end
                stats_results.multi_comp                                    = multcompare(stats_results.stats, 'Display', stat_disp);
            else
                disp('the two conditions are significantly different');
                name = categories(labels);
                [~,stats_results.table,stats_results.stats]                 = kruskalwallis([values', shuffle_values'], {name{1}, 'shuffle'}, stat_disp); % should be replaced by MW
                stats_results.multi_comp                                    = multcompare(stats_results.stats, 'Display', stat_disp);
            end
            
            if ~split_shuffle && N_conditions > 1
                for el = 1:N_conditions
                    col1 = stats_results.multi_comp(:,1);
                    col2 = stats_results.multi_comp(:,2);
                    stats_results.multi_comp(col1 == el, 1) = hb(el).XEndPoints;
                    stats_results.multi_comp(col2 == el, 2) = hb(el).XEndPoints;
                end
            elseif ~split_shuffle && N_conditions == 1
                 stats_results.multi_comp(:,[1,2]) = ceil(stats_results.multi_comp(:,[1,2])/2);
            end

            if rendering
                figure(fig_handle)
                
                
                %non_parametric_chart(values, input_labels, errors)
                
                ns          = stats_results.multi_comp(stats_results.multi_comp(:,6) >= 0.05, [1,2]);
                p0_01       = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.05 & stats_results.multi_comp(:,6) >= 0.01, [1,2,6]);
                p0_001      = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.01 & stats_results.multi_comp(:,6) >= 0.001, [1,2,6]);
                p_strong    = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.001, [1,2,6]);

                if ~isempty(ignore_list)
                    p0_01 = p0_01(~ismember(p0_01(:,[1,2]), ignore_list, 'rows'),:);
                    p0_001 = p0_001(~ismember(p0_001(:,[1,2]), ignore_list, 'rows'),:);
                    p_strong = p_strong(~ismember(p_strong(:,[1,2]), ignore_list, 'rows'),:);
                end
                
                [M,loc]     = max(meanvalue);
                M           = M + sem_values(loc);
                step        = (100 - M)/ (sum(stats_results.multi_comp(:,6) < 0.05)+1);
                for el = 1:size(p0_01, 1)
                    M           = M+step;
                    overbar(p0_01(el,1) ,p0_01(el,2), M,'*'); %  ['p=',num2str(round(p0_01(el,3),3))]
                end
                for el = 1:size(p0_001, 1)
                    M           = M+step;
                    overbar(p0_001(el,1) ,p0_001(el,2), M, '**'); %['p=',num2str(round(p0_001(el,3),3))]
                end
                for el = 1:size(p_strong, 1)
                    M           = M+step;
                    overbar(p_strong(el,1) ,p_strong(el,2), M, '***'); %  ['p=',num2str(round(p_strong(el,3),3))]
                end
            end
        else
            disp(['p = ',num2str(stats_results.p),'; Conditions do not significantly differ from each other'])
        end
    end
end

