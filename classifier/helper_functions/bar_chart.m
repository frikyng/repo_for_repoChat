%% example bar_chart(result, multirun{1}.beh_type)
%% example bar_chart(result, 'beh_type')
%% example bar_chart(result, [strcat('group ', strsplit(num2str(1:numel(groups)),' '))])


function [meanvalue, sem_values, fig_handle, stats_results, value] = bar_chart(result, labels_or_label_fieldname, result_fieldname, additional_handle, condition_labels, rendering, do_stats)
    if nargin < 2 || isempty(labels_or_label_fieldname)
        labels_or_label_fieldname = 'beh_type';
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
        N_behaviours    = numel(result{1}{1}.score);
    end
    
    %% Now, extract the condition
    value = num2cell(NaN(N_behaviours, N_iter, N_conditions));
    for cond_idx = 1:N_conditions
        for iter_idx = 1:N_iter
            for beh_idx = 1:N_behaviours                
                value{beh_idx, iter_idx, cond_idx} = result{cond_idx}{iter_idx}.(result_fieldname){beh_idx};               
            end            
        end        
    end
    
    INVALID = all(cellfun(@isempty, value),2);
    %value(INVALID, : ) = [];
    %result{1}{1}.(labels_or_label_fieldname)(INVALID) = [];

    if ischar(labels_or_label_fieldname)
        labels = result{1}{1}.(labels_or_label_fieldname);
        labels = strrep(labels, '_', '\_');        
    	labels = reordercats(categorical(labels),labels);
    else
        labels % [strcat('group ', strsplit(num2str(1:numel(groups)),' '))])
    end
   
    
    if ~isempty(additional_handle) % qq if non scalar output, need to be adjusted       
        value = cellfun(@(y) additional_handle(y), value);
    end
    
    %value = cell2mat(cellfun(@(x) cell2mat(cellfun(@(y) y.(result_fieldname)), x)'), result, 'UniformOutput', false));

    %% Get mean value if multiple iterations of more than one variable
    meanvalue   = permute(nanmean(value,2),[1,3,2]);
    sem_values  = permute(nanstd(value,[], 2)./sqrt(size(value, 2)),[1,3,2]);


    %% Generate figure
    fig_handle = [];
    if rendering
        fig_handle = figure();
        hb = bar(labels, meanvalue);hold on;
        colororder(viridis(N_conditions))

        %% Add error bars. This needs a trick for catehorical data
        nbars = size(meanvalue, 2);
        x = [];
        for i = 1:nbars
            x = [x ; hb(i).XEndPoints];
            scatter(hb(i).XEndPoints, value(:,:,i), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.3); hold on;
        end
        errorbar(x',meanvalue,sem_values,'k','linestyle','none');hold on;


        %% Adjust plot limits
        ylim([-10,100]);hold on;

        if ~isempty(condition_labels) && numel(condition_labels) == size(meanvalue, 2)
            legend(condition_labels);
        elseif  ~isempty(condition_labels) && numel(condition_labels) ~= size(meanvalue, 2)
            error(['labels must be a string array of ', num2str(size(meanvalue, 2)), ' elements'])
        end
        
        set(fig_handle, 'Color', 'w'); hold on;
    end
    stats_results = {};
    if do_stats
        [stats_results.p,stats_results.table,stats_results.stats] = kruskalwallis(value',[],'off');
        if stats_results.p < 0.05
            disp('some conditions are significantly different from others')
            [stats_results.p,stats_results.table,stats_results.stats]   = kruskalwallis(value', labels);
            stats_results.multi_comp                                    = multcompare(stats_results.stats);
            ns = stats_results.multi_comp(stats_results.multi_comp(:,6) >= 0.05, [1,2]);
            p0_01 = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.05 & stats_results.multi_comp(:,6) >= 0.01, [1,2]);
            p0_001 = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.01 & stats_results.multi_comp(:,6) >= 0.001, [1,2]);
            p0_strong = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.001, [1,2]);
        else
            disp('Conditions do not significantly differ from each other')
        end
    end
end

