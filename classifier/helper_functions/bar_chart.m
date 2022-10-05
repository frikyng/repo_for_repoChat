%% example bar_chart(result, multirun{1}.beh_type)
%% example bar_chart(result, 'beh_type')
%% example bar_chart(result, [strcat('group ', strsplit(num2str(1:numel(groups)),' '))])


function bar_chart(result, labels_or_label_fieldname, result_fieldname, additional_handle)
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

    if ischar(labels_or_label_fieldname)
    	labels = reordercats(categorical(result{1}{1}.(labels_or_label_fieldname)),result{1}{1}.(labels_or_label_fieldname));
    else
        labels % [strcat('group ', strsplit(num2str(1:numel(groups)),' '))])
    end
    
    if ~isempty(additional_handle) % qq if non scalar output, need to be adjusted       
        value = cellfun(@(y) additional_handle(y), value);
    end
    
    %value = cell2mat(cellfun(@(x) cell2mat(cellfun(@(y) y.(result_fieldname)), x)'), result, 'UniformOutput', false));

    %% Get mean value if multiple iterations of more than one variable
    meanvalue   = squeeze(nanmean(value,2));
    sem_values  = squeeze(nanstd(value,[], 2)./sqrt(size(value, 2)));

    %% Generate figure
    figure();hb = bar(labels, meanvalue);hold on; colororder(viridis(N_conditions))
    
    %% Add error bars. This needs a trick for catehorical data
    nbars = size(meanvalue, 2);
    x = [];
    for i = 1:nbars
        x = [x ; hb(i).XEndPoints];
        scatter(hb(i).XEndPoints, value(:,:,i), 'MarkerEdgeColor','none', 'MarkerFaceColor', 'k', 'MarkerFaceAlpha', 0.3); hold on;
    end
    errorbar(x',meanvalue,sem_values,'k','linestyle','none');hold on;
    

    %% Adjust plot limits
    %ylim([-10,50]);hold on;

end

