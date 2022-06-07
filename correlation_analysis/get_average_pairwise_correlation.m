


function [corr_results_sub, n_high_corr, subpeaks_t] = get_average_pairwise_correlation(current_expe, corr_results, data_sm, rendering, condition, corr_results_rand, references)
    if nargin < 4 || isempty(rendering)
        rendering = true;        
    end
    if nargin < 5 || isempty(condition)
        % condition == 1 --> all TP
        % condition == 2 --> bAPs
        % condition == 3 --> btween bAPs
        condition = 1;        
    end
    if nargin < 6 || isempty(corr_results_rand)
        corr_results_rand = [];        
    end
    if nargin < 7 || isempty(references)
        references = 1:size(data_sm,2);        
    end
        
    comb                        = nchoosek(1:size(data_sm,2),2);

   % use_global_events = cond > 1;
    
    if condition == 1 
        global_events_t = true(size(corr_results, 1),1);        
    elseif condition == 2 % then we extract global events
        global_events_t = ismember(1:size(corr_results, 1), [current_expe.event.t_win{current_expe.event.is_global}]);
    elseif condition == 3 % then we exclude global events
        global_events_t = ~ismember(1:size(corr_results, 1), [current_expe.event.t_win{current_expe.event.is_global}]);
    end
    corr_results(~global_events_t,:) = NaN;
    
    subpeaks_t = current_expe.event.t_corr;
    max_corr   = nanmax(current_expe.event.globality_index);

    %% Show mean correlation againt all other ROIs
    n_high_corr = [];
    full_range = 1:size(data_sm,2);
    for key = references
        if iscell(key)
            corr_results_sub    = NaN(size(data_sm));
            for el = full_range(~ismember(full_range, key{1}))
                corr_results_sub(:,el) = nanmean(corr_results(:, (comb(:,2) == el | comb(:,1) == el) & (any(comb(:,2) == key{1},2) | any(comb(:,1) == key{1},2))),2);
            end
            key = 1;
        else
            corr_results_sub        = corr_results(:, comb(:,2) == key | comb(:,1) == key);
            corr_results_sub    = [corr_results_sub(:,1:(key-1)), NaN(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];
        end

        mean_corr           = nanmean(corr_results_sub,1);
        if condition < 3
            n_high_corr(key) = sum(mean_corr > max_corr/2);
        else % check in how many subpeaks this region is involved
            peak_thr        = nanmean(prctile(current_expe.event.fitting.pre_correction_peaks, 10))*1.5; % based on large event detect, set a small value
            min_local_size  = 0;
            local_event     = data_sm(subpeaks_t, key) > peak_thr;
            N_coactive      = sum(data_sm(subpeaks_t(local_event), :) > peak_thr & corr_results_sub(subpeaks_t(local_event), :) > 0.8,1); % at least 0.5 corr amd big enough
            n_high_corr(key) = sum(N_coactive > min_local_size);
            
            t_corr = repmat(subpeaks_t(local_event),1,sum(N_coactive > min_local_size));
            v_corr = data_sm(subpeaks_t(local_event),N_coactive >min_local_size);
            
            if n_high_corr(key) > 1
%                 figure(445);cla();plot(data_sm(:,N_coactive <= min_local_size),'Color',[0.8,0.8,0.8]); hold on;
%                 plot(data_sm(:,N_coactive > min_local_size)); hold on;scatter(t_corr(:),v_corr(:),'kv','filled');
%                 1
            end
            %figure(445);cla();plot(nanmean(real(corr_results_sub),2)); drawnow
           % n_high_corr(key) = sum(sum(corr_results_sub(subpeaks_t, :) > 0.8))/ numel(subpeaks_t); % across all local events, with how many other ROIs it correlates.
            %n_high_corr(key) = mean(sum(corr_results_sub(subpeaks_t, :) > 0.5)); % during local events, how often it belongs to it
        end
        if rendering
            %% TO DISPLAY THE MEAN CORRELATION
                % current_expe.ref.plot_value_tree(mean_corr);
            %% FOR A SPECIFIC TIMEPOINT
                % tp = 711;
                % current_expe.ref.plot_value_tree(corr_results_sub(tp,:))
            %% FOR AN ANIMATION WITH ALL TP
            current_expe.ref.animate_tree(corr_results_sub,1:size(corr_results_sub, 1),'',[0,1])
            title(num2str(key))
            pause(1)
        end
    end
    
    global_events_t = find(global_events_t);
   % subpeaks_t = global_events_t(subpeaks_t);
    
end