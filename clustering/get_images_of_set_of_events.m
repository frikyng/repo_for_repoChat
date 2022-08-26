for beh_idx = 2%1:3
    cl_idx = cluster_idx{beh_idx}  ;  

    rescaled_traces     = obj.rescaled_traces(:, signal_indices);
    cluster_traces      = [];
    for gp = sort(unique(cl_idx))'
        if gp ~= -1
            cluster_traces(gp,:) = nanmean(rescaled_traces(:,cl_idx == gp),2);
        end
    end

    event_index{beh_idx} = {};
    for current_cl_idx = 1%unique(cl_idx(cl_idx ~= -1))'
        
        beh_idx = 2
        current_cl_idx = 5;
        a = nanmedian(source_signal(:, find(cluster_idx{beh_idx} == current_cl_idx)),2);
        a = a-mode(a);
        tp = find(obj.get_tp_for_condition(conditions{2}));
        tp = tp(1);
        figure();findpeaks(a, 'MinPeakHeight',0.5);
        [~, peak_t_new] = findpeaks(a, 'MinPeakHeight',0.5);
        peak_t = peak_t_new(peak_t_new > tp)
        
        %% build difference trace
        a = cluster_traces(current_cl_idx,:)' - nanmean(obj.binned_data.median_traces,2);
        a = movmean(a,10);


        %% Find large upmodulated events
        figure();findpeaks(a, 'MinPeakHeight',0.5);
        [~, peak_t_new] = findpeaks(a, 'MinPeakHeight',0.5);

        %% List all events
        peak_t = vertcat(obj.event.peak_time{:});
        %% Replace follwoing loop by this
        %[event_infos, current_start_tp] = get_trial_nb_from_tp(expe, event_times, assign_to_bAP, merged_events)
        event_index{beh_idx}{current_cl_idx} = [];
        current = 1;
        for record = 1:numel(obj.timescale.tp)
            start = current;
            stop = current + obj.timescale.tp(record) - 1;
            incl = find(peak_t >= start & peak_t <= stop);
            if ~isempty(incl)                    
                for ev = incl'
                    ev_ori = find(cellfun(@(x) any(ismember(x,peak_t(ev))), obj.event.t_win));    
                    if ~isempty(ev_ori) % QQ make sure why this can be empty
                        pt_in_trial = peak_t(ev) - start + 1;
                        event_index{beh_idx}{current_cl_idx}(end+1,:) = [record, ev, ev_ori, pt_in_trial];
                    end
                end
            end   
            current = current + obj.timescale.tp(record);
        end
    end
end

cond = 2;
im{cond} = {};
for rec = unique(event_index{beh_idx}{current_cl_idx}(:,1))'
    subset = find(event_index{beh_idx}{current_cl_idx}(:,1) == rec)';
    im{cond}  = [im{cond} , load_bap(obj, rec, event_index{beh_idx}{current_cl_idx}(subset,4)')];
    close all    
end
figure();imagesc(nanmean(cat(3,im{cond}{:}),3))
im2 = im;
for el = 1:numel(im2{cond} )
    im2{cond} {el} = im2{cond} {el} / prctile(im2{cond} {el}(:), 99.9);
end
figure();imagesc(nanmean(cat(3,im2{cond}{:}),3))
show_stack_GUI(cat(3,im2{cond}{:}))


% %% Find corresponding baps
% for t = peak_t_new'
%     all         = vertcat(obj.event.peak_time{:});
%     [~, loc]    = nanmin(abs(all-t(1))');
%     record      = event_index(loc,1);
%     event       = event_index(loc,2);
%     idx         = find(contains({app.Tree.Children(record).Children.Text}, ['Event ',num2str(event),' t = ']));
%     app.Tree.CheckedNodes =  [app.Tree.CheckedNodes,;app.Tree.Children(record).Children(idx)];
% end
        