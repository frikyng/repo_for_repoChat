function [event_infos, tp_list] = get_trial_nb_from_tp(expe, event_times, assign_to_bAP, merged_events)
    if nargin < 3 || isempty(assign_to_bAP)
        assign_to_bAP = true;
    end
    if nargin < 4 || isempty(merged_events)
        merged_events = true;
    end

    if merged_events
        var = expe.event.t_win;
    else
        var = expe.event.t_win_no_overlap;
    end

    event_infos      = [];
    tp_list          = {};
    current_start_tp = 1;
    for record = 1:numel(expe.timescale.tp)
        start_tp    = current_start_tp;
        stop_tp     = current_start_tp + expe.timescale.tp(record) - 1;
        included    = find(event_times >= start_tp & event_times <= stop_tp);

        if ~isempty(included) % if any event                   
            for ev = included
                if assign_to_bAP
                    containing_bAP  = find(cellfun(@(x) any(ismember(x,event_times(ev))), var));
                    if ~isempty(containing_bAP)
                        containing_bAP  = containing_bAP(end); % if more than one, we take the last one as it is more likely to include a response
                        pt_in_trial     = var{containing_bAP(end)} - start_tp + 1;
                    end
                else
                    containing_bAP  = NaN;
                    pt_in_trial     = event_times(ev) - start_tp + 1;
                end
                
                if ~isempty(containing_bAP) % empty when in assign_to_bAP mode and there is no corresponding bAP
                    event_infos(end+1,:)= [record, ev, containing_bAP];
                    tp_list{end+1}      = pt_in_trial;
                end
            end
        end   
        current_start_tp = current_start_tp + expe.timescale.tp(record);
    end
end