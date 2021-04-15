% example 
% animate_tree(current_expe, repmat(mean_trace, 1, 52), 1:1000)
% animate_tree(current_expe, current_expe.rescaled_traces, 1:1000)
% animate_tree(current_expe, vertcat(current_expe.extracted_traces{:}), 1:1000)

function animate_tree(current_expe, values, tp_range, pause_t)
    if nargin < 3 || isempty(tp_range)
        tp_range = 1:size(values,1);
    end
    if nargin < 4 || isempty(pause_t)
        pause_t = 0;
    end
    
    idx           = current_expe.ref.indices.ROIs_list;
    h = current_expe.ref.plot_value_tree(nanmean(values),'',current_expe.default_handle,'','',1); % 
    t = title('tp = 0');
    duplicated_value = [false, (diff(h.XData(1,:)) == 0 & diff(h.YData(1,:)) == 0 & diff(h.ZData(1,:)) == 0)];
    ok = ~isnan(h.XData(1,:)) %& ~duplicated_value ; % Nan are added for curved scan at the end of each branch

    caxis([0,1])


    for tp = tp_range
        %v = corr_results_sub(tp,idx);
        R = tp:nanmin(tp+0, size(values, 1));
        v = nanmean(values(R,idx),1);  

        h.CData(:, ok) = repmat(v, 2,1);
        t.String = ['tp = ',num2str(tp)];
        drawnow
        pause(pause_t)
    end
end