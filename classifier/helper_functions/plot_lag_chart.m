function [fig_handle, fig_handle_cross] = plot_lag_chart(results, lag_list, sr, behaviours, extracted_beh, lag_smoothing)
    if nargin < 6 || isempty(lag_smoothing)
        lag_smoothing = 0;
    end
    results = cell2mat(results);
    if any(lag_smoothing)
        results = smoothdata(results,2,'gaussian',lag_smoothing);
    end
    x = fliplr(lag_list/sr);
    fig_handle = figure(1);clf();
    tiledlayout('flow');
    for beh = 1:numel(behaviours)
        axes(beh) = nexttile    ;  
        patch([x,fliplr(x)], [results(beh, :), zeros(size(results(beh, :)))], 'k');
        title(strrep(behaviours{beh},'_','\_'));xlabel('Lag (s)');ylabel('Prediction score');
        xline(0,'r--');
    end
    linkaxes(axes, 'xy')
    
    fig_handle_cross = figure(2);clf();
    tiledlayout('flow');
    for beh = 1:size(extracted_beh.raw_behaviours, 1)
        axes(beh) = nexttile    ;
        title(strrep(extracted_beh.original_behaviour_list{beh},'_','\_'));hold on
        plot(linspace(min(x), max(x), max(lag_list)*2+1), xcorr(extracted_beh.raw_behaviours(beh,:),max(lag_list),'normalized')); hold on
        xlabel('Lag (s)');xline(0,'r--')
    end
end

