function show_events(all_traces, events, pre_peak_delay)
    bsl = prctile(all_traces, 1);
    figure(125);clf();hold on;
    ax1 = subplot(1,2,1);title('events')
    ax2 = subplot(1,2,2);title('scaled events')
    for peak = 1:size(events, 3)
        subplot(ax1);cla();
        subplot(ax2);cla();ylim([-0.1,1.1])
        for gp = 1:size(all_traces, 2)
            temp = events(:, gp, peak);
            subplot(ax1); hold on;plot(ax1, temp);
            temp = temp - bsl(gp);
            subplot(ax2); hold on;plot(ax2, temp ./ nanmax(temp(pre_peak_delay-3:pre_peak_delay+3)));
        end
        drawnow
    end
end