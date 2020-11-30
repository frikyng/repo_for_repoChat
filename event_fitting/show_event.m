function show_event(all_traces, peak_times, event_nb)
    pre_peak_delay  = 20;
    post_peak_delay = 50;

    if peak_times(event_nb)+post_peak_delay > size(all_traces, 1)
        post_peak_delay = size(all_traces, 1) -  peak_times(event_nb);
    elseif peak_times(event_nb)-pre_peak_delay < 1
        pre_peak_delay = peak_times(event_nb) - 1;
    end
    event = all_traces((peak_times(event_nb)-pre_peak_delay):(peak_times(event_nb)+post_peak_delay),:);
    figure(123);cla();plot(event)
end

