function fit_data = detect_and_fit_events(all_data, global_timescale, rendering, group_labels)
    if nargin < 3 || isempty(rendering)
        rendering = true;
    end
    fit_data        = {};
    fit_data.events = {};
    refit           = false;
    
    %% Find peaks using wavelet denoising on median trace
    [pk_val, peak_loc, peak_times, widths, denoised_data] = detect_events(nanmedian(all_data, 2), global_timescale');
    
    peak_times = peak_times{1}';
    peak_loc   = peak_loc{1}';
    %[~, peak_loc] = min(abs(peak_times - global_timescale'));
    
    %% Use peak width to help with fitting
    sampling_rate   = median(diff(global_timescale));
    peak_av_width   = nanmedian([widths{:}]); % that's the average peak width at half height, in points (median is in case there are a few trials with avery different sr
    % Based on it, we define a reasonable fitting window
    pre_peak_delay  = round(peak_av_width * 3);
    post_peak_delay = round(peak_av_width * 7);

    %% We'll do the fitting on denoised data
    all_data        = wavelet_denoise(all_data);
% 
%     %% QQ likely to mess up the rest of the analysis
%     if isempty(peak_times)
%         return
%     end
    
    %% Extract all events
    if isempty(peak_loc)
        peak_loc = NaN; % to avoid default detection
    end
    events = extract_events(all_data, peak_loc, pre_peak_delay, post_peak_delay);
    
    %% Get global fit settings to have a rough estimate of the fitting window required. Use isolated events (or end-of-burst)
    to_use = find(diff(peak_loc) > (peak_av_width * 5)); % keep events spaced by at least 5 average peak width from the next.
    figure(1019);cla();hold on;
    fit_data.global_fit = fit_decay(squeeze(nanmean(nanmedian(events(:,:,to_use), 2),3)), pre_peak_delay, [], rendering, [], peak_av_width);
    legend({'mean event','fit result'});
    if rendering
        figure(123);cla();
    end
    arrangefigures([1,2]); 
%     if rendering
%         show_events(all_data, events, pre_peak_delay)
%     end

    %% Now fit all events with fairly free parameters
    spmd
        warning('off')
    end
    [fits, all_peak_values, fit_data] = fit_all_events(events, all_data, peak_loc, pre_peak_delay, fit_data, rendering, [], peak_av_width);
    
    %% Make some figures about fit propeties
    fits_free = 1/permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
    figure(1014);clf();plot(fits_free(:,:,3)) ; hold on; plot(nanmedian(fits_free(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');xticklabels(group_labels);xtickangle(45);set(gcf,'Color','w');
    figure(1015);clf();plot(fits_free(:,:,3)'); hold on; plot(nanmedian(fits_free(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');
    fit_data.pre_refit_events = fit_data.events;
    arrangefigures([1,2]); 
    

    %% Re-Fit events using constraint for their group
    if refit
        %% Now refit but keep tau values within reason (+/- 1 mad)
        test = permute(cat(3, fit_data.events{:}),[1,3,2]);
        tau1 = nanmedian(test(:,:,3),2);
        tau2 = nanmedian(test(:,:,5),2);
        A1   = nanmedian(test(:,:,2),2);
        A2   = nanmedian(test(:,:,4),2);
        guesses = [zeros(size(tau1)), A1, tau1, A2, tau2];

        %show_events(all_data, events)

        [fits, all_peak_values, fit_data] = fit_all_events(events, all_data, peak_loc, pre_peak_delay, fit_data, ~rendering, guesses, peak_av_width);
    end
    
    %% Finally, correct the events amplitudes and store results
    fit_data.pre_correction_peaks = all_peak_values;
    all_peak_values_pre           = all_peak_values;
    [all_peak_values, correction] = correct_traces(all_data, events, fits, all_peak_values, round(peak_loc/sampling_rate));
    fit_data.post_correction_peaks= all_peak_values;
    fit_data.peak_times           = peak_times;
    fit_data.peak_pos             = peak_loc;
    R                             = max(range(all_data));
    for gp = 1:size(all_peak_values, 2)
        figure(10020 + gp);cla();hold on;title(['Traces and fits group ',num2str(gp),' (',group_labels{gp},')']);set(gcf,'Color','w');
        plot(global_timescale, all_data(:, gp),'k');hold on;ylim([-R/20,R + R/20])
        scatter(peak_times, nanmedian(all_peak_values_pre(:,gp), 2), 'bv', 'filled');hold on;
        scatter(peak_times, nanmedian(all_peak_values(:,gp), 2), 'kv', 'filled');hold on;
        for el = 1:size(fits, 1)
            tp = fits{el,gp}(:,1);
            tp = tp > 0 & tp <= size(global_timescale, 2); 
            plot(global_timescale(fits{el,gp}(tp,1)),fits{el,gp}(tp,2),'r--'); hold on;
        end
        drawnow
    end

    figure(1020);cla();hold on;plot(global_timescale, nanmedian(all_data, 2),'r','LineWidth', 1.5); hold on;plot(repmat(peak_times, size(all_peak_values, 2), 1)', double(all_peak_values));set(gcf,'Color','w');xlabel('Time (s)');ylabel('Amplitude');title('Event Amplitude trajectory');
    %figure(1002);hold on;scatter(peak_times, nanmean(all_peak_values, 2), 'bo', 'filled');hold on;plot(global_timescale, nanmean(correction, 2),'r--');set(gcf,'Color','w');

    %% Make some figures about fit propeties
    fits = permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
    figure(1014);clf();plot(fits(:,:,3)) ; hold on; plot(nanmedian(fits(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');set(gcf,'Color','w');ylabel('tau 1')%xticklabels(bin_legend);xtickangle(45)
    
    
    
%     
% figure(1014);
% ax = gca
% Xs = [];
% Ys = [];
% for el = 1:numel(ax.Children)
%     Xs = [Xs, ax.Children(el).XData];
%     Ys = [Ys, ax.Children(el).YData];
% end
%  
% Ys = sr./Ys;
% 
% Ym = [];
% Ysd= [];
% YN = [];
% Xm = unique(Xs);
% for bin = Xm
%     Ym(bin) = nanmean(Ys(Xs == bin));  
%     Ysd(bin) = nanstd(Ys(Xs == bin));  
%     YN(bin) = numel(Ys(Xs == bin));  
% end
% 
% figure();plot(Xm,Ym, 'ko-'); hold on
% errorbar(Xm, Ym, Ysd ./ sqrt(YN), 'ko-') ;
% title(['mean tau per group']);set(gcf,'Color','w');
% xlabel('bin')
% ylabel('tau decay')
% xticklabels(bin_legend);xtickangle(45);
% 
%     
%     
%     
    
    figure(1015);clf();plot(fits(:,:,3)'); hold on; plot(nanmedian(fits(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');ylabel('tau 1')
    
    tau1 = fits(:,:,3);
    if isnan(tau1) % only when no even was detected
        tau1 = [];
    end
    figure(1016);clf();scatter(all_peak_values(:),tau1(:)); hold on; title('tau vs event amplitude'); xlabel('Amplitude'); set(gcf,'Color','w');
end