function fit_data = detect_and_fit_events(all_data, global_timescale, rendering, group_labels, peak_thr)
    if nargin < 3 || isempty(rendering)
        rendering = true;
    end
    if nargin < 5 || isempty(peak_thr)
        peak_thr = 1
    end
    
    fit_data                = {};
    fit_data.events         = {};
    fit_data.group_labels   = group_labels;
    refit                   = true
    gcamp_tau               = 0.2
    fit_data.fast_fit       = true
    rendering = 1;
    
%% qq bring fit constraints up here
%     
       %all_data        = all_data(1:300,:);
     %  global_timescale = global_timescale(1:300);
    all_data        = all_data - prctile(nanmean(all_data'),1);
    
    %% Fitting adjustements with low SR can be a bit funny when cells are bursting, but upsampling data a bit does the job
    if (nanmedian(diff(global_timescale))/2) < gcamp_tau
        original_data =  all_data;
        original_timescale =  global_timescale;       
        resampling = ceil(gcamp_tau / (nanmedian(diff(global_timescale))/2));
    else
        resampling = 1;
    end
    fit_data.resampling_factor      = resampling;
    
    
    %% Filter signal (won't work if signal is already upsampled)
    global_median       = nanmedian(all_data, 2);
    [denoised, missing] = wavelet_denoise(global_median);    
    bsl                 = global_median(~missing) - denoised(~missing);
    %thr                 = rms(bsl) * thr_factor;    
    %figure(25);cla();envelope(bsl,20,'peak');
    [a, ~] = envelope(bsl,20,'peak');
    peak_thr = rms(a)*peak_thr;    
    %thr = (prctile(a-b,99) - prctile(a-b,1))*peak_thr;
    fprintf('PLEASE REVIEW THRESHOLDING CRITERIA BEFORE FINAL ANALYSIS\n') 
    denoised(missing)   = NaN; % restore NaNs

    if resampling > 1
        global_timescale    = interpolate_to(global_timescale, numel(global_timescale) * resampling); 
        denoised            = interpolate_to(denoised        , numel(global_timescale));
        try
            all_data            = interpolate_to(all_data, numel(global_timescale), 'cubic')'; % failed with 2019-11-06_exp_1  
        catch
            all_data            = interpolate_to(all_data, numel(global_timescale))';                
        end
    end

    %% Find peaks using wavelet denoising on median trace
    [pk_val, peak_loc, peak_times, widths] = detect_events(denoised', global_timescale', peak_thr);

    fit_data.peak_times             = peak_times{1}';
    fit_data.peak_pos               = peak_loc{1}';

    %% Use peak width to help with fitting
    peak_av_width   = nanmedian([widths{:}]); % that's the average peak width at half height, in points (median is in case there are a few trials with avery different sr
    
    % Based on it, we define a reasonable fitting window
    pre_peak_delay  = round(peak_av_width * 3);
    post_peak_delay = round(peak_av_width * 15);

    %% We'll do the fitting on denoised data
    [all_data, missing] = wavelet_denoise(all_data);

    %% Extract all events
    if isempty(fit_data.peak_pos)
        fit_data.peak_pos = NaN; % to avoid default detection
    end
    events = extract_events(all_data, fit_data.peak_pos, pre_peak_delay, post_peak_delay);
    
    %% Get global fit settings to have a rough estimate of the fitting window required. Use small, isolated events (or end-of-burst)
    to_use = find(diff(fit_data.peak_pos) > (peak_av_width * 5)); % keep events spaced by at least 5 average peak width from the next, to avoid contamination of decay
    figure(1019);cla();hold on;plot(reshape(events,size(events,1),[]), 'Color',[0.9,0.9,0.9]); hold on;title('Individual median events amd median fit');hold on;
    fit_data.global_fit = fit_decay(squeeze(nanmean(nanmedian(events(1:round(peak_av_width * 5),:,to_use), 2),3)), pre_peak_delay, [], true, [], peak_av_width, [], 'exp2'); % the fitting window can be shorter as we don't have bursts in this subselection
    xlabel('points');legend({'mean event','fit result'}); hold on;
    arrangefigures([1,2]); 

    %     if rendering
    %         show_events(all_data, events, pre_peak_delay)
    %     end
    rendering = false
    %% In fast_fit mode, we force tau and don't need a refit
    if fit_data.fast_fit
        guesses     = NaN(size(events, 2), size(fit_data.global_fit,2));
        for gp = 1:size(events, 2)
            guesses(gp, :) = fit_decay(squeeze(nanmean(events(1:round(peak_av_width * 7),gp,to_use),3)), pre_peak_delay, [], false, [], peak_av_width, [], 'exp2');
        end
        forced_tau  = true;
        refit       = false;        
    else
        guesses     = repmat(fit_data.global_fit,size(events, 2),1);
        forced_tau  = false;        
    end
    
    %% Now fit all events with fairly free parameters.
    % initial guesses is set by the median tau 
    upsampled_sr = 1/median(diff(global_timescale));
    jitter_win = round([0.2*upsampled_sr,0.2*upsampled_sr]); %200 ms jitter
    
    [all_peak_values, fit_data] = fit_all_events(events, all_data, fit_data.peak_pos, pre_peak_delay, fit_data, rendering, guesses, peak_av_width, 'exp2', forced_tau, jitter_win);

    %% Make some figures about fit propeties
    fits_free = 1/permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
    figure(1014);clf();plot(fits_free(:,:,3)) ; hold on; plot(nanmedian(fits_free(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');xticklabels(group_labels);xtickangle(45);set(gcf,'Color','w');
    figure(1015);clf();plot(fits_free(:,:,3)'); hold on; plot(nanmedian(fits_free(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');
    fit_data.pre_refit_events = fit_data.events;
    arrangefigures([1,2]); 

    %% Re-Fit events using constraint for their group
    fit_data.pre_correction_peaks   = all_peak_values;
    if refit
        %% Now refit but keep tau values within reason (+/- 1 mad)
        test = permute(cat(3, fit_data.events{:}),[1,3,2]);
        tau1 = nanmedian(test(:,:,3),2);
        A1   = nanmedian(test(:,:,2),2);
        if size(test, 3) > 3
            tau2 = nanmedian(test(:,:,5),2);
            A2   = nanmedian(test(:,:,4),2);
            guesses = [zeros(size(tau1)), A1, tau1, A2, tau2];
        else
            guesses = [zeros(size(tau1)), A1, tau1];
        end

        %show_events(all_data, events)
        [all_peak_values, fit_data] = fit_all_events(events, all_data, fit_data.peak_pos, pre_peak_delay, fit_data, rendering, guesses, peak_av_width, 'exp1',true);
        refitted_events = 1/permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
        figure(1014);clf();plot(refitted_events(:,:,3)) ; hold on; plot(nanmedian(refitted_events(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');xticklabels(group_labels);xtickangle(45);set(gcf,'Color','w');
        figure(1015);clf();plot(refitted_events(:,:,3)'); hold on; plot(nanmedian(refitted_events(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');  
        fit_data.pre_correction_peaks   = all_peak_values;
    end
    
    %% Correct the traces before downsampling
    %[all_peak_values, correction]   = correct_traces(all_data, events, fit_data.window, all_peak_values, fit_data.peak_pos, global_timescale);
    fit_data.post_correction_peaks  = all_peak_values;

    %% Revert the effect of resampling (so the timepoints match the real trace)
    if fit_data.resampling_factor ~= 1
        all_data = original_data;
        fit_data = fix_resampling(fit_data);
    end
    
    %% Finally, correct the events amplitudes and store results    


    %if rendering
    show_results(fit_data, all_data, group_labels, original_timescale)
    %end
end

function fit_data = fix_resampling(fit_data)
    %% If the initial trace was resampled, some elements needs to be readjusted
    fit_data.peak_pos = round(fit_data.peak_pos / fit_data.resampling_factor);
    for el = 1:numel(fit_data.window)
        fit_data.window{el} = interpolate_to(fit_data.window{el}, ceil(size(fit_data.window{el},1)/fit_data.resampling_factor))';
        fit_data.window{el}(:,1) = round(fit_data.window{el}(:,1)/fit_data.resampling_factor);
    end
    % .t_peak % ?
end

function show_results(fit_data, traces, group_labels, global_timescale)
    R                             = max(range(traces));
    peak_values_raw             = fit_data.pre_correction_peaks;
    peak_values_corrected       = fit_data.post_correction_peaks;
    for gp = 1:size(peak_values_corrected, 2)
        figure(10050 + gp);cla();hold on;title(['Traces and fits group ',num2str(gp),' (',group_labels{gp},')']);set(gcf,'Color','w');
        plot(global_timescale, traces(:, gp),'k');hold on;ylim([-R/20,R + R/20])
        scatter(fit_data.peak_times, nanmedian(peak_values_raw(:,gp), 2), 'bv', 'filled');hold on;
        scatter(fit_data.peak_times, nanmedian(peak_values_corrected(:,gp), 2), 'kv', 'filled');hold on;
        for el = 1:size(fit_data.window, 1)
            tp = fit_data.window{el,gp}(:,1);
            tp = tp > 0 & tp <= size(global_timescale, 2); 
            plot(global_timescale(fit_data.window{el,gp}(tp,1)),fit_data.window{el,gp}(tp,2),'r--'); hold on;
        end
        drawnow
    end

    %figure(1002);hold on;scatter(peak_times, nanmean(all_peak_values, 2), 'bo', 'filled');hold on;plot(global_timescale, nanmean(correction, 2),'r--');set(gcf,'Color','w');

    %% Make some figures about fit properties
    fits = permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
    figure(1014);clf();plot(fits(:,:,3)) ; hold on; plot(nanmedian(fits(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');set(gcf,'Color','w');ylabel('tau 1')%xticklabels(bin_legend);xtickangle(45)
        
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
    
    figure(1015);clf();plot(fits(:,:,3)'); hold on; plot(nanmedian(fits(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');ylabel('tau 1')
    
    tau1 = fits(:,:,3);
    if isnan(tau1) % only when no even was detected
        tau1 = [];
    end
    figure(1016);clf();scatter(fit_data.post_correction_peaks(:),tau1(:)); hold on; title('tau vs event amplitude'); xlabel('Amplitude'); set(gcf,'Color','w');
end