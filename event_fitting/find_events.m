function [fit_data] = find_events(all_data, peak_thr, global_timescale)
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
end

