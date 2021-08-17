function [all_peak_values, fit_summary] = fit_all_events(events, all_data, peak_locs, pre_peak_delay, fit_summary, rendering, guesses, peak_av_width, fit_type, forced_tau, jitter_win)
    if nargin < 7 || isempty(guesses)
        guesses = [];
    end
    if nargin < 8 || isempty(guesses)
        peak_av_width = [];
    end
    if nargin < 9 || isempty(fit_type)
        fit_type = 'exp1';
    end
    if nargin < 10 || isempty(forced_tau)
        forced_tau = false;
    end
    if nargin < 11 || isempty(jitter_win)
        jitter_win = [2, 4];
    end
    
    spmd
        warning('off')
    end    
    if rendering
        figure(123);cla(); % temp figure
    end 

    if isempty(events)
        all_peak_values    = [];
        fit_summary.events = {NaN(1,3,1)};
        fit_summary.window = [];
        fit_summary.t_peak = [];
        return
    end
    
    all_peak_values = all_data(peak_locs, :);
    all_peak_locs   = repmat(peak_locs',1,size(all_data,2));
    event_fit       = cell(size(events, 3), 1);
    ncores          = (feature('numcores') * double(~rendering));
    temp_val        = cell(1, size(events, 3));
    temp_t          = cell(1, size(events, 3));
    
rendering = 0
forced_tau = -1
   parfor (peak = 1:size(events, 3), ncores)     
   %for peak = 1:size(events, 3)   
        fit_params{peak}            = [];
        t_peak_offset{peak}         = []; % relative to the median peak time
        
        %% Get median peak loc to help with bad fits QQ (could be passed as an input)
        median_event        = nanmedian(events(:, : ,peak), 2);
        [pks, locs] = findpeaks(median_event);
        
        difference = abs(locs - pre_peak_delay);
        [m, closest] = min(difference);
        if sum(difference == m) > 1 && locs(closest) < pre_peak_delay %sometimes, you may detect 2 peaks at the same distace. Then we'll favour the one happening after peak time
            closest = closest + 1;
        end

        median_event_loc    = locs(closest);
        median_event_amp    = median_event(median_event_loc);

        for gp = 1:size(events, 2)
            if rendering
                figure(123);cla(); % temp figure
                plot([pre_peak_delay, pre_peak_delay],[0, max(pks)], 'k--'); hold on;
            end

            %% Get contraints for peak search (between previous and next peak, or trace limits)
            [prev, next]    = get_limits(peak, peak_locs, all_data(:, gp)); 
            
            %% Adjust search window                
            peak_win = adjust_limits(jitter_win, prev, next, pre_peak_delay);

            %% Fit                
            [flag, x_fit_axis, y_fit_axis, bkp, detected_peak_loc] = fit_decay(events(:, gp ,peak), pre_peak_delay, guesses(gp, :), rendering, peak_win, peak_av_width, peak_locs - peak_locs(peak) + pre_peak_delay, fit_type, true);   

            %% If fit failed despite our attempts (usually because of the absence of peak or really strange kinetics, we set correction to 0
            if any(isnan(flag)) && ~any(isnan(flag))
                flag    = bkp;
            elseif any(isnan(flag))
                x_fit_axis       = (1:numel(events(:, gp ,peak)))';
                y_fit_axis       = zeros(size(x_fit_axis));
            elseif strcmp(fit_type, 'exp2') && numel(flag) == 3
                flag    = [flag, 0, 0];   
            end

            %% Real location of the event in the whole trace (in timepoints)
            event_t_start       = peak_locs(peak)-pre_peak_delay - 1; % real t_start of the event fitting segment
            t_fit               = x_fit_axis + event_t_start;
            y_fit_axis(end)     = NaN;
            event_fit{peak}{gp} = [t_fit, y_fit_axis];
            if isnan(flag)
                flag = NaN(1, 5);
            end
            fit_params{peak}    = [fit_params{peak} ; flag];
            
            %% Now, the fit doesn't necessarily start at the peak time, so we need to find that one
            t_peak_offset{peak} = [t_peak_offset{peak} ; detected_peak_loc - pre_peak_delay]; % t_peak in window (after t_peak delay)
            temp_val{peak}{gp}  = events(detected_peak_loc, gp ,peak); % ultimately all_peak_values(peak, gp)
            temp_t{peak}{gp}    = detected_peak_loc - pre_peak_delay;
            
            if rendering     
                rescaled_median_event = (median_event / median_event_amp) * temp_val{peak}{gp};
                scatter(median_event_loc, rescaled_median_event(median_event_loc),300, 'v', 'filled','MarkerFaceColor',[0.6,0.6,0.6]);
                scatter(temp_t{peak}{gp} + pre_peak_delay, temp_val{peak}{gp},300, 'v', 'MarkerFaceColor','b')
                hold on;plot(rescaled_median_event, '--','Color',[0.6,0.6,0.6]);
                drawnow;
            end
        end
    end
    
    %% Can't put it in parfor
    for peak = 1:size(events, 3)
        for gp = 1:size(events, 2)
            all_peak_values(peak, gp) = temp_val{peak}{gp};
            %all_peak_times(peak, gp) = temp_t{peak}{gp};
        end
    end
    
    fit_summary.events              = fit_params;
    fit_summary.window              = vertcat(event_fit{:, 1});
    fit_summary.peak_jitter         = t_peak_offset;
end

function [prev, next] = get_limits(peak, peak_times, all_data)
    tp = size(all_data,1);
    %% Identify next peak location make sure not to fit beyond
    if peak < numel(peak_times)
        [~, low_point] = nanmin(all_data(peak_times(peak):peak_times(peak+1), :)); %if there's a peak following, don't go further han the lowest point between the 2 peaks
        next = low_point;%diff([peak_times(peak),peak_times(peak)+low_point]);
        if low_point < 10
            %% You need enough points for fitting. if you don't have them, wel'll use another strategy
            next = diff([peak_times(peak), tp]);
        end
    else
        next = diff([peak_times(peak), tp]);
    end

    %% Make sure we don't use points before 1 
    if peak > 1
        prev = diff(peak_times(peak-1:peak));
    else
        prev = diff([1, peak_times(peak)]);
    end
end

function peak_win = adjust_limits(peak_win, prev, next, pre_peak_delay)   
    if peak_win(1) >= pre_peak_delay
        peak_win(1) = pre_peak_delay - 1;
    end
    
    if peak_win(1) >= prev
        peak_win(1) = prev - 2;
    end
    if peak_win(2) >= next
        peak_win(2) = next - 2;
    end
end