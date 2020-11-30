function [event_fit, all_peak_values, fit_summary] = fit_all_events(events, all_data, peak_locs, pre_peak_delay, fit_summary, rendering, guesses, peak_av_width)
    if nargin < 7 || isempty(guesses)
        guesses = repmat(fit_summary.global_fit,size(events, 2),1);
    end
    if nargin < 8 || isempty(guesses)
        peak_av_width = [];
    end
    
    if isempty(events)
        [event_fit, all_peak_values] = deal([]);
        fit_summary.events = {NaN(1,3,1)};
        fit_summary.window = [];
        fit_summary.t_peak = [];
        return
    end

    all_peak_values = all_data(peak_locs, :);
    all_peak_locs   = repmat(peak_locs',1,size(all_data,2));
    event_fit       = cell(size(events, 3), 1);
    ncores          = (feature('numcores') * double(~rendering));
    temp_val = cell(1, size(events, 3));
    temp_t = cell(1, size(events, 3));
    %parfor (peak = 1:size(events, 3), ncores) %fliplr(event_to_correct)
    for peak = 1:size(events, 3)
        fit_params{peak}    = [];
        t_peak{peak}        = [];
        
        %% Get median peak loc to help with bad fits QQ (could be passed as an input)
        median_event        = nanmedian(events(:, : ,peak), 2);
        [pks, locs] = findpeaks(median_event);
        [~, closest] = min(abs(locs - pre_peak_delay));
        Mloc = locs(closest);
        M    = median_event(Mloc);

        for gp = 1:size(events, 2)
            if rendering
                figure(123);cla();
                plot([pre_peak_delay, pre_peak_delay],[0, max(pks)], 'k--'); hold on;
            end

            %% Get contraints for peak search (between previous and next peak, or trace limits)
            [prev, next]    = get_limits(peak, peak_locs, all_data(:, gp)); 
            
            %% Adjust search window                
            peak_win = adjust_limits([4, 4], prev, next, pre_peak_delay);

            %% Fit                
            [flag, x, y, bkp] = fit_decay(events(:, gp ,peak), pre_peak_delay, guesses(gp, :), rendering, peak_win, peak_av_width, peak_locs - peak_locs(peak) + pre_peak_delay);%Mloc, );      

            %% If fit failed despite our attempts (usually because of the absense of peak or really strange kinetics, we set correction to 0
            if any(isnan(flag)) && ~any(isnan(flag))
                flag    = bkp;
            elseif any(isnan(flag))
                x       = (1:numel(events(:, gp ,peak)))';
                y       = zeros(size(x));
            end

            t = x + (peak_locs(peak)-pre_peak_delay - 1);

            y(end)  = NaN;
            event_fit{peak}{gp} = [t, y];
            if isnan(flag)
                flag = NaN(1, 5);
            end
            fit_params{peak} = [fit_params{peak} ; flag];
            t_peak{peak} = [t_peak{peak} ; x(1) - pre_peak_delay];
            [M2, M2loc] = nanmax(events(pre_peak_delay-3:pre_peak_delay+3, gp ,peak));  
            temp_val{peak}{gp} = M2;
            temp_t{peak}{gp} = M2loc;
            if rendering     
                temp2 = (median_event / M) * M2;
                scatter(Mloc, temp2(Mloc), 'v', 'filled','MarkerFaceColor',[0.6,0.6,0.6]);
                scatter(x(1), y(1), 'v', 'filled','MarkerFaceColor','b')
                hold on;plot(temp2, '--','Color',[0.6,0.6,0.6]);
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
    
    event_fit = vertcat(event_fit{:, 1});    
    fit_summary.events = fit_params;
    fit_summary.window = event_fit;
    fit_summary.t_peak = t_peak;
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