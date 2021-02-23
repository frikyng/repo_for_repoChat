%% Fit event decay. Return NaN for the most problematic cases
function [p, x_fit, y_fit, bkp, detected_peak] = fit_decay(trace, pre_peak_delay, fit_params, rendering, jitter_win, peak_av_width, other_peaks, fit_type, forced_tau)
    if nargin < 3 || isempty(fit_params)
        fit_params = [];
    end
    if nargin < 4 || isempty(rendering)
        rendering = true;
    end 
    if nargin < 5 || isempty(jitter_win)
        jitter_win = [2,4];
    end    
    if nargin < 6 || isempty(peak_av_width)
        peak_av_width = [];
    end  
    if nargin < 7 || isempty(other_peaks)
        other_peaks = [];
    end 
    if nargin < 8 || isempty(fit_type)
        fit_type = 'exp1';
    end 
    if nargin < 9 || isempty(forced_tau)
        forced_tau = false;
    end     
    if isempty(fit_params) && rendering
        figure(1019);hold on;
    end
    
    %% When there was no detected events
    if isempty(trace)
        [p, x_fit, y_fit, bkp] = deal([]);
        return
    end

    %% Prepre fit options
    persistent opts
    if isempty(opts)
       opts = optimoptions('lsqcurvefit','Display','off');
    end
    
    %% Prepare variables
    t_ax = (1:numel(trace))';
    peak = trace(pre_peak_delay - jitter_win(1):pre_peak_delay+jitter_win(end));
    if rendering
        plot(t_ax, trace); hold on; % plot before cropping
    end
    
    %% Try fitting process. If any failure, return NaNs
    %% If you have multiple peaks in a burst, there is a chance that you fitted the wrong event. That's ok and is handled lower down
    %% Now fit trace
    [x_fit, y_fit, p, detected_peak] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, nanmax(peak), '', rendering, peak_av_width, fit_type, forced_tau && ~forced_tau == -1);
    
    if isempty(x_fit)
        x_fit = [];
        y_fit = [];
        p     = default_p(fit_type);
        bkp   = p;
    elseif ~forced_tau || forced_tau ==-1
        colortype = 'b-';

        %% Check quality.
        [failed, error_estimate, anchor_point] = quality_check(trace, x_fit, y_fit, other_peaks);
        
        if any(ismember(failed, 4))
            failed = 0;
            colortype = 'm';
        end
        
        
        %% For the quick correction mode
        if forced_tau == -1
            if any(failed)
                [x_fit, y_fit, p, detected_peak] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, nanmax(peak), '', rendering, peak_av_width, fit_type, true);
                colortype = 'k-';
                bkp   = p;
                if (x_fit(1) - pre_peak_delay) > jitter_win(2)
                    corr = x_fit(1) - (detected_peak) - 4; %4pt is arbitrary - we need to find the delay between peak and decay
                    x_fit = x_fit - corr;
                    colortype = 'r-';
                end
                if rendering
                    plot(x_fit,y_fit, colortype, 'LineWidth',2);
                end
                return
            else
                bkp   = p;
                if (x_fit(1) - pre_peak_delay) > jitter_win(2)
                    corr = x_fit(1) - (detected_peak) - 4; %4pt is arbitrary - we need to find the delay between peak and decay
                    x_fit = x_fit - corr;
                    colortype = 'r-';
                end
                if rendering
                    plot(x_fit,y_fit, colortype, 'LineWidth',2);
                end 
                return
            end
        else
            
        
            bkp         = [];
            %% Show issue if any
            if rendering && any(failed)
                %% Show bad fit
                recover = plot(x_fit, y_fit, 'k--'); hold on; 
                if numel(error_estimate) > 1
                    hold on;scatter(x_fit, error_estimate, 'kx');
                end
            end

            %% If slope is too shallow, add weights
            if any(ismember(failed,3))             
                [x_fit, y_fit, p, detected_peak] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, nanmax(peak), x_fit(~isnan(error_estimate)), rendering, peak_av_width, fit_type, forced_tau);
                if rendering
                    recover.XData = x_fit;
                    recover.YData = y_fit;
                    colortype = 'g-';
                end

                if isempty(x_fit)
                    p       = default_p(fit_type);
                    bkp     = p;
                    return
                end

                %% Now, rechecking
                [failed, error_estimate, anchor_point] = quality_check(trace, x_fit, y_fit, other_peaks);
                if rendering
                    plot(x_fit,y_fit,colortype, 'LineWidth',2);
                end
            end

            %% If we did cut through the peak, we'll fit the first clean decay that follows, and then shift the trace (if there are enough points for a fit).
            if any(ismember(failed,2)) && ((pre_peak_delay + anchor_point) < (t_ax(end) - 10))
                colortype = 'r-';

                x_off = (x_fit(1)-detected_peak);
                if ~x_off
                    x_off = 1;
                end

                range   = 5:-1:-anchor_point; % qq why 5
                err     = NaN(size(range));
                parfor a = 1:numel(range)
                    pt = range(a);
                    [x_fit_temp, y_fit, p, ~, x, y] = fit_trace(t_ax, trace, pre_peak_delay + anchor_point + pt, jitter_win, fit_params, nanmax(peak), x_fit(~isnan(error_estimate)), rendering, peak_av_width, fit_type, forced_tau);
                    if ~isempty(x)
                        plot(x_fit_temp,y_fit,'r--', 'LineWidth',0.5);                
                        x_fit_range     = detected_peak:x(end-1);
                        if numel(p) == 3
                            y_fit_range     = p(1)+p(2)*exp(-((x_fit_range - x_fit_range(1))*p(3)));
                        elseif numel(p) == 5
                            y_fit_range     = p(1)+p(2)*exp(-((x_fit_range - x_fit_range(1))*p(3))) +p(4)*exp(-((x_fit_range - x_fit_range(1))*p(5)));
                        end
                        trace_for_fit   = trace(x_fit_range)';
                        err(a)          = mean(abs(trace_for_fit - y_fit_range));
                    else
                        err(a)          = NaN;
                    end
                end
                [~, a] = nanmin(err);
                [x_fit, y_fit, p, ~] = fit_trace(t_ax, trace, pre_peak_delay + anchor_point + range(a), jitter_win, fit_params, nanmax(peak), x_fit(~isnan(error_estimate)), rendering, peak_av_width, fit_type, forced_tau);
                plot(x_fit,y_fit,'r--', 'LineWidth',3);
            end

            [failed, error_estimate, anchor_point] = quality_check(trace, x_fit, y_fit, other_peaks);

            %% Check if we ended fitting the correct stuff
            if ~isempty(bkp)
                p = bkp;
            end

            if rendering
                plot(x_fit,y_fit,colortype, 'LineWidth',2);
            end
        end
    else
        colortype = 'b-';
        if (x_fit(1) - pre_peak_delay) > jitter_win(2)
            corr = x_fit(1) - (detected_peak) - 4; %4pt is arbitrary - we need to find the delay between peak and decay
            x_fit = x_fit - corr;
            colortype = 'r-';
        end
        bkp         = [];
        if rendering
            plot(x_fit,y_fit,colortype, 'LineWidth',2);
        end
    end
end


function p = default_p(fit_type)
    if strcmp(fit_type, 'exp2')
        p       = NaN(1, 5);
    elseif strcmp(fit_type, 'exp1')
        p       = NaN(1, 3);
    end
end