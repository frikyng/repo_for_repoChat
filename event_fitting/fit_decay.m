%% Fit event decay. Return NaN for the most problematic cases
function [p, x_fit, y_fit, bkp] = fit_decay(trace, pre_peak_delay, fit_params, rendering, jitter_win, peak_av_width, other_peaks)
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
    if isempty(fit_params) && rendering
        figure(1019);cla(); hold on;
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
    peak = trace(pre_peak_delay -jitter_win(1):pre_peak_delay+jitter_win(end));
    if rendering
        plot(t_ax, trace); hold on; % plot before cropping
    end
    
    %% Try fitting process. If any failure, return NaNs
    %% If you have multiple peaks in a burst, there is a change that you fitted the wrong event. That's ok and handled lower down
    %% Now fit trace
    [x_fit, y_fit, p, detected_peak] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, nanmax(peak), '', rendering, peak_av_width);
    
    if isempty(x_fit)
        x_fit = [];
        y_fit = [];
        p     = NaN(1, 5);
        bkp   = p;
    else
        colortype = 'b-';

        %% Check quality.
        [failed, error_estimate, anchor_point] = quality_check(trace, x_fit, y_fit, other_peaks);
        bkp         = [];
        %% Show issue if any
        if rendering && any(failed)
            %% Show bad fit
            recover = plot(x_fit, y_fit, 'k--'); hold on; 
            if numel(error_estimate) > 1
                hold on;scatter(x_fit, error_estimate, 'k.');
            end
        end
        
        %% If slope is too shallow, add weights
        if any(ismember(failed,3))             
            [x_fit, y_fit, p, detected_peak] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, nanmax(peak), x_fit(~isnan(error_estimate)), rendering, peak_av_width);
            if rendering
                recover.XData = x_fit;
                recover.YData = y_fit;
                colortype = 'g-';
            end
            
            if isempty(x_fit)
                p     = NaN(1, 5);
                bkp   = p;
                return
            end
        end
        
        
        %% If we did cut through, we'll fit the first clean decay that follows, and then shift the trace (if there are neough points for a fit_
        if any(ismember(failed,2)) && ((pre_peak_delay + anchor_point) < (t_ax(end) - 10))
            colortype = 'r-';
            
            [x_fit, y_fit, p, ~] = fit_trace(t_ax, trace, pre_peak_delay + anchor_point, jitter_win, fit_params, nanmax(peak), x_fit(~isnan(error_estimate)), rendering, peak_av_width);
            if rendering
                recover.XData = x_fit;
                recover.YData = y_fit;
            end
            
            if isempty(x_fit)
                p     = NaN(1, 5);
                bkp   = p;
                return
            end
            
            %% Shift trace to align with corrrect peak
            x_off = (x_fit(1)-detected_peak);
            x_fit = x_fit - x_off;
            x_mid_point_1 = find(y_fit < y_fit(1)/2,1,'first')+x_off+detected_peak;
            if isempty(x_mid_point_1)
                x_mid_point_1 = y_fit(end);
            end
            mid_point_1 = [x_mid_point_1  , y_fit(1)/2];                        % First Point

            %% Adjust weight to get correct peak amplitude
            p([2, 4]) = p([2, 4]) * (nanmax(peak) / (p(2) + p(4)));
            t0 = x_fit(1);
            y_fit = p(1)+p(2)*exp(-((x_fit - t0)*p(3))) + p(4)*exp(-((x_fit - t0)*p(5)));
            x_mid_point_2 = find(y_fit < y_fit(1)/2,1,'first')+detected_peak;
            if isempty(x_mid_point_2)
                x_mid_point_2 = y_fit(end);
            end
            mid_point_2 = [x_mid_point_2 , y_fit(1)/2];                         % Second Point

            if rendering
                mid_point_1(2)          = mid_point_2(2);                
                diff_mid_points         = mid_point_2 - mid_point_1;                         % Difference
                hold on;quiver(mid_point_1(1),mid_point_1(2),diff_mid_points(1),diff_mid_points(2),0, 'LineWidth', 1, 'Color', 'k');
            end

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
end

function [x_fit, y_fit, p, pre_peak_delay] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, peak, extra_w, rendering, peak_av_width)   
    y     = double(trace);
    x_fit = [];
    y_fit = [];
    p     = NaN(1, 5);
    if nargin >= 8 && ~isempty(extra_w)
        extra_w = [t_ax(extra_w) , y(extra_w)];
    else
        extra_w = [];
    end
    
    %% Clip trace if any breaking (NaN) point is found (before or after peak), since there could be huge time gap between the 2
    if all(isnan(trace))
        return
    end
    break_point = find(isnan(trace(1:pre_peak_delay)));
    if any(break_point)
        t_ax = t_ax((break_point(end)+1):end);
        y    = y((break_point(end)+1):end);
    end
    break_point = find(isnan(y));
    if any(break_point)
        t_ax = t_ax(1:break_point(1)-1);
        y    = y(1:break_point(1)-1);
    end

    %% Clear all data before peak
    test_peak   = y;
    test_peak([1:pre_peak_delay-jitter_win(1),pre_peak_delay+jitter_win(2):end]) = NaN;
    [~, locs]   = findpeaks(test_peak,'MinPeakWidth', 2);
    if isempty(locs)
        [~, pre_peak_delay] = nanmax(test_peak); % make sure we have the local max. This prevents further complications
    else
        pre_peak_delay = locs(end);
    end
    fit_range   = pre_peak_delay:numel(t_ax);
    x           = t_ax(fit_range);
    y           = y(fit_range);
    
    %% Remove the first 20% of decay to help fit
    eighty  = max(y) - (20 * max(y)/100);
    eighty  = find(y < eighty,1,'first');
    x       = x(eighty:end);
    y       = y(eighty:end);
    
    %% If the peak was very close to the end of trial, we may have a problem here. We interupt fit
    if isempty(x) || numel(x) < 3
        return 
    end

    %% Remove any secondary peak in the decay (clip at inflection point)
    thr = y(1)/10;
    if thr < 0
       thr = rms(diff(y)) * 10;
    end
    [val, secondary_peak] = findpeaks(y*-1, x, 'MinPeakProminence', thr);
    if ~isempty(secondary_peak)
        y = y(1:find(x == secondary_peak(1)));
        x = x(1:find(x == secondary_peak(1)));
        if rendering
            plot(x, y, 'r'); hold on % plot what's used
        end

        %% Add trace lowest point since the fit should go through that one (prevent underfitting)
        val = val * -1;
        [~, low] = min(val);
        x = [x; secondary_peak(low)];
        y = [y; val(low)];
        if rendering
            scatter(secondary_peak(low),  val(low), 'filled')
        end
        
%         win = (secondary_peak(low)-2):(secondary_peak(low)+2);
%         win(win < 1) = [];
%         win(win > x(end)) = [];
%         win = ismember(x, win);
%         y = [y; nanmean(y(win))]; % local minimum
%         x = [x; secondary_peak(low)];
%         if rendering
%             scatter(secondary_peak(low),  y(end), 'filled')
%         end
        
        
        %% Give stronger weight to the end of the fit region, and any extra point
        w = ones(size(y));
        w(end) = numel(trace);
    else
        w = ones(size(y));
        if rendering
            plot(x, y, 'r'); hold on % plot what's used
        end
    end
    
    if ~isempty(extra_w)
        x = [x; extra_w(:,1)];
        y = [y; extra_w(:,2)];
        w(end+1:end+size(extra_w, 1)) = numel(w);
        robust = 'Bisquare';
    else
        robust = 'Off';
    end
    
    %% Adjust contraints for t0 and x0 (decay baseline)
    t0 = x(1);
    x0 = 0;
    
    %% Now fit from scratch or use params and scale (if there are neought points
    if (isempty(fit_params) || numel(x) >= 4) && ~any(isnan(fit_params))
        if isempty(fit_params) % free robust fit
            f = fit((x-t0), y, 'exp2',...
                'Lower',[0,-Inf,0,-Inf],...
                'Upper',[peak,-1e-3,peak,-1e-3],...
                'Weights',w,...
                'Robust','Bisquare');
        else % guided fit
            if peak < peak2peak(y), peak = peak2peak(y); end
            fit_params_low      = [0,-10/peak_av_width,0,-10/peak_av_width]; % (-1/peak_av_width) width is a reasonable starting approximation for tau (i.e 67% decay)
            fit_params_high     = [peak,-0.1/peak_av_width,peak,-0.1/peak_av_width];
   
            f = fit((x-t0),y,'exp2',...
                'Lower',fit_params_low,...
                'Upper',fit_params_high,...
                'Weights',w,...
                'StartPoint',fit_params(2:end),...
                'Robust',robust);
        end
        
        %% Make sure we have fast tau first, slow tau second
        if f.b > f.d
            p = [x0, f.c, f.d*-1, f.a, f.b*-1];
        else
            p = [x0, f.a, f.b*-1, f.c, f.d*-1];
        end 
    elseif numel(x) >= 2  && ~any(isnan(fit_params))% minimal fit
        f = fit((x-t0),y,'exp1','Lower',[0,-Inf],'Upper',[max(0, peak),-1e-3],'Weights',w,'StartPoint',fit_params(2:3));% ,'Robust','Bisquare'
        p = [x0, f.a, f.b*-1, 0, 0];
    else
        return
    end

    %% Generate fit curve    
    x_fit = (x(1):t_ax(end))';
    y_fit = p(1)+p(2)*exp(-((x_fit - t0)*p(3))) + p(4)*exp(-((x_fit - t0)*p(5)));
    % hold on;plot(x_fit, y_fit)
    % hold on;plot(p(2)*exp(-((x_fit - t0)*p(3))), 'k'); hold on; plot(p(4)*exp(-((x_fit - t0)*p(5))),'r'); % to plot tau 1 and tau 2
end

function [failed, error_estimate, half_decay] = quality_check(trace, x_fit, y_fit, other_peaks)
    failed          = 0;
    error_estimate  = [];
    half_decay      = find(y_fit < y_fit(1)/2, 1, 'first');
    
    %% Check for empty fit
    if isempty(x_fit)
        failed  = 1;
        return
    end
    
    if x_fit(end) > numel(trace)
        y_fit = y_fit(x_fit < numel(trace));
        x_fit = x_fit(x_fit < numel(trace));
    end
    
    %% Compute difference between fit and trace
    error_estimate     = trace(x_fit(1):x_fit(end)) - y_fit; %difference between fit and trace
    thr     = -nanmax(y_fit)/20; % loor for disprencies > 5% peak
    error_estimate(error_estimate > thr) = NaN;

    if ~isempty(other_peaks)
        other_peaks = [];
        next_peak = [];%other_peaks(find(other_peaks > x_fit(1), 1, 'first'));
    end
    
    %% Check if fit cut through signal (wrong peak estimate, need offset)
    if sum(trace(x_fit(1):x_fit(half_decay)) - y_fit(1:half_decay)) > (y_fit(1)/2)
        %% Walk along the slope until we reach the orginal peak height if possible
        last = min((x_fit(1) + half_decay), x_fit(end));
        better = find(trace(x_fit(1):last) > y_fit(1), 1, 'last');
        if ~isempty(better)
            half_decay = better;
        end
        if isempty(other_peaks) || last < next_peak
            failed  = [failed, 2];
        end
    end
    
    %% Check if at least 5 bad conscutive points ("underfit, curve to shallow")
    win = 10;
    if any(~isnan(movmean(error_estimate, win,'Endpoints', 'discard'))) 
        failed  = [failed, 3];
    end
end
