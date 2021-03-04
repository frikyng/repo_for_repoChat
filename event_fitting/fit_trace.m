function [x_fit, y_fit, p, pre_peak_delay, x_used, y_used] = fit_trace(t_ax, trace, pre_peak_delay, jitter_win, fit_params, peak, extra_w, rendering, peak_av_width, fit_type, forced_tau)   
    [x_fit, y_fit, x_used, y_used]  = deal([]);    
    p                               = default_p(fit_type);
    if nargin >= 8 && ~isempty(extra_w)
        extra_w = [t_ax(extra_w) , double(trace(extra_w))];
    else
        extra_w = [];
    end
    
    %rendering = 1;
    
    %% Clip trace if any breaking point is found (i.e. NaN before or after peak)
    % There no guarentee of timscale continuity across these breaking points, so we should't fit across
    if all(isnan(trace))
        return
    end
    [t_ax, y] = adjust_for_breakpoints(trace, pre_peak_delay, t_ax, double(trace));

    %% Clear all data before peak, preceding the allowed jitter_win
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
    figure(124);cla();
    plot(x,y); hold on;
    
    thr = y(1)/20;
    if thr < 0
       thr = rms(diff(y)) * 20;
    end
    [val, secondary_peak] = findpeaks(y*-1, x, 'MinPeakProminence', thr);
    
    if ~isempty(secondary_peak)
        y = y(1:find(x == secondary_peak(1)));
        x = x(1:find(x == secondary_peak(1)));
        if rendering
            plot(x, y, 'r'); hold on % plot what's used
        end

        %% Add trace lowest point since the fit should go through that one (prevent underfitting)
        if forced_tau
            w = ones(size(y)+1);
            val = val * -1;
            [~, low] = min(val);
            x = [x; secondary_peak(low)];
            y = [y; val(low)];
            if rendering
                scatter(secondary_peak(low),  val(low), 'filled')
            end
            
            %% Give stronger weight to the end of the fit region, and any extra point    
            w(end) = numel(trace);
        else
            w = ones(size(y));
        end

                
    else
        w = ones(size(y));
        if rendering
            plot(x, y, 'g'); hold on % plot what's used
        end
    end
    
    %% For function output
    x_used = x;
    y_used = y;
    
    if ~isempty(extra_w) && forced_tau
        x = [x; extra_w(:,1)];
        y = [y; extra_w(:,2)];
        w(end+1:end+size(extra_w, 1)) = numel(w);
        %robust = 'Bisquare';
    else
        %robust = 'Off';
    end
    
    %% Adjust contraints for t0 and x0 (decay baseline)
    t0 = x(1);
    x0 = 0; 
    peak = nanmax(y)*4; % need to be more than peak, because we fit the decay part only, so y at x0 is higher than peak
    if peak < 0
        return
    end
    robust = 'Off';
    
    %% Now fit from scratch or use params and scale (if there are enough points)
    if (strcmp(fit_type, 'exp1') || numel(x) < 4) && numel(x) >= 2
        if strcmp(fit_type, 'exp2') && ~isempty(fit_params)
            w1 = fit_params(2)/(fit_params(2)+fit_params(4));
            fit_params = [fit_params(1), fit_params(2)+fit_params(4), fit_params(3)*w1 +fit_params(5)*(1-w1)];
        end
        
        if isempty(fit_params) % free robust fit
            f = fit((x-t0), y,  'exp1'                      ,...
                                'Lower'     ,[0,-Inf]       ,...
                                'Upper'     ,[peak,-1e-3]   ,...
                                'Weights'   ,w              ,...
                                'Robust'    ,'Off'     ); %
        else % guided fit
            if ~forced_tau
                fit_params_low      = [0,-10/peak_av_width,0]; % (-1/peak_av_width) width is a reasonable starting approximation for tau (i.e 67% decay)
                fit_params_high     = [peak,-0.1/peak_av_width,peak];
            else
                fit_params_low      = [peak-peak/5,-fit_params(3)]; % peak estimate can vary by 20%
                fit_params_high     = [peak+peak/5,-fit_params(3)]; % [peak+peak/5,-fit_params(3)+fit_params(3)/5]
            end
            f = fit((x-t0),y,   'exp1'                           ,...
                                'Lower'     ,single(fit_params_low(1:2)) ,...
                                'Upper'     ,single(fit_params_high(1:2)),...
                                'Weights'   ,w                   ,...
                                'StartPoint',fit_params(2:3)   ,...
                                'Robust'    ,robust);
        end
        
        if strcmp(fit_type, 'exp1')
            p = [x0, f.a, f.b*-1];
        elseif strcmp(fit_type, 'exp2')
            p = [x0, f.a, f.b*-1, 0, 0];
        end
    elseif strcmp(fit_type, 'exp2') && numel(x) >= 2
        MIN_TAU = 0.01;
        if isempty(fit_params) % free robust fit
            f = fit((x-t0), y , 'exp2'                              ,...
                                'Lower'     ,[0,-Inf,0,-Inf]        ,...
                                'Upper'     ,[peak,-MIN_TAU,peak,-MIN_TAU],...
                                'Weights'   ,w                      ,...
                                'Robust'    ,'Off'             );
        else % guided fit
            if ~forced_tau
                fit_params_low      = [0    ,-10 /peak_av_width,0   ,-10 /peak_av_width]; % (-1/peak_av_width) width is a reasonable starting approximation for tau (i.e 67% decay)
                fit_params_high     = [peak ,-1  /peak_av_width,peak,-1  /peak_av_width];
                func                = 'exp2';
                
            else % exp2 with fixed tau and constant amplitude ratio
                R = fit_params(2)/ (fit_params(2) + fit_params(4));
                func                = @(a, b, c, x) a*R*exp( b*x ) + a*(1-R)*exp( c*x );
                fit_params_low      = [0    ,-fit_params(3),-fit_params(5)]; 
                fit_params_high     = [peak ,-fit_params(3),-fit_params(5)];
                fit_params(4)       = [];                
            end
            try
            f = fit((x-t0),y   , func                             ,...
                                'Lower'     ,single(fit_params_low)         ,...
                                'Upper'     ,single(fit_params_high)        ,...
                                'Weights'   ,w                      ,...
                                'StartPoint',fit_params(2:end)      ,...
                                'Robust'    ,robust                 );
            catch
                1
            end
            if forced_tau 
                temp = {};
                temp.a = f.a*R;
                temp.b = f.b;
                temp.c = f.a*(1-R);
                temp.d = f.c;
                f = temp;
            end                            
        end
        
        %% Make sure we have fast tau first, slow tau second
        if f.b > f.d
            p = [x0, f.c, f.d*-1, f.a, f.b*-1];
        else
            p = [x0, f.a, f.b*-1, f.c, f.d*-1];
        end 
    else
        return
    end

    %% Generate fit curve    
    x_fit = (x(1):t_ax(end))';
    if numel(p) == 5
        y_fit = p(1)+p(2)*exp(-((x_fit - t0)*p(3))) + p(4)*exp(-((x_fit - t0)*p(5)));
    elseif numel(p) == 3
        y_fit = p(1)+p(2)*exp(-((x_fit - t0)*p(3)));
    end
    
    if rendering
        plot(x_fit, y_fit, 'k--');
    end
    % hold on;plot(x_fit, y_fit)
    % hold on;plot(p(2)*exp(-((x_fit - t0)*p(3))), 'k'); hold on; plot(p(4)*exp(-((x_fit - t0)*p(5))),'r'); % to plot tau 1 and tau 2
end

function p = default_p(fit_type)
    if strcmp(fit_type, 'exp2')
        p       = NaN(1, 5);
    elseif strcmp(fit_type, 'exp1')
        p       = NaN(1, 3);
    end
end

function [t_ax, y] = adjust_for_breakpoints(trace, pre_peak_delay, t_ax, y)
    %% Find breakpoints preceding the event
    break_point = find(isnan(trace(1:pre_peak_delay)));
    if any(break_point)
        y(1:break_point(end)) = -Inf;  % temporary inf range, later replaced by NaN
    end
    
    %% Find breakpoints following the eent
    break_point = find(isnan(y));
    if any(break_point)
        t_ax = t_ax(1:break_point(1)-1);
        y    = y(1:break_point(1)-1);
    end
    
    %% Restore inf to NAN.
    % We don't want to clip this region as this would shift the
    % coordinates. We just hide the data
    y(isinf(y)) = NaN;
end