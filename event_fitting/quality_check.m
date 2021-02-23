function [failed, error_estimate, half_fit_decay] = quality_check(trace, x_fit, y_fit, other_peaks)
    failed          = 0;
    error_estimate  = [];
    half_fit_decay      = find(y_fit < y_fit(1)/2, 1, 'first');
    
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
    thr     = -nanmax(y_fit)/20; % look for disprencies > 5% peak
    error_estimate(error_estimate > thr) = NaN;

    if ~isempty(other_peaks)
        other_peaks = [];
        next_peak = [];%other_peaks(find(other_peaks > x_fit(1), 1, 'first'));
    end
    
    %% Check if fit cut through signal (wrong peak estimate, need offset)
    if sum(trace(x_fit(1):x_fit(half_fit_decay)) - y_fit(1:half_fit_decay)) > (y_fit(1)/2) % see figure();plot(trace(x_fit(1):x_fit(half_decay)));hold on;plot(y_fit(1:half_decay))
        %% Walk along the slope until we reach the orginal peak height if possible
        last = min((x_fit(1) + half_fit_decay), x_fit(end));
        better = find(trace(x_fit(1):last) > y_fit(1), 1, 'last');
        if ~isempty(better)
            half_fit_decay = better;
        end
        if isempty(other_peaks) || last < next_peak
            failed  = [failed, 2];
        end
    end
    
    %% Check if at least 5 bad conscutive points ("underfit, curve to shallow")
    win = 5;
    if any(~isnan(movmean(error_estimate, win,'Endpoints', 'discard'))) 
        failed  = [failed, 3];
    end
    
    if all(isnan(error_estimate)) && any(failed)
        failed = [failed, 4];
        return
    end    
end