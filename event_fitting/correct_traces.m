function [all_peak_values, correction] = correct_traces(all_data, events, fits, all_peak_values, peak_times)
    
    if isempty(events)
        [all_peak_values, correction] = deal([]);
        return
    end
    
    
    %% Finally, correct the events amplitudes
    correction = NaN(size(all_data));
    valid_gp   = find(~all(isnan(all_data))); % QQ not sure when we have full NaN series, but if any, the group will be empty so we skip ths part of the loop
    
    
%     %% Plot all trace and decay fits
%     for gp = 1:size(fits, 2)
%         figure();plot(all_data(:,gp),'k'); hold on;
%         for el = 1:size(fits, 1)
%             plot(fits{el,gp}(:,1),fits{el,gp}(:,2),'r--'); hold on;
%         end
%     end
%     

    figure(6664);cla();plot(nanmedian(all_data, 2),'k'); hold on; title('median fit');set(gcf,'Color','w');xlabel('frames')
    all_mean_events = {};
    for peak = 1:size(fits, 1)
        events = fits(peak, :);
        R = min(cellfun(@(x) nanmin(x(:,1)), events)):max(cellfun(@(x) nanmax(x(:,1)), events));
        mean_event = NaN(size(fits, 2), numel(R));
        for gp = 1:size(fits, 2)
            mean_event(gp, fits{peak, gp}(:,1) - R(1) + 1) =  fits{peak, gp}(:,2);
        end
        all_mean_events{peak} = [R', nanmean(mean_event)'];
        plot(R', nanmean(mean_event)','r--'); hold on
    end
    
    %% Now correct
    for peak = 1:size(events, 3)
        for gp = valid_gp
            t = fits{peak,gp}(:, 1);
            y = fits{peak,gp}(:, 2);
            if any(t > size(correction, 1))
                lim = find(t == size(correction, 1));
                t = t(1:lim);
                y = y(1:lim);
            elseif any(t < 1)
                y(t < 1) = [];
                t(t < 1) = [];
            end
            overlap = find(~isnan(correction(t,gp)), 1, 'first');
            if isempty(overlap)
                %% no peak within the fitting window
                overlap = size(t, 1);
            else
                %% One or more peak to correct
                to_correct = find(ismember(peak_times, t));
                if ~isempty(to_correct) && sum(to_correct ~= peak)
                    to_correct = to_correct(to_correct ~= peak); % fitler the peak itself
                    for pk = to_correct
                        residual = y(t == peak_times(pk));
                        if isnan(residual)
                            residual = y(find(t == peak_times(pk))-1);
                        end
                        all_peak_values(pk, gp) = all_peak_values(pk, gp) - residual;
                        
                    end                  
                end
            end
            correction(t(1:overlap), gp) = y(1:overlap);   
            correction(t(overlap),   gp) = NaN         ;
        end
    end
end