function [all_peak_values, correction] = correct_traces(all_data, events, fits, all_peak_values, peak_times, timescale)
    
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
% %     

    N_max = numel(timescale);
    for peak = 1:size(fits, 1)
        for gp = 1:size(fits, 2)
            invalid = fits{peak, gp}(:, 1) > N_max;
            fits{peak, gp}(invalid, :) = []; 
        end
    end


    figure(1005);cla();plot(timescale, nanmedian(all_data, 2),'k'); hold on; title('median fit');set(gcf,'Color','w');xlabel('frames')
    all_mean_event_decays = cell(1,size(fits, 1));
    for peak = 1:size(fits, 1)
        events              = fits(peak, valid_gp);
        R                   = min(cellfun(@(x) nanmin(x(:,1)), events)):max(cellfun(@(x) nanmax(x(:,1)), events));
        mean_event_decay    = NaN(size(fits, 2), numel(R));
        for gp = valid_gp
            mean_event_decay(gp, fits{peak, gp}(:,1) - R(1) + 1) =  fits{peak, gp}(:,2);
        end
        all_mean_event_decays{peak} = [R', nanmean(mean_event_decay)'];
%         mean_event_decay(:,R > numel(timescale)) = [];
%         R(:,R > numel(timescale)) = [];
        
        plot(timescale(R'), nanmean(mean_event_decay)','r--'); hold on
    end
    
    %% Now correct
    for gp = valid_gp
        figure(10050 + gp);cla();plot(timescale, all_data(:,gp),'k'); hold on; title(['corrected amplitudes for group ',num2str(gp)]);set(gcf,'Color','w');xlabel('frames')
        for peak = 1:size(fits, 1)
            R = all_mean_event_decays{peak}(:,1);
            
%             all_mean_event_decays{peak}(R > numel(timescale),:) = [];
%             fits{peak, gp}(R > numel(timescale),:) = [];
%             R(R > numel(timescale)) = [];
%             
            
            dec = NaN(1, numel(R));
            dec(fits{peak, gp}(:,1) - R(1) + 1) =  fits{peak, gp}(:,2);
            plot(timescale(R), dec','r--'); hold on
            
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
                    to_correct = to_correct(to_correct ~= peak); % filter the peak itself
                    for pk = to_correct
                        t_next_peak_within_window = find(t == peak_times(pk));
                        residual = y(t_next_peak_within_window);
                        if isnan(residual)
                            residual = y(find(t == peak_times(pk))-1);
                        end
                        quiver(timescale(peak_times(pk)),all_peak_values(pk, gp),0,-residual,'b'); hold on
                        all_peak_values(pk, gp) = all_peak_values(pk, gp) - residual;                        
                    end                  
                end
            end
            correction(t(1:overlap), gp) = y(1:overlap);   
            correction(t(overlap),   gp) = NaN         ;
        end
    end
end