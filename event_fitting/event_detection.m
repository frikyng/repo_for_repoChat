classdef event_detection < handle
    %% subclass of arboreal_scan_experiment 
    properties
        thr_for_detection   = 0.2;
        thr_for_global      = 0.5;
        thr_for_bad_ROIs    = 0.2;
    end

    methods
        function [bad_ROIs, mean_corr_with_others] = find_events(obj, idx_filter, thr_for_detection, method, thr_for_global)
            if nargin < 2 || isempty(idx_filter)
                idx_filter = 1:obj.n_ROIs;
            elseif ischar(idx_filter) && strcmp(idx_filter, 'soma')
                idx_filter = obj.ref.indices.somatic_ROIs;
            end
            if nargin < 3 || isempty(thr_for_detection)
                thr_for_detection = obj.thr_for_detection;
            else
                obj.thr_for_detection = thr_for_detection;
            end
            if nargin < 4 || isempty(method)
                method = 'corr';
            end
            if nargin < 5 || isempty(thr_for_global)
                thr_for_global = obj.thr_for_global;
            else
                obj.thr_for_global = thr_for_global;
            end            

            obj.disp_info(['Now detecting events with tree-wide pairwise correlation at least > ',num2str(thr_for_detection),' %%'],1)

            %% Get original traces (unscaled)
            raw_traces              = obj.extracted_traces_conc(:, idx_filter);

            %% Get original traces
            [obj.event, correlation_res] = detect_events(raw_traces, obj.t, method, thr_for_detection, [], obj.rendering,'', thr_for_global);
            
%             sz = vertcat(obj.event.peak_value{:});
%             thr1 = max(sz) * 0.8;
%             thr2 = max(sz) * 0.2;
%             mask = cellfun(@(x) mean(x) < thr1 & mean(x) > thr2, obj.event.peak_value);            
%             for fn = fieldnames(obj.event)'                
%                  if numel(obj.event.(fn{1})) == numel(mask)
%                      obj.event.(fn{1}) = obj.event.(fn{1})(mask);
%                  end
%             end

            %% Identify and log poorly correlated ROIs
            [bad_ROIs, mean_corr_with_others] = obj.find_bad_ROIs(correlation_res, idx_filter);
            obj.disp_info('Events were detected and obj.event field has been populated.',1)
        end
        
        function fit_data = fit_events(obj, tau_decay, event_win_size, tolerance)
            if nargin < 2 || isempty(tau_decay)
                tau_decay               = 0.2; % To make sure that we have enough data point to fit events (otherwise, we resample data)
            end
            if nargin < 3 || isempty(event_win_size)
                event_win_size          = [1,5]; % Expressed as # of event average widths before and after peak time
            end
            if nargin < 4 || isempty(tolerance)
                tolerance               = [0.2,0.2]; % pre and post peak toleance (in s) for jitter across regions
            end
            
            %%% #### STEP 1 - DETECTION SETTINGS AND PREPARATION - ######
            
            %% Make sure events were detected
            if isempty(obj.event)
                obj.find_events();
            end
            
            traces                  = obj.binned_data.median_traces;
            peak_pos                = unique(sort(vertcat(obj.event.peak_time{obj.event.is_global})))';
            peak_av_width           = nanmean(vertcat(obj.event.peak_width{obj.event.is_global}));
            dt                      = nanmedian(obj.timescale.sr);
            global_timescale        = obj.t;

            %% Build event fitting holder (assigned as obj.events.fit_data only at the end)
            fit_data                = {};
            fit_data.events         = {};
            fit_data.group_labels   = obj.binned_data.bin_legend;
            
            %% TBD
            refit                   = true
            fit_data.fast_fit       = true

            %% We do the fitting on denoised data - THIS MUST BE DONE BEFORE RESAMPLING
            %% ## DEBUG % figure();plot(all_data,'r'); [test_denoised, ~] = wavelet_denoise(all_data);hold on;plot(test_denoised,'k');
            [traces, ~] = wavelet_denoise(traces);
            
            %% QQ NOT GREAT
            %traces = tweak_scaling(traces, peak_pos);
            
            %% Fitting adjustements with low SR can be a bit funny when cells are bursting, but upsampling data a bit does the job            
            if (dt/2) < tau_decay
                fit_data.resampling_factor          = ceil(tau_decay / (dt/2));        
            else
                fit_data.resampling_factor          = 1;
            end 
            
            if fit_data.resampling_factor > 1
                global_timescale    = interpolate_to(global_timescale, numel(global_timescale) * fit_data.resampling_factor); 
                peak_pos            = round(peak_pos * fit_data.resampling_factor);
                peak_av_width       = round(peak_av_width*fit_data.resampling_factor);
                dt                  = dt/fit_data.resampling_factor;
                try
                    traces            = interpolate_to(traces, numel(global_timescale), 'cubic')'; % failed with 2019-11-06_exp_1  
                catch
                    traces            = interpolate_to(traces, numel(global_timescale))';                
                end
            end

            %% Get adjusted peak_times and width following resampling
            fit_data.peak_times = global_timescale(peak_pos);
            fit_data.peak_pos   = peak_pos;
            
            %% We define a reasonable fitting window around each event
            pre_peak_delay  = round(peak_av_width * event_win_size(1));
            post_peak_delay = round(peak_av_width * event_win_size(2));
            jitter_win      = round([tolerance*dt,tolerance*dt]); % 200 ms jitter around peaks for event detection
            
            %%% #### STEP 2 - EXTRACTION - ######
            
            %% Extract all events
            if isempty(fit_data.peak_pos)
                fit_data.peak_pos = NaN; % to avoid default detection
            end
            events = obj.extract_events(traces, fit_data.peak_pos, pre_peak_delay, post_peak_delay); %% QQ SEE WHAT TO DO WITH TH EXTRA INPUT "exact_t"

            %% Get global fit settings to have a rough estimate of the fitting window required. Use small, isolated events (or end-of-burst)
            MIN_SPACING = post_peak_delay;
            to_use = find(diff(fit_data.peak_pos) > MIN_SPACING); % keep events spaced by at least 5 average peak width from the next, to avoid contamination of decay
            figure(1019);cla();hold on;plot(reshape(events,size(events,1),[]), 'Color',[0.9,0.9,0.9]); hold on;title('Individual median events amd median fit');hold on;
            fit_data.global_fit = fit_decay(squeeze(nanmean(nanmedian(events(1:round(MIN_SPACING),:,to_use), 2),3)), pre_peak_delay, [], true, [], peak_av_width, [], 'exp2'); % the fitting window can be shorter as we don't have bursts in this subselection
            xlabel('points');legend({'mean event','fit result'}); hold on;
            arrangefigures([1,2]); 
            
            %% In fast_fit mode, we force tau to its global value (across regions) and don't need a refit
            if fit_data.fast_fit
                guesses     = NaN(size(events, 2), size(fit_data.global_fit,2));
                for gp = 1:size(events, 2)
                    guesses(gp, :) = fit_decay(squeeze(nanmean(events(1:round(peak_av_width * 5),gp,to_use),3)), pre_peak_delay, [], false, [], peak_av_width, [], 'exp2');
                end
                forced_tau  = true;
                refit       = false;        
            else
                guesses     = repmat(fit_data.global_fit,size(events, 2),1);
                forced_tau  = false;        
            end

            %% Now fit all events with fairly free parameters.
            % initial guesses is set by the median tau 
            fit_rendering = false
            [all_peak_values, fit_data] = fit_all_events(events, traces, fit_data.peak_pos, pre_peak_delay, fit_data, fit_rendering, guesses, peak_av_width, 'exp2', forced_tau, jitter_win);
            fit_data.pre_refit_events   = fit_data.events;
            
            %% Make some figures about fit properties
            fits_free = 1/permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
            if obj.rendering
                figure(1014);clf();plot(fits_free(:,:,3)) ; hold on; plot(nanmedian(fits_free(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');xticklabels(fit_data.group_labels);xtickangle(45);set(gcf,'Color','w');
                figure(1015);clf();plot(fits_free(:,:,3)'); hold on; plot(nanmedian(fits_free(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');
                arrangefigures([1,2]); 
            end

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
               fit_rendering = false;
                [all_peak_values, fit_data] = fit_all_events(events, traces, fit_data.peak_pos, pre_peak_delay, fit_data, fit_rendering, guesses, peak_av_width, 'exp1',true);
                refitted_events                 = 1/permute(cat(3, fit_data.events{:}),[1,3,2]); % gp, event, parameter
                figure(1014);clf();plot(refitted_events(:,:,3)) ; hold on; plot(nanmedian(refitted_events(:,:,3),2), 'k', 'LineWidth',2); title('median tau 1 per group'); xlabel('Groups');xticklabels(fit_data.group_labels);xtickangle(45);set(gcf,'Color','w');
                figure(1015);clf();plot(refitted_events(:,:,3)'); hold on; plot(nanmedian(refitted_events(:,:,3),1), 'k', 'LineWidth',2); title('median tau 1 per event'); xlabel('Events');set(gcf,'Color','w');  
                fit_data.pre_correction_peaks   = all_peak_values;
            end

            %% Correct the traces before downsampling
            %[all_peak_values, correction]   = correct_traces(all_data, events, fit_data.window, all_peak_values, fit_data.peak_pos, global_timescale);
            fit_data.post_correction_peaks      = all_peak_values;

            %% Revert the effect of fit_data.resampling_factor (so timepoints match the real trace)
            if fit_data.resampling_factor ~= 1
                traces = obj.binned_data.median_traces;
                fit_data = obj.fix_resampling(fit_data);
            end

            %% Finally, store and show results  
            obj.event.fitting = fit_data;
            if obj.rendering
                obj.show_results(fit_data, traces, fit_data.group_labels, obj.t);arrangefigures([1,2]);
            end
        end

        function fit_data = fix_resampling(~, fit_data)
            %% If the initial trace was resampled, some elements needs to be readjusted
            fit_data.peak_pos = round(fit_data.peak_pos / fit_data.resampling_factor);
            for el = 1:numel(fit_data.window)
                fit_data.window{el}         = interpolate_to(fit_data.window{el}, ceil(size(fit_data.window{el},1)/fit_data.resampling_factor))';
                fit_data.window{el}(:,1)    = round(fit_data.window{el}(:,1)/fit_data.resampling_factor);
            end
            % .peak_times is ok
            % .peak_jitter TO FIX
            % .global_fit & events TO CHECK --> possibly tau to fix
        end
        
        function events = extract_events(obj, data, peak_loc, pre_peak_delay, post_peak_delay, exact_t_event)
            %% WARNING : if used in standalone mode, there is no data resampling
            if nargin < 2 || isempty(data) 
                data                    = obj.binned_data.median_traces;
            end
            if nargin < 3 || isempty(peak_loc)
                obj.events              = obj.find_events();
                peak_loc                = unique(sort(vertcat(obj.event.peak_time{obj.event.is_global})))';
            end
            if nargin < 4 || isempty(pre_peak_delay) 
                peak_av_width           = nanmean(vertcat(obj.event.peak_width{obj.event.is_global}));
                pre_peak_delay          = round(peak_av_width * event_win_size(1));
            end
            if nargin < 5 || isempty(post_peak_delay)
                peak_av_width           = nanmean(vertcat(obj.event.peak_width{obj.event.is_global}));
                post_peak_delay         = round(peak_av_width * event_win_size(2));
            end
            if nargin < 6 || isempty(exact_t_event)
                exact_t_event = [];
            end

            events = [];    

            %% In case you are in a case with no event > thr
            if isempty(peak_loc) || all(isnan(peak_loc))
                return
            end

            %% For 99% of the other cases
            for peak = 1:size(peak_loc, 2)

                pr_pd = pre_peak_delay;
                po_pd = post_peak_delay;

                %% Make sure window doesn't exceed timescale
                if (peak_loc(peak)-pr_pd) <= 1 
                    pr_pd = peak_loc(peak) - 1;
                elseif (peak_loc(peak)+po_pd) > size(data, 1)
                    po_pd = size(data, 1)- peak_loc(peak);
                end

                %% extract events
                if isempty(exact_t_event)
                    % standard case, centred around peak location of the median
                    % or mean trace
                    new = data((peak_loc(peak)-pr_pd):(peak_loc(peak)+po_pd),:);
                    if pr_pd ~= pre_peak_delay
                       new = [NaN(pre_peak_delay - pr_pd, size(new, 2)); new];
                    elseif po_pd ~= post_peak_delay
                       new = [new; NaN(post_peak_delay - po_pd, size(new, 2))];
                    end
                else
                    % custom case, obtained after fitting individual events,
                    % and adjusted for jitter
                    new = repmat((peak_loc(peak)-pr_pd):(peak_loc(peak)+po_pd), size(data, 2), 1) + exact_t_event{peak};                
                end

                events = cat(3, events, new);
            end
            
            %% DEBUG % show_events(data, events, pre_peak_delay);
        end

        function show_results(obj, fit_data, traces, group_labels, global_timescale)
            R                             = max(range(traces));
            peak_values_raw             = fit_data.pre_correction_peaks;
            peak_values_corrected       = fit_data.post_correction_peaks;
            for gp = 1:size(peak_values_corrected, 2)
                figure(10050 + gp);cla();hold on;title(['Traces and fits group ',num2str(gp),' (',fit_data.group_labels{gp},')']);set(gcf,'Color','w');
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
    end
end

