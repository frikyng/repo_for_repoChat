function [global_scal, global_offset, best_ind_scal, best_ind_offset] = scale_every_recordings(all_traces_per_rec, demo, pk_locs, bsl_range, smoothing)
    if nargin < 2 || isempty(demo)
        demo = 0;
    end
    func1                   = @nanmedian;
    func2                   = @nanmean;
    options                 = optimset;
    
    %% For each recording, get the best scaling factor per bin matching the cell-wide median
    best_ind_scal           = {};
    best_ind_offset         = {};

   % all_traces_concat       = vertcat(all_traces_per_rec{:});
%     ref_signal              = nanmedian(all_traces_concat, 2)';
    
    if nargin < 3 && ~isempty(pk_locs)
        %% Autoestimate peak detection level

        %% Filter signal
        representative_trace    = fillmissing(representative_trace,'linear');
        ref_signal              = wdenoise(double(representative_trace));

        noise                   = representative_trace - ref_signal;

        %% For threshold debugging
    %         figure(25);cla();plot(representative_trace)
    %         hold on;plot(denoised,'r')
    %         hold on;plot(noise,'b')
    %         hold on; envelope(noise,50,'rms');

        [top_noise, ~]          = envelope(noise,50,'rms'); % used to be 20, 'peak'
        initial_thr             = mean(top_noise)*5;    
        %initial_thr             = (prctile(top_noise-bottom_noise,99) - prctile(top_noise-bottom_noise,1))*5;
        [pks, pk_locs]             = findpeaks(wavelet_denoise(representative_trace')', 'MinPeakProminence', initial_thr);
    end
    
%     if nargin < 4 || isempty(pk_range)
%         pk_range                = arrayfun(@(x) (x-20):(x+20), locs', 'uni', false);
%         pk_range                = unique([horzcat(pk_range{:})]);
%         pk_range                = pk_range(pk_range > 0 & pk_range <= numel(ref_signal));
%     end
%     ref_signal(~ismember(1:numel(ref_signal), pk_range)) = NaN;
%     


    if nargin < 5 || isempty(smoothing)
        % no smoothing
    else
        all_traces_per_rec                  = cellfun(@(x) smoothdata(x,'gaussian',smoothing), all_traces_per_rec, 'UniformOutput', false);
                        %[~, t_peak_all]         = findpeaks(wavelet_denoise(nanmedian(all_trace,2))', 'MinPeakProminence', obj.event.peak_thr);

    end


    ncores                  = (feature('numcores') * double(~(demo == 2)));
    tp                      = [cellfun(@(x) size(x, 1), all_traces_per_rec), inf];
    
    %bsl_percentile          = sum(tp(1:end-1)) / numel(bsl_range);
    bsl_percentile = 50;% since we alrady use a range
    parfor (rec = 1:numel(all_traces_per_rec), ncores)
    %for rec = 1:numel(all_traces_per_rec)
        all_traces_in_rec       = all_traces_per_rec{rec};
        representative_trace    = func1(all_traces_in_rec, 2);
        figure(123);cla();plot(representative_trace);

        tp_range = [(sum(tp(1:rec-1))+1),sum(tp(1:rec))];
        pk_tp_current = pk_locs(pk_locs > tp_range(1) & pk_locs < tp_range(2)) - tp_range(1)+1;
        bsl_tp_current = bsl_range(bsl_range > tp_range(1) & bsl_range < tp_range(2)) - tp_range(1)+1;
%         
%         
%         tp_range                      = 1:size(all_traces_in_rec, 1);
%         tp_range(bsl_range(bsl_range < 1490)) = NaN;
%         pk_tp = isnan(tp_range);
%         bsl_tp = find(~pk_tp);
%         pk_tp = find(pk_tp);
        

        best_ind_scal{rec}    = [];
        best_ind_offset{rec}    = [];
        for trace_idx =  1:size(all_traces_in_rec, 2) 
            if all(isnan(all_traces_in_rec(:,trace_idx)))
                best_ind_offset{rec}(trace_idx) = NaN;
                best_ind_scal{rec}(trace_idx) = NaN;
            else
                subset_for_demo = trace_idx; % can change that to a fixed index to display only one trace                
                demo = 0
                %% Adjust offset so both traces baselines are at 0            
                best_ind_offset{rec}(trace_idx)   = fminbnd(@(f) find_traces_offset_func(f, representative_trace(bsl_tp_current), all_traces_in_rec(bsl_tp_current,trace_idx), demo == 2 && trace_idx == subset_for_demo, bsl_percentile), 0.01, 99.99, options); % find best percentile to normalize traces
                [~ ,offset_median, current_trace] = find_traces_offset_func(best_ind_offset{rec}(trace_idx), representative_trace, all_traces_in_rec(:,trace_idx)); % apply value

                %% Get the best scaling factor between the 2 curves
%                 if rec > 1
%                     peak_times_in_record = pk_locs - sum(tp(1:rec-1));
%                 else
%                     peak_times_in_record = pk_locs;
%                 end
%                 peak_times_in_record = peak_times_in_record(peak_times_in_record > 0 & peak_times_in_record < (tp(rec)+1));
%                 
                demo = 0;
                try
                    best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, offset_median, current_trace, demo == 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                catch % very rare ;  not sure why
                   % peak_times_in_record = [peak_times_in_record;peak_times_in_record]; % otherwise following function sem to crash
                    best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, offset_median, current_trace, demo == 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                end
                
                %% #### DEBUG
                %% Show why we need to do scaling on events --> bsl noise scaling way too variable
                %   figure();plot(offset_median .\ current_trace)
                %% Show why we canot take the entir event for scaling --> decays perturbate measurment
                %   figure();plot(offset_median(pk_tp) .\ current_trace(pk_tp))
                
                %figure(123);cla();plot(offset_median);hold on; plot(current_trace/best_ind_scal{rec}(trace_idx))
            end
        end
    end
    clear data_f_idx % just to remove parfor warning

    %% All recordings, compute the best scaling factor per bin matching the overall median
    global_scal      = func2(vertcat(best_ind_scal{:}), 1);
    %best_ioffset    = func2(vertcat(ideal_offset{:}), [], 1);
    global_offset    = func2(vertcat(best_ind_offset{:}), 1);
end