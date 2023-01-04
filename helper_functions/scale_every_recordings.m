%% all_traces_per_rec
% 1 x N trials cell array of t x M ROI matrices to rescale
%% demo
% if 1 ...
% if 2 ...
% pk_locs
%% If provided, rescaling is done using  peak location
% bsl_range
%% If provided, offt is computed on these datapoints only
% smoothing
%% If ~0 , presmooth traces before scaling 
% weighted
% If true, offset and peaks estimate will be weighted based on the number
% of datapoints / events (respectively)


function [global_scaling, global_offset, best_ind_scal, best_ind_offset, N_pks, N_pt_bsl] = scale_every_recordings(all_traces_per_rec, demo, pk_locs, bsl_range, smoothing, weighted)
    if nargin < 2 || isempty(demo)
        demo = 0;
    end
    if nargin < 3 || isempty(pk_locs)
        %error('to check and revise')
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
    if nargin < 4 || isempty(bsl_range)
        error('to do')    
        %bsl_percentile          = sum(tp(1:end-1)) / numel(bsl_range);
        bsl_percentile = '?';
    else
        bsl_percentile = 50; % since we alrady use a range for the baseline, we just stick to the median value
    end
    if nargin < 5 || isempty(smoothing) || ~any(smoothing)
        % no smoothing, pass
    else
        all_traces_per_rec                  = cellfun(@(x) smoothdata(x,'gaussian',smoothing), all_traces_per_rec, 'UniformOutput', false);
        %% CONSIDER --> %[~, t_peak_all]         = findpeaks(wavelet_denoise(nanmedian(all_trace,2))', 'MinPeakProminence', obj.event.peak_thr);
    end
    if nargin < 6 || isempty(weighted)
        weighted = true;
    end
    
    %% Defines how the reference trace is obtained from individual recordings.
    consensus_trace_func    = @nanmedian;
    
    %% Optimization function values (N iteration etc...)
    options                 = optimset;
    
    %% For each recording, get the best scaling factor per bin matching the cell-wide median
    best_ind_scal           = {};
    best_ind_offset         = {};
    ncores                  = (feature('numcores') * double(~(demo >= 2)));
    tp                      = [cellfun(@(x) size(x, 1), all_traces_per_rec), inf];
    N_pks                   = zeros(size(all_traces_per_rec))'; % for weights, if you use them
    N_pt_bsl                = zeros(size(all_traces_per_rec))'; % for weights, if you use them

    
    parfor (rec = 1:numel(all_traces_per_rec), ncores)
    %for rec = 1:numel(all_traces_per_rec)
        all_traces_in_rec       = all_traces_per_rec{rec};
        representative_trace    = consensus_trace_func(all_traces_in_rec, 2);
        figure(123);cla();plot(representative_trace);

        
        tp_range        = [(sum(tp(1:rec-1))+1),sum(tp(1:rec))];
        pk_tp_current   = pk_locs(pk_locs > tp_range(1) & pk_locs < tp_range(2)) - tp_range(1)+1;
        bsl_tp_current  = bsl_range(bsl_range > tp_range(1) & bsl_range < tp_range(2)) - tp_range(1)+1;        
        N_pks(rec)      = numel(pk_tp_current);
        N_pt_bsl(rec)   = numel(bsl_tp_current);

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
                best_ind_scal{rec}(trace_idx)   = NaN;
            else
                subset_for_demo = trace_idx; % can change that to a fixed index to display only one trace                
                
                %demo = 2 * double(rec == 2); 
                 
                %% Adjust offset so both traces baselines are at 0     
                if isempty(bsl_tp_current)
                    bsl_tp_current = prctile(representative_trace, 10);
                    bsl_tp_current = find(representative_trace < bsl_tp_current);
                end
                best_ind_offset{rec}(trace_idx)   = fminbnd(@(f) find_traces_offset_func(f, representative_trace(bsl_tp_current), all_traces_in_rec(bsl_tp_current,trace_idx), demo == 2 && trace_idx == subset_for_demo, bsl_percentile), 0.01, 99.99, options); % find best percentile to normalize traces
                [~ ,offset_median, current_trace] = find_traces_offset_func(best_ind_offset{rec}(trace_idx), representative_trace, all_traces_in_rec(:,trace_idx)); % apply value

                if ~isempty(pk_tp_current)
                    try
                        best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, offset_median, current_trace, demo >= 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                    catch % very rare ;  not sure why
                       % peak_times_in_record = [peak_times_in_record;peak_times_in_record]; % otherwise following function sem to crash
                        best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, offset_median, current_trace, demo >= 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                    end
                    if demo >= 3
                        uiwait(figure(6663))
                    end
                else                    
                    best_ind_scal{rec}(trace_idx) = NaN;
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
    if ~weighted
        global_scaling   = nanmedian(vertcat(best_ind_scal{:}), 1);
        global_offset    = nanmedian(vertcat(best_ind_offset{:}), 1);
        N_pks(:)         = NaN;
        N_pt_bsl(:)      = NaN;
    else
        global_scaling = w_mean(vertcat(best_ind_scal{:}), N_pks);
        global_offset   = w_mean(vertcat(best_ind_offset{:}), N_pt_bsl); 
    end    
end

function out = w_mean(var, weights)
    dim = 1;
    out = nansum(weights.*var,dim)./nansum(weights,dim);
    out(out == 0) = NaN;
end