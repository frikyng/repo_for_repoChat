function [best_isf, best_ioffset] = scale_every_recordings(all_traces_per_rec, demo)
    if nargin < 2 || isempty(demo)
        demo = 0;
    end
    func1                   = @nanmedian;
    func2                   = @nanmax;
    options                 = optimset;
    
    %% For each recording, get the best scaling factor per bin matching the cell-wide median
    ideal_scal_f            = {};
    ideal_offset            = {};
       
    %% Autoestimate peak detection level
    all_traces_concat       = vertcat(all_traces_per_rec{:});
    representative_trace    = nanmedian(all_traces_concat, 2)';
    
    %% Filter signal
    missing                 = isnan(representative_trace);
    representative_trace    = fillmissing(representative_trace,'linear');
    denoised                = wdenoise(double(representative_trace));
    
    noise                   = representative_trace - denoised;
    
    %% For threshold debugging
%         figure(25);cla();plot(representative_trace)
%         hold on;plot(denoised,'r')
%         hold on;plot(noise,'b')
%         hold on; envelope(noise,50,'rms');
    
    [top_noise, ~]          = envelope(noise,50,'rms'); % used to be 20, 'peak'
    initial_thr             = mean(top_noise)*5;    
    %initial_thr             = (prctile(top_noise-bottom_noise,99) - prctile(top_noise-bottom_noise,1))*5;
    [pks, locs]             = findpeaks(wavelet_denoise(representative_trace')', 'MinPeakProminence', initial_thr);
    
    pk_range                = arrayfun(@(x) (x-20):(x+20), locs', 'uni', false);
    pk_range                = unique([horzcat(pk_range{:})]);
    pk_range                = pk_range(pk_range > 0 & pk_range <= numel(denoised));
    denoised(~ismember(1:numel(denoised), pk_range)) = NaN;
    
    ncores                  = (feature('numcores') * double(~(demo == 2)));
    tp                      = [cellfun(@(x) size(x, 1), all_traces_per_rec), inf];
    
    parfor (rec = 1:numel(all_traces_per_rec), ncores)
    %for rec = 1:numel(all_traces_per_rec)
        all_traces_in_rec       = movmean(all_traces_per_rec{rec}, 3); % was 50
        representative_trace    = func1(all_traces_in_rec, 2);
        figure(123);cla();plot(representative_trace)

        ideal_scal_f{rec}    = [];
        ideal_offset{rec}    = [];
        for trace_idx =  1:size(all_traces_in_rec, 2) 
            if all(isnan(all_traces_in_rec(:,trace_idx)))
                ideal_offset{rec}(trace_idx) = NaN;
                ideal_scal_f{rec}(trace_idx) = NaN;
            else
                subset_for_demo = trace_idx; % can change that to a fixed index to display only one trace
                
                %% Adjust offset so both traces baselines are at 0            
                ideal_offset{rec}(trace_idx)   = fminbnd(@(f) find_traces_offset_func(f, representative_trace, all_traces_in_rec(:,trace_idx), demo == 2 && trace_idx == subset_for_demo), 0.01, 99.99, options); % find best percentile to normalize traces
                [~ ,offset_median, current_trace] = find_traces_offset_func(ideal_offset{rec}(trace_idx), representative_trace, all_traces_in_rec(:,trace_idx)); % apply value

                %% Get the best scaling factor between the 2 curves
                if rec > 1
                    peak_times_in_record = locs - sum(tp(1:rec-1));
                else
                    peak_times_in_record = locs;
                end
                peak_times_in_record = peak_times_in_record(peak_times_in_record > 0 & peak_times_in_record < (tp(rec)+1));
                ideal_scal_f{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, offset_median, current_trace, demo == 2 && trace_idx == subset_for_demo, peak_times_in_record), 1e-3, 100, options); % scaling factor must be > 0
            end
        end
    end
    clear data_f_idx % just to remove parfor warning

    %% All recordings, compute the best scaling factor per bin matching the overall median
    best_isf        = func1(vertcat(ideal_scal_f{:}), 1);
    %best_ioffset    = func1(vertcat(ideal_offset{:}), [], 1);
    best_ioffset    = func1(vertcat(ideal_offset{:}), 1);

    n_gp = numel(ideal_scal_f{1});
    n_rec = numel(ideal_scal_f);
    figure(1012);clf();hold on;
    col = mat2cell(viridis(n_gp), ones(1,n_gp), 3);
    p = plot(1:n_rec, vertcat(ideal_scal_f{:}),'o-');set(gcf,'Color','w');
    arrayfun(@(x, y) set(x, 'Color', y{1}), p, col);
    set(gca, 'YScale', 'log')
    %legend(bin_legend);
    xticks(1:n_rec)
    xlim([0.8,n_rec+0.2])
    ylabel('Scaling Factor');
    xlabel('Recording #');
    title('Scaling factor per ROI, per trial');

    figure(1013);clf();hold on;set(gcf,'Color','w');
    xlim([1,n_gp]);hold on;
    title('global scaling factor and offset per ROI');
    xlabel('ROI');
    plot(best_isf,'ko-'); hold on;ylabel('Consensus Scaling Factor');    
    yyaxis right;plot(best_ioffset,'ro-');ylabel('Consensus Baseline Percentile'); 
    ax = gca;
    ax.YAxis(1).Color = 'r';
    ax.YAxis(2).Color = 'k';
    
    %all_traces_concat_scaled = vertcat(all_traces_per_rec{:});
end