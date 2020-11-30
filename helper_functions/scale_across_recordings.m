function [all_traces_per_bin, all_traces] = scale_across_recordings(all_traces_per_bin, timescale, bin_legend, demo)

    options = optimset;%('TolX',20); % Prevents oversubstraction of baseline
    
    func = @nanmedian;
    
    %% For each recording, get the best scaling factor per bin matching the cell-wide median
    ideal_scal_f = {};
    ideal_offset = {};
    all_traces   = all_traces_per_bin;

    %% For testing scaling criteria
%     estimates = NaN(numel(all_traces_per_bin), size(all_traces_per_bin{1}, 2));
%     for gp = 1:size(all_traces_per_bin{1}, 2)
%         test = cellfun(@(x) x{gp} ,all_traces_per_bin, 'UniformOutput', false); % 'end' or any group you want to test
%     	test = cellfun(@(x) nanmedian(x, 2) ,test, 'UniformOutput', false);
%         try
%             test = cellfun(@(x) wavelet_denoise(x) ,test, 'UniformOutput', false);
%             estimates(:, gp) = cellfun(@(x) iqr(x) ,test);
%         end
%     end
%     estimates = nanmedian(estimates) * 5;
%     % figure();val = cellfun(@(x) iqr(x) ,test);plot(val); val

    %% Autoestimate peak detection level
    all_concat = vertcat(all_traces_per_bin{:});
    all_concat = cell2mat(cellfun(@(x) nanmedian(x, 2)', all_concat, 'UniformOutput', false)');
    findpeaks(wavelet_denoise(nanmedian(all_concat, 1)')', 'MinPeakProminence', rms(nanmedian(all_concat, 1)'));
    [pks, locs]= findpeaks(wavelet_denoise(nanmedian(all_concat, 1)')', 'MinPeakProminence', rms(nanmedian(all_concat, 1)'));
    thr_estimates  = prctile(all_concat(:, locs), 10, 2);
    thr_estimates(thr_estimates < 0) = min(thr_estimates(thr_estimates > 0)); % no negative values

    ncores          = (feature('numcores') * double(~(demo == 2)));
    %parfor (data_f_idx = 1:numel(all_traces_per_bin), ncores)
    for data_f_idx = 1:numel(all_traces_per_bin)
        current_rec_per_bin = all_traces_per_bin{data_f_idx};
        
        %% Fix a rare case where some bins are empty.
        current_rec_per_bin = fill_missing(current_rec_per_bin);
              
        median_ref_per_bin  = cell2mat(cellfun(@(x) func(x, 2) , current_rec_per_bin, 'UniformOutput', false));
        median_ref_per_bin  = wavelet_denoise(median_ref_per_bin);
        median_trace        = func(median_ref_per_bin, 2);
        
        %% try to get a good estimate of the noise level for peak-scaling
        %peak_estimate = nanmedian(iqr(median_trace)) * 3;
        %peak_estimate = max(peak_estimate, mad(median_trace)*2);
        %[peak_estimate, rms(median_trace), mad(median_trace), mode(median_trace)]
        
        % median_trace        = detrend(median_trace);
        ideal_scal_f{data_f_idx}    = [];
        ideal_offset{data_f_idx}    = [];
        for bin_idx = 1:size(median_ref_per_bin, 2)
            if any(~isnan(median_ref_per_bin(:,bin_idx)))
                %% Adjust offset so both traces baselines are at 0
                ideal_offset{data_f_idx}(bin_idx)   = fminbnd(@(f) find_traces_offset_func(f, median_trace, median_ref_per_bin(:,bin_idx), demo == 2), 0.1, 99.9, options); % find best percentile to normalize traces
                [~ ,modified_median, current_trace] = find_traces_offset_func(ideal_offset{data_f_idx}(bin_idx), median_trace, median_ref_per_bin(:,bin_idx)); % apply value
                
                %% Get the best scaling factor between the 2 curves
                ideal_scal_f{data_f_idx}(bin_idx) = fminbnd(@(f) scale_trace_func(f, modified_median, current_trace, demo == 2, thr_estimates(bin_idx)), 1e-3, 20, options); % scaling factor must be > 0
                [~,~,~,n_peaks] = scale_trace_func(ideal_scal_f{data_f_idx}(bin_idx), modified_median, current_trace, demo == 2, thr_estimates(bin_idx));
%                 if ~n_peaks
%                     ideal_scal_f{data_f_idx}(bin_idx) = NaN;
%                 end
            else
                ideal_offset{data_f_idx}(bin_idx) = NaN;
                ideal_scal_f{data_f_idx}(bin_idx) = NaN;
            end
        end
    end
    clear data_f_idx % just to remove parfor warning

    %% All recordings, compute the best scaling factor per bin matching the overall median
    best_isf        = func(vertcat(ideal_scal_f{:}), 1);
    best_ioffset    = func(vertcat(ideal_offset{:}), 1);
    
    n_gp = numel(ideal_scal_f{1});
    n_rec = numel(ideal_scal_f);
    figure(1012);cla();hold on;
    col = mat2cell(viridis(n_gp), ones(1,n_gp), 3);
    p = plot(1:n_rec, vertcat(ideal_scal_f{:}),'o-');set(gcf,'Color','w');
    arrayfun(@(x, y) set(x, 'Color', y{1}), p, col);
    legend(bin_legend);
    xticks(1:n_rec)
    xlim([0.8,n_rec+0.2])
    ylabel('intra-group Scaling Factor');
    xlabel('Recording #');
    title('Scaling factor per group, per trial');

    figure(1013);clf();hold on;set(gcf,'Color','w');
    xlim([1,n_gp]);hold on;
    title('global scaling factor and offset per group');
    xlabel('Groups');
    xticklabels(bin_legend);
    xtickangle(45);
    xticks(1:n_gp)
    plot(best_isf,'ko-'); hold on;ylabel('Consensus Scaling Factor');    
    yyaxis right;plot(best_ioffset,'ro-');ylabel('Consensus Baseline Percentile'); 
    xlim([0.8,n_gp+0.2])
    ax = gca;
    ax.YAxis(1).Color = 'r';
    ax.YAxis(2).Color = 'k';
    

    %% Now rescale each bin with a single value for all recordings.
    % so for a given bin, it's the same scaling factor for all
    % recordings
    for gp = 1:n_gp
        %% Do median then scale ( tried to do indivual trace scaling then median, but it doesnt work because that's not how th scaling factor was computed
        test = cellfun(@(x) func(x{gp},2) ,all_traces, 'UniformOutput', false); % 'end' or any group you want to test
        test = cellfun(@(x) x - prctile(x, best_ioffset(gp)), test, 'UniformOutput', false);
        test = cellfun(@(x) x / best_isf(gp), test, 'UniformOutput', false);
        for data_f_idx = 1:numel(all_traces_per_bin)
            all_traces_per_bin{data_f_idx}{gp} = test{data_f_idx};
            all_traces_per_bin{data_f_idx}{gp}(1:2,:) = NaN; % first 2 points are sometimes a bit lower than the rest
        end
        
        %% Now do the same on individual traces
        test = cellfun(@(x) x{gp} ,all_traces, 'UniformOutput', false); % 'end' or any group you want to test
        test = cellfun(@(x) x - prctile(x, best_ioffset(gp)), test, 'UniformOutput', false);
        test = cellfun(@(x) x / best_isf(gp), test, 'UniformOutput', false);
        for data_f_idx = 1:numel(all_traces_per_bin)
            all_traces{data_f_idx}{gp} = test{data_f_idx};
            all_traces{data_f_idx}{gp}(1:2,:) = NaN; % first 2 points are sometimes a bit lower than the rest
        end
    end
    all_traces_per_bin = cellfun(@(x) cell2mat(x) ,all_traces_per_bin, 'UniformOutput', false);
end

function current_rec_per_bin = fill_missing(current_rec_per_bin)
    pb = cellfun(@isempty, current_rec_per_bin);
    if any(pb)
        ref = size(current_rec_per_bin{find(~pb, 1, 'first')});
        current_rec_per_bin(pb) = repmat({single(NaN(ref(1),1))}, 1, sum(pb));
    end
end

