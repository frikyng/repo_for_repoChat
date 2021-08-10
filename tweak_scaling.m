function [traces, ideal_offset, ideal_gain] = tweak_scaling(traces, peak_pos, smoothing)
    traces                  = smoothdata(traces,'gaussian',smoothing);

    %% Nan'd traces will be ignored
    valid = ~all(isnan(traces),1);

    %% Set range for 
    offset_range = linspace(0,50,500);
    gain_range   = logspace(-1,2,1000);
    
    %% Traces baseline should not be negative.
    % As we do not know what the optimal gain is, we try different gain
    % with different prctile values for the baseline. The baseline estimate

    best_offset  = NaN(numel(offset_range),size(traces, 2));
    parfor offset_idx = 1:numel(offset_range)
        offset_pct = offset_range(offset_idx);
        temp_traces             = traces - prctile(traces, offset_pct);        
        best_offset(offset_idx,:) = nanmin(temp_traces);
    end
    
    best_offset = smoothdata(diff(diff(best_offset)),1,'gaussian',size(best_offset,1)/20);
    ideal_offset = NaN(1,size(traces, 2));
    for trace_idx = 1:size(traces, 2)
        if valid(trace_idx)
            ideal_offset(trace_idx) = offset_range(find(best_offset(:,trace_idx) < 0,1,'first'));
            %figure(123);cla();plot(traces(:,trace_idx),'r');hold on;
            traces(:,trace_idx) =  traces(:,trace_idx) - prctile(traces(:,trace_idx), ideal_offset(trace_idx));
            %plot(traces(:,trace_idx),'k');pause(0.2)
        end
    end
    
    peak_ref = nanmedian(traces(peak_pos, :),2); %median_peak_amps
    ideal_gain = peak_ref \ traces(peak_pos, :);
    normalized              = traces ./ ideal_gain;
end