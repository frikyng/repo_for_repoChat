R = {};
analyse_case = 'peaks'; % 'peaks' or 'rise'

for exp_idx = numel(meta.experiments):-1:1
    exp = meta.experiments(exp_idx);
    
    try
        
        
   % exp.detect_events
        
        
    peak_times = vertcat(exp.event.peak_time{:});
    peak_widths =  vertcat(exp.event.peak_width{:});
    [~, unique_peak] = unique(peak_times);
    tp = peak_times(unique_peak);
    peak_widths = peak_widths(unique_peak);
    
    %figure();hist(vertcat(exp.event.peak_width{:}),0:1:100)
    %win = floor(median(vertcat(exp.event.peak_width{:})));
    %valid = 
    
    
    rise = get_rise_time(exp);
    
    %% Somatic trace and events
    soma = exp.ref.indices.somatic_ROIs;
    soma_trace = nanmean(exp.rescaled_traces(:,soma),2)';
    if strcmp(analyse_case, 'peaks')
    	tp1 = find_local_max(soma_trace, tp, peak_widths);
    elseif strcmp(analyse_case, 'rise')
        tp1 = find(rise);
    end
    amps_soma = soma_trace(:,tp1);
    
%     rise = movmax(rise, 30);
%     starts = find(diff(rise) == 1);
%     stops = find(diff(rise) == -1);
%     amp = [];
%     duration = [];
%     for event = 1:numel(starts) 
%         amp(event) = nanmax(soma_trace(starts(event):stops(event)));
%         duration(event) = numel(starts(event):stops(event));
%     end
%    figure();scatter(duration ,amp, 'ko','filled')
    
    %% Distal median group
    temp_R = [];
    for bin = 1:numel(exp.binned_data.groups)
        distal_rois = exp.binned_data.groups{bin};
        distal_trace = exp.binned_data.median_traces(:,bin)';
        if strcmp(analyse_case, 'peaks')
        	tp2 = find_local_max(distal_trace, tp, peak_widths);
        elseif strcmp(analyse_case, 'rise')
            tp2 = tp1;
        end
        amps_test_gp = distal_trace(:,tp2);

        valid = ~isnan(amps_soma) & ~isnan(amps_test_gp);
        
        temp = corrcoef(amps_soma(valid), amps_test_gp(valid));
        temp_R(bin) = temp(2,1);
        
        name = exp.source_folder;
        name = strsplit(name, '/');
        name = name{end-1};
        name = strrep(name, '_', '\_');
        
        amp_difference = [amps_soma'-amps_test_gp'];
        if any(amp_difference)
        
            figure(123); clf();
            subplot(2,2,1); scatter(amps_soma,  amps_test_gp, [], 1:numel(amps_soma),'filled');hold on;colormap(viridis);xlabel('soma/proximal');ylabel('distal');title(num2str(temp_R(bin)))
            subplot(2,2,[2,4]);plot([soma_trace',distal_trace']); hold on; scatter([tp1,tp2],[amps_soma', amps_test_gp'],'filled');
            title([name, ' ; soma vs [',exp.binned_data.bin_legend{bin},' ]']); hold on;

            hold on;scatter(tp1(amp_difference < 0),amp_difference(amp_difference < 0)-max(amp_difference) ,'filled','vr');
            hold on;scatter(tp1(amp_difference > 0),amp_difference(amp_difference > 0)-max(amp_difference) ,'filled','^g');
            
            subplot(2,2,3);cla();hist(amp_difference, linspace(nanmin(amp_difference),nanmax(amp_difference), 10))
            pause(0.1);
        end
        
    end
    R{exp_idx} = temp_R;
    
    %P(exp_idx)

    end
end

% for expe_idx = 94:-1:1
% figure(123);cla();plot(R{expe_idx});
% pause(0.5);
% end

function local_max = find_local_max(current_signal, peak_locs, win)
    local_max = peak_locs;
    for pk = 1:numel(peak_locs)
        R = (peak_locs(pk)-floor(win(pk)/2)):(peak_locs(pk)+ceil(win(pk)/2));
        R = R(R > 0 & R < numel(current_signal));
        [~, local_max(pk)] = nanmax(current_signal(R));
        local_max(pk) = local_max(pk) + R(1) - 1;
    end
end