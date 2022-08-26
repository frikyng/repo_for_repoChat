R = {};

for exp_idx = numel(meta.experiments):-1:1
    exp = meta.experiments(exp_idx);
    
    try
    tp = unique(vertcat(exp.event.peak_time{:}));
    
 
    
    %% Somatic trace and events
    soma = exp.ref.indices.somatic_ROIs;
    soma_trace = nanmean(exp.rescaled_traces(:,soma),2)';
    %tp1 = find_local_max(soma_trace, tp, win);
    tp1 = find(rise);
    amps_soma = soma_trace(:,tp1);
    figure(123);cla();hist(exp.event.width, 0:5:200)
    
    
    %% Distal median group
    temp_R = [];
    for bin = 1:numel(exp.binned_data.groups)
        distal_rois = exp.binned_data.groups{bin};
        distal_trace = exp.binned_data.median_traces(:,bin)';
        %tp2 = find_local_max(distal_trace, tp, win);
        tp2 = tp1;
        amps_test_gp = distal_trace(:,tp2);

        valid = ~isnan(amps_soma) & ~isnan(amps_test_gp);
        
        temp = corrcoef(amps_soma(valid), amps_test_gp(valid));
        temp_R(bin) = temp(2,1);
    end
    R{exp_idx} = temp_R;
    
    %P(exp_idx)
%     figure(123); clf(); subplot(1,2,1); scatter(amps_soma,  amps_distal, [], 1:numel(amps_soma),'filled');hold on;colormap(viridis);xlabel('soma/proximal');ylabel('distal')
%     subplot(1,2,2);plot([soma_trace',distal_trace']); hold on; scatter([tp1,tp2],[amps_soma', amps_distal'],'filled');
%     hold on;scatter(tp1,[amps_soma'-amps_distal']-50,'filled','k')
%     pause(0.5);
    end
end