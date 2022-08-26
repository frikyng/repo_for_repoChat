function cluster_idx = get_phate_on_events(Y_PHATE_3D, phates_idx, current_signal, frac, valid_ROIs, obj, tp)
    figure();
    for w = 1:size(Y_PHATE_3D,2)
        L                       = Y_PHATE_3D(:,w);
        [~, loc]                = sort(L);
        %plot(nanmean(signal(loc(end-100:end),:)) - nanmedian(signal)); hold on
        plot(nanmean(current_signal(loc(end-100:end),:))); hold on
%         %L(L<0.2) = 0;
%         all_weights{w}          = L/sum(Y_PHATE_3D(:,w));
%         weighted_averages(w, :) = nanmean(signal.* all_weights{w}, 1);
%         plot(weighted_averages(w, :)); hold on
    end
    
    
    L                       = Y_PHATE_3D(:,phates_idx);
    [sorted_L, loc]         = sort(L);
    
    data_subset = current_signal(loc(sorted_L > max(sorted_L) * frac),:); % loc(end-100:end)

    T_PHATE_3D = phate(data_subset', 'ndim', 9, 't', []);
    
    
    
    
    
%     sub_PHATE_3D = phate(data_subset, 'ndim', 3, 't', []);
%     figure();scatter3(sub_PHATE_3D(:,1), sub_PHATE_3D(:,2), sub_PHATE_3D(:,3), 30, 'filled'); hold on;
    
    
%     kD = pdist2(T_PHATE_3D,T_PHATE_3D,'euc','Smallest',5);
%     kd_sorted = sort(kD(:))';
%     slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
%     [~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
%     epsilon = kd_sorted(minloc)  ;
    
    
    for ep = logspace(-1,2,30)
        cluster_idx     = dbscan(T_PHATE_3D , ep, 10); 
        figure(19);clf();
        scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, cluster_idx, 'filled'); hold on;
        title(num2str(ep))
    end
    
    cluster_idx     = dbscan(T_PHATE_3D , 2.4, 10); 
    [gp, n_in_gp, ic ] = unique(cluster_idx);
    thr_color = jet(numel(gp));
    thr_color = thr_color(randperm(numel(gp), numel(gp)), :);
    thr_color = thr_color(ic, :);
    thr_color(ic == 1, :) = repmat([0.8,0.8,0.8],sum(ic == 1),1);
    
    [GC,GR] = groupcounts(cluster_idx);
    [GC, idx] = sort(GC);
    GR = GR(idx);    
    for el = GR(GC > 1000 & GR > 0)'
        thr_color(ic == el+1, :) = repmat([0.8,0.8,0.8],sum(ic == el+1),1);
    end

    
    %% if this crashes, comment out this block ##########
    
    f = figure();%scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, C, 'filled'); hold on;
    ax1 = subplot(3,1,1);
    plot(1:size(current_signal, 2), nanmean(data_subset),'k'); hold on;
    %thr_color(thr_color > 1) = 2;
    %thr_color = 1:size(current_signal, 2)
    sc = scatter(1:size(current_signal, 2), nanmean(data_subset), 30, thr_color, 'filled'); hold on
    ax2 = subplot(3,1,2);
    plot(data_subset');
    ax3 = subplot(3,1,3);
    [loco, loco2, loco3] = obj.get_behaviours('encoder');
    plot(loco3.value);
    linkaxes([ax1, ax2, ax3],'xy')
    
%     f = figure();%scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, C, 'filled'); hold on;
%     ax1 = subplot(3,1,1);
%     plot(1:size(current_signal, 2), nanmean(data_subset),'k'); hold on;
%     %thr_color(thr_color > 1) = 2;
%     %thr_color = 1:size(current_signal, 2)
%     sc = scatter(1:size(current_signal, 2), nanmean(data_subset), 30, thr_color, 'filled'); hold on
%     ax2 = subplot(3,1,2);
%     plot(data_subset');
%     ax3 = subplot(3,1,3);
%     plot(loco3.value);
%     linkaxes([ax1, ax2, ax3],'xy')
    %%  ##########

    %% to NOT LOAD DATA, comment out this block ##########
%     f = figure();%scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, C, 'filled'); hold on;
%     plot(1:size(current_signal, 2), nanmean(data_subset),'k'); hold on;
%     ROI_size = unique(obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1));
%     if numel(ROI_size) > 1
%         error('variable length not implemented yet')
%     end
%     actual_original_ROIs = find(valid_ROIs);    
%     strong_original_ROIs = actual_original_ROIs(loc(sorted_L > max(sorted_L) * frac));
%     weak_original_ROIs = actual_original_ROIs(loc(sorted_L < min(sorted_L) * frac));
%     strong_ROIs_to_load = {unique(ceil(strong_original_ROIs / ROI_size))};
%     weak_ROIs_to_load = {unique(ceil(weak_original_ROIs / ROI_size))};
%     load_specific_trial(obj, tp, [strong_ROIs_to_load, weak_ROIs_to_load])  
    %%  ##########
end