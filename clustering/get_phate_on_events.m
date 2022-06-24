function get_phate_on_events(Y_PHATE_3D, phates_idx, current_signal, n_pts)
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
    [sorted_L, loc]                = sort(L);
    
    data_subset = current_signal(loc(sorted_L > max(sorted_L)/2),:); % loc(end-100:end)
    T_PHATE_3D = phate(data_subset, 'ndim', 9, 't', []);
    cluster_idx     = dbscan(T_PHATE_3D , epsilon, 10); 
    
    
    colors = 1:size(current_signal, 2);
    colors(:) = 0;
    colors(7000:8000) = 1
    figure();%scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, C, 'filled'); hold on;
    plot(1:size(current_signal, 2), nanmean(data_subset),'k'); hold on;
    thr_color = cluster_idx;
    thr_color(thr_color > 1) = 2;
    scatter(1:size(current_signal, 2), nanmean(data_subset), 30, thr_color', 'filled'); 
    %scatter3(T_PHATE_3D(:,1), T_PHATE_3D(:,2), T_PHATE_3D(:,3), 30, cluster_idx, 'filled'); hold on;
end