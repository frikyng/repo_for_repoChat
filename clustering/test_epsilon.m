 function suggested = test_epsilon(obj, Y_PHATE_3D, current_signal, all_ROIs, valid_ROIs)
    %% to test a range of epsilon value and identify best number    
    
    rendering   = false
    range       = logspace(-1,2,30);
    
    v1 = [];
    v2 = [];
    
    current_ep = 1
    gp = [];
    while (numel(gp) > 1 || current_ep <=1) && current_ep < 30        
        ep = range(current_ep);
        fprintf(['testing epsilon = ',num2str(ep), '\n'])
        cluster_idx = dbscan(Y_PHATE_3D , ep, 3);

        %colors = 1:sum(cluster_idx > 0);
        [gp ,~, indices] = unique(cluster_idx(cluster_idx > 0));
        colors = viridis(numel(gp));
        colors = colors(randperm(size(colors,1), size(colors,1)),:);
        colors = colors(indices, :);
        
        v1(current_ep) = numel(gp);
        v2(current_ep) = sum(cluster_idx <= 0);      
        

        if rendering
            figure(1);clf();title(num2str(ep));
            subplot(1,2,1);
            scatter3(Y_PHATE_3D(cluster_idx <= 0,1), Y_PHATE_3D(cluster_idx <= 0,2), Y_PHATE_3D(cluster_idx <= 0,3), 30, 'MarkerFaceColor' , [0.8,0.8,0.8], 'MarkerEdgeColor' , 'none'); hold on;
            scatter3(Y_PHATE_3D(cluster_idx > 0,1), Y_PHATE_3D(cluster_idx > 0,2), Y_PHATE_3D(cluster_idx > 0,3), 30, colors, 'filled'); hold on;
        end
        

        blank_v = NaN(size(all_ROIs));
        blank_v(valid_ROIs) = cluster_idx;
        blank_v(blank_v <= 0) = NaN;
        values = split_values_per_voxel(blank_v, obj.ref.header.res_list(:,1),find(~all(isnan(current_signal),2)));
        updated_colors = repmat([0.8, 0.8, 0.8], numel(all_ROIs),1);
        updated_colors(~(isnan(blank_v)), :) = colors;
        if rendering
            s = subplot(1,2,2);
            obj.ref.plot_value_tree(values, '' ,'','','',s,'classic',updated_colors);
            %obj.ref.plot_value_tree(values, find(~all(isnan(current_signal),2)),'','','',s,'classic',colors);
            title(num2str(ep));
            pause(2)
        end
        
        current_ep = current_ep + 1;
    end
    
%     figure();plot(v1); hold on;
%     plot(v2)
    suggested = knee_pt(range(1:numel(v2)), v1);
 end

 
 
%  
%     figure(); hold on
%     epsilon = [];
%     for m = 1:50
%         
%         % out = squareform(pdist(Y_PHATE_3D));
%         kD = pdist2(Y_PHATE_3D,Y_PHATE_3D,'euc','Smallest',m);
%         %kD = pdist2(Y_PHATE_2D(:,1),Y_PHATE_2D(:,2),'euc','Smallest',m);
%         %         plot(sort(kD(end,:)));
%         %         title('k-distance graph')
%         %         xlabel('Points sorted with 50th nearest distances')
%         %         ylabel('50th nearest distances')
%         %         grid 
%         
%         kd_sorted = sort(kD(:))';
%         slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
%         [~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
%         epsilon(m) = kd_sorted(minloc)    
%     end
%     plot(epsilon)
%     P = polyfit(epsilon,1:numel(epsilon),1);
%     P = P(1)
%     
%     
%     figure();
%     plot(sort(kD(:)));
%     title('k-distance graph')
%     xlabel('Points sorted with 50th nearest distances')
%     ylabel('50th nearest distances')
%     grid 
%     