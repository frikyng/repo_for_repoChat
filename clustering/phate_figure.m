function cluster_idx = phate_figure(obj, Low_D_Data, epsilon, original_Data, Fig_number, ROIs, method)
    if nargin < 7 || isempty(method)
        method = 'dbscan';
    end

    MIN_CLUSTER_SIZE = 4 % default (Ester et al., 1996), although we may want 2x NDim for High dimesnional data  (Sander et al., 1998)

    %% Get cluster
    if strcmp(method, 'dbscan')
        cluster_idx     = dbscan(Low_D_Data , epsilon, MIN_CLUSTER_SIZE);  
    elseif strcmp(method, 'hdbscan')
        clusterer = HDBSCAN(Low_D_Data);
        %% clusterer.run_hdbscan(minpts, minclustsize, minClustNum, dEps, outlierThresh, plotResults)
        clusterer.run_hdbscan(MIN_CLUSTER_SIZE,MIN_CLUSTER_SIZE,1,[],epsilon,false);%whitebg('w'); hold on;set(gcf,'Color','w')
        cluster_idx = clusterer.labels;
    else
        error('method not implemented')
    end

    %% Get valid pts and cmap      
    outliers        = cluster_idx <= 0;
    valid_points    = cluster_idx > 0;
    current_cmap    = jet(numel(unique(cluster_idx(valid_points))));
    
    RANDOMIZE_COLORS = true;
    if RANDOMIZE_COLORS
        current_cmap = current_cmap(randperm(size(current_cmap, 1)),:);
    end
    %% Add gray
    current_cmap = [0.8,0.8,0.8;current_cmap];
    
    %% Prepare Figure
    figure(Fig_number);clf();set(gcf,'Units','normalized','Position',[0.05 0.05 0.9 0.9]);title(num2str(epsilon));
    
    %% Scatter subplot
    s1 = subplot(2,2,1);        
    scatter3(Low_D_Data(valid_points,1), Low_D_Data(valid_points,2), Low_D_Data(valid_points,3), 30, cluster_idx(cluster_idx > 0), 'filled'); hold on;
    colormap(current_cmap(2:end,:)); hold on;
    scatter3(Low_D_Data(outliers,1), Low_D_Data(outliers,2), Low_D_Data(outliers,3), 30, 'MarkerFaceColor' , [0.8,0.8,0.8], 'MarkerEdgeColor', 'none');
    
    %% Tree subplot
    % Tree can be low D or Hd depending on the type of analysis
    s2 = subplot(2,2,2);
    valid = ROIs(valid_points);
    if ~obj.use_hd_data
        %% For Low_D tree
        obj.ref.plot_value_tree(cluster_idx(valid_points), valid,'','cluster tree (one value per ROI)','',s2,'curved',current_cmap(2:end,:));
        %obj.ref.plot_value_tree(cluster_idx, find(~all(isnan(signal),2)),'','','',s2,'classic','jet');
    else
        %% FOR HD tree
        values = split_values_per_voxel(cluster_idx, obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), ROIs);
        obj.ref.plot_value_tree(values, '','','cluster tree (one value per voxel)','',s2,'curved',current_cmap);
    end
   
    %% Plot traces
    s3 = subplot(2,2,[3,4]);
    title('median traces per group')
    for gp = unique(cluster_idx(valid_points)')
       % if sum(cluster_idx == gp) > 50
            plot(nanmedian(original_Data(:,ROIs(cluster_idx == gp)),2));hold on;
        %end
    end
end