 function [suggested, n_clust] = test_epsilon(obj, Low_Dim_Data, current_signal, all_ROIs, valid_ROIs, n_clust)
    if nargin < 6 || isempty(n_clust)
        n_clust = [];
    end

    %% to test a range of epsilon value and identify best number
    method      = 'dbscan';
    rendering   = false;
    if strcmp(method, 'dbscan')
        test_range       = logspace(-1,2,150);
    else
        test_range       = 1:20;
    end

    n_gp        = [];
    n_noise_pt  = [];

    current_ep  = 1;
    MIN_CLUSTER_SIZE = obj.dbscan_min_gp_size %2*size(Y_PHATE_3D,2); % default (Ester et al., 1996), although we may want 2x NDim for High dimesnional data  (Sander et al., 1998)
    gp          = [];
    tolerance   = 3;
    verbose = false;
    while current_ep < numel(test_range)
        try
            test_value  = test_range(current_ep);

            %% Get cluster
            if strcmp(method, 'dbscan')
                if verbose
                    fprintf(['testing epsilon = ',num2str(test_value), '\n']);
                end
                cluster_idx = dbscan(Low_Dim_Data , test_value, MIN_CLUSTER_SIZE); % Minpts from (Sander et al., 1998)
            elseif strcmp(method, 'hdbscan')
                if verbose
                    fprintf(['testing Min cluster size = ',num2str(test_value), '\n'])
                end
                clusterer = HDBSCAN(Low_Dim_Data);
                %% clusterer.run_hdbscan(minpts, minclustsize, minClustNum, dEps, outlierThresh, plotResults)
                clusterer.run_hdbscan(2,1,1,test_value,0,false);%whitebg('w'); hold on;set(gcf,'Color','w')
                cluster_idx = clusterer.labels;
            else
                error('method not implemented')
            end

            %colors = 1:sum(cluster_idx > 0);
            prev_gp = gp;
            [gp ,~, indices] = unique(cluster_idx(cluster_idx > 0));
            if numel(gp) < numel(prev_gp) && all(gp == 1) % to avoid stpping at a local min
                tolerance = tolerance - 1;
            end
            if numel(gp) < numel(prev_gp) && all(gp == 1) && ~tolerance
                break
            end

            colors = viridis(numel(gp));
            colors = colors(randperm(size(colors,1), size(colors,1)),:);
            colors = colors(indices, :);

            n_gp(current_ep) = numel(gp);
            n_noise_pt(current_ep) = sum(cluster_idx <= 0);


            if rendering
                figure(1);clf();title(num2str(test_value));
                subplot(1,2,1);
                scatter3(Low_Dim_Data(cluster_idx <= 0,1), Low_Dim_Data(cluster_idx <= 0,2), Low_Dim_Data(cluster_idx <= 0,3), 30, 'MarkerFaceColor' , [0.8,0.8,0.8], 'MarkerEdgeColor' , 'none'); hold on;
                scatter3(Low_Dim_Data(cluster_idx > 0,1), Low_Dim_Data(cluster_idx > 0,2), Low_Dim_Data(cluster_idx > 0,3), 30, colors, 'filled'); hold on;
            end


            blank_v = NaN(size(all_ROIs));
            blank_v(valid_ROIs) = cluster_idx;
            blank_v(blank_v <= 0) = NaN;
            if obj.use_hd_data
                values = split_values_per_voxel(blank_v, obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1),find(~all(isnan(current_signal),2)));
            else
                values = blank_v;
            end
            updated_colors = repmat([0.8, 0.8, 0.8], numel(all_ROIs),1);
            updated_colors(~(isnan(blank_v)), :) = colors;
            if rendering
                s = subplot(1,2,2);
                obj.ref.plot_value_tree(values, '' ,'','','',s,'classic',updated_colors);
                %obj.ref.plot_value_tree(values, find(~all(isnan(current_signal),2)),'','','',s,'classic',colors);
                title(num2str(test_value));
                pause(2)
            end
        end
        current_ep = current_ep + 1;
    end

    if isempty(n_clust)
        max_loc = find((n_noise_pt/size(Low_Dim_Data, 1)) < 0.05, 1, 'first');
        n_clust = test_range(max_loc);
        %[n_clust, max_loc] = max(n_gp);
    else
        max_loc = find(n_gp == abs(n_clust), 1, 'last');
        if isempty(max_loc)
            [n_clust, max_loc] = max(n_gp);
        end
    end
%     suggested = knee_pt(test_range(max_loc:numel(n_noise_pt)), n_gp(max_loc:end));

    suggested = test_range(max_loc);
    if isempty(suggested)
        suggested = test_range(1);
    end
 end
