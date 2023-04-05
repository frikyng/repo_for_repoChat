classdef cluster_analysis < handle
    %% Subclass of arboreal_scan_experiment
    properties
        dbscan_min_gp_size = 4
    end

    methods
        function cluster_factors(obj, clust_meth, N_clust)
            if nargin < 2 || isempty(clust_meth)
                clust_meth                    = 'hierarchical';
            end
            if nargin < 3 || isempty(N_clust)
                N_clust                    = [];
            end
            
            %% Clear previous results
            obj.dimensionality.cluster_idx      = [];
            obj.dimensionality.clust_meth       = clust_meth;
            obj.dimensionality.N_clust        	= N_clust;
            obj.dimensionality.clust_groups     = {};
            obj.dimensionality.epsilon       	= [];
            obj.dimensionality.labels           = [];

            %% If N cluster is 0, skip clustering
            if obj.dimensionality.N_clust == 0
                obj.dimensionality.cluster_idx  = (1:size(obj.dimensionality.LoadingsPM, 1))';
                obj.dimensionality.sorted_idx   = (1:size(obj.dimensionality.LoadingsPM, 1))';
                obj.dimensionality.clust_groups = {1:size(obj.dimensionality.LoadingsPM, 1)};
                return
            end

            %% If N cluster is unknow, try to guess
            if isempty(obj.dimensionality.N_clust) || isnan(obj.dimensionality.N_clust)
                R = 1:20;
                if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'kmeans','silhouette','KList',R);
                    obj.dimensionality.N_clust = eva.OptimalK;
                    obj.disp_info(['Optimal number of clusters is ',num2str(obj.dimensionality.N_clust)],1)
                elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'linkage','silhouette','KList',R);
                    obj.dimensionality.N_clust = eva.OptimalK;
                    obj.disp_info(['Optimal number of clusters is ',num2str(obj.dimensionality.N_clust)],1)
                elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                    obj.dimensionality.epsilon = test_epsilon(obj, obj.dimensionality.LoadingsPM);
                    obj.disp_info(['Optimal espilon s ',num2str(obj.dimensionality.epsilon)])
                elseif strcmp(obj.dimensionality.clust_meth, 'strongest')
                    obj.dimensionality.N_clust = NaN;
                else
                    obj.disp_info('No auto auto-determination of the number of cluster for this method',4)
                end
            elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                if obj.dimensionality.N_clust > 0
                    obj.dimensionality.epsilon = obj.dimensionality.N_clust;
                    obj.dimensionality.N_clust = [];
                else
                    [obj.dimensionality.epsilon, obj.dimensionality.N_clust] = test_epsilon(obj, obj.dimensionality.LoadingsPM,[],[],[],obj.dimensionality.N_clust);
                end
            end
            
            if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                obj.disp_info('Clustering done using KMEANS',1)
                cluster_idx = kmeans(obj.dimensionality.LoadingsPM , obj.dimensionality.N_clust);
                figure(3);clf(); silhouette(obj.dimensionality.LoadingsPM,cluster_idx)
            elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical') % see https://fr.mathworks.com/help/stats/hierarchical-clustering.html
                obj.disp_info('Clustering done using hierarchical clustering',1)
                cluster_idx = clusterdata(obj.dimensionality.LoadingsPM,'Linkage', 'ward', 'MAXCLUST', obj.dimensionality.N_clust);%, 'Criterion','distance' 'MAXCLUST', 40)
                figure(3);clf(); silhouette(obj.dimensionality.LoadingsPM,cluster_idx)
                %Y = pdist(obj.dimensionality.LoadingsPM ,'euclidean');Z = linkage(Y,'ward');figure();dendrogram(Z);                                
            elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                obj.disp_info('Clustering done using dbscan',1)
                cluster_idx                 = dbscan(obj.dimensionality.LoadingsPM, obj.dimensionality.epsilon, obj.dbscan_min_gp_size, 'Distance', 'euclidean');
                obj.dimensionality.N_clust  = numel(unique(cluster_idx(cluster_idx > 0)));
                if any(cluster_idx <= 0)
                    col = UNASSIGNED_ROI_COLOR;
                else
                    col = [];
                end
                col                         = [col ; jet(obj.dimensionality.N_clust)];
                obj.dimensionality.labels   = [cluster_idx, cluster_idx, cluster_idx];
                count                       = 1;
                for v = unique(cluster_idx)'
                    subset                                  = obj.dimensionality.labels(:,1) == v;
                    obj.dimensionality.labels(subset, :)    = repmat(col(count,:),sum(subset),1);
                    count                                   = count + 1;
                end
                figure(1111);clf();
                for offset = 0:(min(size(obj.dimensionality.LoadingsPM, 2)-3, 2))
                    subplot(2,2,offset+1); hold on; grid; hold on
                    scatter3(obj.dimensionality.LoadingsPM (:,1+offset),obj.dimensionality.LoadingsPM (:,2+offset),obj.dimensionality.LoadingsPM (:,3+offset),20,obj.dimensionality.labels, 'filled');
                end               
            elseif strcmp(obj.dimensionality.clust_meth, 'strongest')
                obj.disp_info('No clustering, but assigning each ROI the value of its strongest component',1)
                %% To assign to strongest component
                for row = 1:size(LoadingsPM,1)
                    [~, maxloc]                     = max(LoadingsPM(row, :));
                    LoadingsPM(row, :)              = 0;
                    LoadingsPM(row, maxloc)         = 1;
                end
            else
                obj.dimensionality.cluster_idx      = [];
                obj.dimensionality.clust_groups     = {};
            end
            
            %% Sort clusters by number of elements
            if obj.dimensionality.N_clust > 0
                obj.disp_info('Cluster number is sorted by number of ROI',1)
                [~, gp] = sort(hist(cluster_idx,unique(cluster_idx)), 'descend');
            else
                gp = 1;
            end
            a = unique(cluster_idx)';
            gp = a(gp);
            idx_sorted = NaN(size(cluster_idx));
            count = 1;
            for gp_idx = gp(gp > 0)
                idx_sorted(cluster_idx == gp_idx) = count;
                count = count + 1;
            end
            idx_sorted(isnan(idx_sorted))       = 0;
            obj.dimensionality.cluster_idx      = idx_sorted;       % stored group ids (reordered by number of element)
            [~, obj.dimensionality.sorted_idx]  = sort(idx_sorted); % get ROI index to reorder the data by group

            %% Now, build groups
            obj.dimensionality.clust_groups     = {};
            for el = 1:obj.dimensionality.N_clust
                obj.dimensionality.clust_groups{el} = find(obj.dimensionality.cluster_idx == el)';
            end
            
            %% Plot clusters
            if obj.rendering                
                obj.plot_cluster_tree();
            end
            obj.disp_info({['Clustering of dimensionality-reduced data done using ',obj.dimensionality.clust_meth],...
                           [num2str(obj.dimensionality.N_clust), ' clusters have been identified']},1)
        end
    end
end
