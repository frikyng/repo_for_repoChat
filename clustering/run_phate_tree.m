%% Running
%     {'encoder'          }
%     {'RT3D_MC'          }
%     {'BodyCam_Laser'    }
%     {'BodyCam_Eye'      }
%     {'BodyCam_L_whisker'}
%     {'BodyCam_R_whisker'}
%     {'BodyCam_Trunk'    }
%     {'EyeCam_Wheel'     }
%     {'EyeCam_Perioral'  }
%     {'EyeCam_Whiskerpad'}
%     {'EyeCam_L_forelimb'}
%     {'EyeCam_R_forelimb'}

Fig_count = 1000;

%% list behaviours to test ('' for all timepoints)
%conditions = {'', 'active', 'quiet',  'encoder', '~encoder', 'BodyCam_L_whisker', '~BodyCam_L_whisker','BodyCam_L_whisker','~BodyCam_L_whisker','BodyCam_Trunk','~BodyCam_Trunk','EyeCam_L_forelimb','EyeCam_R_forelimb'}
conditions = {'', 'encoder', '~encoder'}%, 'encoder_peaks', '~encoder_peaks', 'encoder_active', '~encoder_active'}
conditions = {'encoder', '~encoder', 'quiet', 'active'}

%% Timefiltering if required
obj.filter_win = [10, 0];

%% sect source (RAW or rescaled traces)
source = obj.rescaled_traces;;%obj.extracted_traces_conc;%

%% Don't keep any Inf
source(isinf(source)) = NaN;

%% If using all voxels, remove bad ROIs list because it is designed for full segments
obj.bad_ROI_list = [];

%% Flag ROIs that have NaN vaues at one point as they may mess up later computations
bad_ROI_list                    = find(any(isnan(source),1));
signal_indices                  = true(1, obj.n_ROIs); %% ROIs or voxels, depending on the data source
signal_indices(bad_ROI_list)    = false;
signal_indices                  = find(signal_indices);

%% Now, for each condition, classify and cluster the data
for el = 1:numel(conditions)
    %% Get timpoints for condition
    tp = obj.get_tp_for_condition(conditions{el});
    
    %% Get signal
    signal = double(source)';
   
    %% Filter out unrequired timepoints
    signal(:, ~tp) = NaN;
    
    %% Filter out bad ROIs/voxels
    signal(bad_ROI_list, :) = NaN; 
    
    %% Remove Nans
    signal = signal(~all(isnan(signal),2),:);
    signal = signal(:,~all(isnan(signal),1));

    %% Assign color code
    C = 1:size(signal, 1);

    %% PCA
    %     Y_PCA = svdpca(signal, 2, 'random');
    %     figure()
    %     scatter(Y_PCA(:,1), Y_PCA(:,2), 10, C, 'filled');
    %     colormap(viridis); set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
    %     axis tight; xlabel 'PCA1'; ylabel 'PCA2'
    %     drawnow


    %% tSNE
    %     Y_tSNE = tsne(signal,'Theta',0.5,'Verbose',2, 'perplexity', 20);

    %% PHATE 3D
    Y_PHATE_3D = phate(signal, 'ndim', 20, 't', [], 'pot_eps', 1e-5);
    close(gcf)
    scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30); hold on;

    %     for e = logspace(0.5,1.5,20)
    %         figure(1);clf();title(num2str(e));
    %         subplot(1,2,1);
    %         cluster_idx = dbscan(Y_PHATE_3D , e, 3);
    %         temp_scatter = Y_PHATE_3D(cluster_idx ~= -1,:);        
    %         scatter3(temp_scatter(:,1), temp_scatter(:,2), temp_scatter(:,3), 30, cluster_idx(cluster_idx ~= -1), 'filled'); hold on;
    %         colormap(jet);hold on;
    %         scatter3(Y_PHATE_3D(cluster_idx == -1,1), Y_PHATE_3D(cluster_idx == -1,2), Y_PHATE_3D(cluster_idx == -1,3), 30, 'MarkerFaceColor' , [0.8,0.8,0.8]); hold on;
    %         s = subplot(1,2,2);
    %         obj.ref.plot_value_tree(cluster_idx, find(~all(isnan(signal),2)),'','','',s,'classic','jet');
    %         title(num2str(e));
    %         pause(2)
    %     end

    
    
    %% (h)DBScan  clustering   
    
    phate_figure(obj, Y_PHATE_3D, epsilon, source(tp,:), Fig_count, signal_indices);
    hold on;sgtitle(['Cluster for condition : ',strrep(conditions{el},'_','\_')])

    %% Hierarchical clustering     
    figure(Fig_count + 1000);clf();
    s1 = subplot(3,2,1); 
    eva = evalclusters(Y_PHATE_3D,'linkage','silhouette','KList',1:100);
    cluster_idx = clusterdata(Y_PHATE_3D,'Linkage', 'ward', 'MAXCLUST', eva.OptimalK);%, 'Criterion','distance' 'MAXCLUST', 40)
    cluster_idx = clusterdata(Y_PHATE_3D,'Linkage', 'ward', 'MAXCLUST', 40);%, 'Criterion','distance' 'MAXCLUST', 40)
    scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30, cluster_idx, 'filled'); hold on;
    Z = linkage(Y_PHATE_3D,'ward');
    s2 = subplot(3,2,2); dendrogram(Z,0,'ColorThreshold','default','Orientation','left');
    s3 = subplot(3,2,3); 
    values = split_values_per_voxel(cluster_idx, obj.ref.header.res_list(:,1), signal_indices);
    obj.ref.plot_value_tree(values, '','','cluster tree (one value per voxel)','',s3,'curved',current_cmap);
    s4 = subplot(3,2,4);cla();silhouette(Y_PHATE_3D,cluster_idx);hold on;axis fill
    s5 = subplot(3,2,[5,6]);
    for gp = unique(cluster_idx')
        plot(nanmedian(source(:,signal_indices(cluster_idx == gp)),2));hold on;
    end

    Fig_count = Fig_count + 1    
end


    

%     figure(); hold on
%     epsilon = [];
%     for m = 1:50
%         
%         out = squareform(pdist(Y_PHATE_3D));
%         kD = pdist2(out(:,1),out(:,2),'euc','Smallest',m);
%         %kD = pdist2(Y_PHATE_2D(:,1),Y_PHATE_2D(:,2),'euc','Smallest',m);
%         %         plot(sort(kD(end,:)));
%         %         title('k-distance graph')
%         %         xlabel('Points sorted with 50th nearest distances')
%         %         ylabel('50th nearest distances')
%         %         grid 
%         
%         kd_sorted = sort(kD(end,:));
%         slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
%         [~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
%         epsilon(m) = kd_sorted(minloc)    
%     end
%     plot(epsilon)
%     P = polyfit(epsilon,1:numel(epsilon),1);
%     P = P(1)
%     
%     

%     plot(sort(kD(end,:)));
%     title('k-distance graph')
%     xlabel('Points sorted with 50th nearest distances')
%     ylabel('50th nearest distances')
%     grid 
    
    
    
    
%     cluster_idx = dbscan(Y_PHATE_3D , nanmedian(epsilon), 5);
%     f = figure();
%     subplot(1,2,1);
%     cluster_idx = dbscan(Y_PHATE_3D , nanmedian(epsilon), 3);
%     scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30, cluster_idx, 'filled');
%     colormap(jet);
%     s = subplot(1,2,2);
%     obj.ref.plot_value_tree(cluster_idx, find(~all(isnan(signal),2)),'','','',s,'classic','jet');
%     title(num2str(e));
    

%     eva = evalclusters(Y_PHATE_3D,'kmeans','silhouette','KList',1:20);
% 	R = 0.01:0.01:5
%     ok = []
%     for v = 1:numel(R)
%         cluster_idx = dbscan(Y_PHATE_3D , R(v), 5);
%         ok(v) = sum(cluster_idx == -1);
%     end

%     close all



