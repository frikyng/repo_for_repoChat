
% Load OBJECT when already correct
%% load('C:\Users\vanto\Documents\MATLAB\extracted_arboreal_scans 2\arboreal_scans_thin_mask.mat')

%% Rebuild object with HD data in it 
obj = arboreal_scan_experiment('C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans_2\extracted_arboreal_scans\2019-11-05_exp_3',true)

%% ####### run once #######
%% run once to get event time from low D
%% Define wether to use HD data or LD data
obj.use_hd_data = false;

%% Quick processing to have proper event detection. This is required for filtering based on activity
obj.prepare_binning({'depth',50});
obj.find_events();
obj.rescale_traces();
obj.set_median_traces();
obj.compute_similarity();

%% Define wether to use HD data or LD data
obj.use_hd_data = true;
try
    obj.rescale_traces(); %% run twice .. to fix
end
obj.rescale_traces(); %% run once
%% ##############

%% Define key Phate parameters
N_Dim = 9

%% Define if extraction done across ROIs, or across timepoints
analysis_mode = 'space'
if ~any(strcmp(analysis_mode, {'time','space'}))
    error('analysis_mode must be "time" or "space"')
end

%% list behaviours to test ('' for all timepoints)
%conditions = {'', 'active', 'quiet',  'encoder', '~encoder', 'BodyCam_L_whisker', '~BodyCam_L_whisker','BodyCam_L_whisker','~BodyCam_L_whisker','BodyCam_Trunk','~BodyCam_Trunk','EyeCam_L_forelimb','EyeCam_R_forelimb'}
% conditions = {'', 'encoder', '~encoder'};%, 'encoder_peaks', '~encoder_peaks', 'encoder_active', '~encoder_active'}
% conditions = {'encoder', '~encoder', 'quiet', 'active'};
conditions = {''}

%% List of typical conditions.
%% type obj.behaviours.types' to get valid entries
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


%% Time filtering of traces if required. 
%% ## ! ## This introduces temporal correlation
obj.filter_win = [10, 0];

%% If using all voxels, remove bad_ROIs_list field because it is designed for full segments
if obj.use_hd_data    
    obj.bad_ROI_list = [];
end

%% select signal source (RAW or rescaled traces)
source_signal = obj.rescaled_traces;
%source_signal = obj.extracted_traces_conc;

%% Don't keep any Inf values, if any (only happens in HD case, occasionally)
source_signal(isinf(source_signal)) = NaN;

Fig_count = 1000;

%% Flag ROIs that have NaN vaues at one point as they may mess up later computations
bad_ROI_list                    = find(any(isnan(source_signal),1));
bad_ROI_list(bad_ROI_list > obj.n_ROIs) = [];
signal_indices                  = true(1, obj.n_ROIs); %% ROIs or voxels, depending on the data source
signal_indices(bad_ROI_list)    = false;
signal_indices                  = find(signal_indices);

    
%% Now, for each condition, classify and cluster the data
for el = 1:numel(conditions)
    %% Get timpoints for current behavioural condition
    tp = obj.get_tp_for_condition(conditions{el});
    
    %% Get signal
    current_signal = double(source_signal(:, 1:obj.n_ROIs));
   
    %% Filter out unrequired timepoints
    current_signal(~tp, :) = NaN;
    
    %% Filter out bad ROIs/voxels
    current_signal(:, bad_ROI_list) = NaN; 

    %% Remove Nans
    all_ROIs        = 1:size(current_signal, 2);
    valid_ROIs      = ~any(isnan(current_signal),1);
    all_ROIs(all_ROIs > obj.n_ROIs) = [];
    all_ROIs(all_ROIs > obj.n_ROIs) = [];
    current_signal  = current_signal(:,valid_ROIs);
    current_signal  = current_signal(~any(isnan(current_signal),2),:);
    
    %% Define if we will extract infor along space or time
    if strcmp(analysis_mode, 'space')
        current_signal = current_signal';
    end

    %% Assign timepoint color code
    colors = 1:size(current_signal, 1);

    %% PCA
    %     Y_PCA = svdpca(signal, 2, 'random');
    %     figure()
    %     scatter(Y_PCA(:,1), Y_PCA(:,2), 10, colors, 'filled');
    %     colormap(viridis); set(gca,'xticklabel',[]); set(gca,'yticklabel',[]);
    %     axis tight; xlabel 'PCA1'; ylabel 'PCA2'
    %     drawnow


    %% tSNE
    %     Y_tSNE = tsne(signal,'Theta',0.5,'Verbose',2, 'perplexity', 20);

    %% PHATE 3D
    Y_PHATE_3D = phate(current_signal, 'ndim', N_Dim, 't', []);
    close(gcf); figure(Fig_count + 3000);clf(); title('Phate first 3 dimensions scatter plot')
    scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30); hold on;


    %% Display Phates on tree
    n_row = floor(sqrt(size(Y_PHATE_3D, 2)));
    n_col = ceil(sqrt(size(Y_PHATE_3D, 2)));
    figure(Fig_count + 2000);clf()
    for dim = 1:size(Y_PHATE_3D, 2)
        sub = subplot(n_row,n_col,dim);
        obj.ref.plot_value_tree(split_values_per_voxel(Y_PHATE_3D(:,dim), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',sub,'curved','viridis');
    end

    %% (h)DBScan  clustering      
    kD = pdist2(Y_PHATE_3D,Y_PHATE_3D,'euc','Smallest',5);
    kd_sorted = sort(kD(:))';
    slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
    [~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
    epsilon = kd_sorted(minloc)*2;
    
    %% Uncomment to see how clustering evolve with epsilon
    %     suggested = test_epsilon(obj, Y_PHATE_3D, current_signal, all_ROIs, valid_ROIs); 
    %     epsilon = 6

    %% Show clusters, push on tree, show traces
    phate_figure(obj, Y_PHATE_3D, epsilon, source_signal(tp,:), Fig_count, signal_indices);
    hold on;sgtitle(['Cluster for condition : ',strrep(conditions{el},'_','\_')])
    
    %% Rerun phate along time
    which_phate = 3
    loc = get_phate_on_events(Y_PHATE_3D, which_phate, current_signal, 0.5) % half max phate

    %% Hierarchical clustering for reference    
    %     figure(Fig_count + 1000);clf();
    %     s1 = subplot(3,2,1); 
    %     eva = evalclusters(Y_PHATE_3D,'linkage','silhouette','KList',1:100);
    %     cluster_idx = clusterdata(Y_PHATE_3D,'Linkage', 'ward', 'MAXCLUST', eva.OptimalK);%, 'Criterion','distance' 'MAXCLUST', 40)
    %     cluster_idx = clusterdata(Y_PHATE_3D,'Linkage', 'ward', 'MAXCLUST', 40);%, 'Criterion','distance' 'MAXCLUST', 40)
    %     scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30, cluster_idx, 'filled'); hold on;
    %     Z = linkage(Y_PHATE_3D,'ward');
    %     s2 = subplot(3,2,2); dendrogram(Z,0,'ColorThreshold','default','Orientation','left');
    %     s3 = subplot(3,2,3); 
    %     values = split_values_per_voxel(cluster_idx, obj.ref.header.res_list(:,1), signal_indices);
    %     obj.ref.plot_value_tree(values, '','','cluster tree (one value per voxel)','',s3,'curved',current_cmap);
    %     s4 = subplot(3,2,4);cla();silhouette(Y_PHATE_3D,cluster_idx);hold on;axis fill
    %     s5 = subplot(3,2,[5,6]);
    %     for gp = unique(cluster_idx')
    %         plot(nanmedian(source_signal(:,signal_indices(cluster_idx == gp)),2));hold on;
    %     end

    Fig_count = Fig_count + 1    
end


    

    
    
    
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



