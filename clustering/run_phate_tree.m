path            = '';
use_hd_data     = false;
time_filter     = 10;
type_of_trace   = 'subtracted_peaks'; %raw, rescaled, subtracted

%% Load object
if exist('obj', 'var')
    [obj, source_signal, signal_indices] = prepare_phate_analysis(obj, use_hd_data, time_filter, type_of_trace);
else
    [obj, source_signal, signal_indices] = prepare_phate_analysis(path, use_hd_data, time_filter, type_of_trace);
end

%% ##############
%% Define key Phate parameters
N_Dim = 9; % can set to # of ROIs in the object. Ideally, use a square number to have a nicer layout

%% Define if extraction done across ROIs, or across timepoints
analysis_mode = 'space';
if ~any(strcmp(analysis_mode, {'time','space'}))
    error('analysis_mode must be "time" or "space"')
end

%% list behaviours to test ('' for all timepoints)
% conditions = {'', 'active', 'quiet',  'encoder', '~encoder', 'BodyCam_L_whisker', '~BodyCam_L_whisker','BodyCam_L_whisker','~BodyCam_L_whisker','BodyCam_Trunk','~BodyCam_Trunk','EyeCam_L_forelimb','EyeCam_R_forelimb'}
% conditions = {'', 'encoder', '~encoder'};%, 'encoder_peaks', '~encoder_peaks', 'encoder_active', '~encoder_active'}
% conditions = {'encoder', '~encoder', 'quiet', 'active'};
%conditions = {'encoder_active','BodyCam_L_whisker_active'};
%conditions = {''};%{'BodyCam_L_whisker'}%{'encoder','~encoder'};

conditions = {'peaks'};%{'peaks_~trigger[5,5]','peaks_trigger[0,5]', 'peaks_trigger[5,0]'}

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


%% Now, for each condition, classify and cluster the data
Fig_count = 1000;
Y_PHATE_3D = {};
cluster_idx = {};
for el = 1:numel(conditions)
    %% Set/Restore whole signal
    current_signal = source_signal;
    
    %% Get timpoints for current behavioural condition
    if contains(type_of_trace, 'subtracted')&&  ~contains(conditions{el}, 'peaks')
        error_box('You should not use subtracted on full signals as kinetics may vary between ROIs. You should use peaks')
    end
    tp = obj.get_tp_for_condition(conditions{el});
    %     tp = diff(tp);
    %     tp(tp < 0) = 0;
    %     tp = [0, tp];

    %% Filter out unrequired timepoints
    current_signal(~tp, :) = NaN;

    %% Remove Nans
    all_ROIs        = 1:size(current_signal, 2);
    valid_ROIs      = ~all(isnan(current_signal),1);
    all_ROIs(all_ROIs > obj.n_ROIs) = [];
    all_ROIs(all_ROIs > obj.n_ROIs) = [];
    current_signal  = current_signal(:,valid_ROIs);
    current_signal  = current_signal(~all(isnan(current_signal),2),:);
    
    %% Define if we will extract info along space or time
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
    Y_PHATE_3D{el} = phate(current_signal, 'ndim', N_Dim, 't', []);
    close(gcf); 
    figure(Fig_count + 3000);clf(); title('Phate first 3 dimensions scatter plot')
    scatter3(Y_PHATE_3D{el}(:,1), Y_PHATE_3D{el}(:,2), Y_PHATE_3D{el}(:,3), 30); hold on;


    %% Display Phates on tree
    n_row = floor(sqrt(N_Dim));
    n_col = ceil(sqrt(N_Dim)) ;
    figure(Fig_count + 2000);clf()
    for dim = 1:N_Dim
        sub = subplot(n_row,n_col,dim);
        if obj.use_hd_data
            obj.ref.plot_value_tree(split_values_per_voxel(Y_PHATE_3D{el}(:,dim), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',sub,'curved','viridis');
        else
            obj.ref.plot_value_tree(Y_PHATE_3D{el}(:,dim), find(valid_ROIs),'',['phate #',num2str(dim),' Loadings (per ribbon)'],'',sub,'curved','viridis');
        end
    end
    
    figure();plot(Y_PHATE_3D{el})
    figure();imagesc(Y_PHATE_3D{el});colorbar;%caxis([-1 1]);
    
    % ESTIMATES for the maximum # of PHATE loadings to use
%     var_PHATE = var(Y_PHATE_3D{el});
%     figure();subplot(2,2,1); plot(var_PHATE); title('Var of each Loading')
%     
%     cumsum_PHATE = cumsum(mean(abs(Y_PHATE_3D{el}),1));
%     subplot(2,2,2); plot(cumsum_PHATE);title('Cum Sum of each mean(abs(Loading))')
%         
%     range_PHATE = range(Y_PHATE_3D{el});
%     subplot(2,2,3);plot(range_PHATE);title('Range for each Loading')
%     
%     cumsum_range_PHATE = cumsum(range(Y_PHATE_3D{el}));
%     subplot(2,2,4); plot(cumsum_range_PHATE); title('Cum Sum of each range(Loading)')
    
        %%testing if we can use specified PHATE dimensions only, and then do
        %%clustering to identify ROIs that are difficult to see between conditions
        % determine which PHATE dimensions 
        
        %New_PHATE_3D = Y_PHATE_3D{el}(:,[9,11,13,14,16,25]);

    %% estimate epsilon for subsequent (h)DBScan clustering      
    %original line
    kD = pdist2(Y_PHATE_3D{el},Y_PHATE_3D{el},'euc','Smallest',5);
    
    %for testing new PHATE (comment out if not using)
    %kD = pdist2(New_PHATE_3D,New_PHATE_3D,'euc','Smallest',5);
    
    kd_sorted = sort(kD(:))';
    slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
    [~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
    epsilon = kd_sorted(minloc);
    
    %% estimate epsilon from 'knee' ## Uncomment second line to see how clustering evolves with epsilon
    %epsilon = 1;
    %suggested = 1;
    
    suggested = test_epsilon(obj, Y_PHATE_3D{el}, current_signal, all_ROIs, valid_ROIs); %open test_epsilon and change rendering to true to see it
    %suggested = 3
    %for testing new PHATE (comment out if not using)
    %suggested = test_epsilon(obj, Y_PHATE_3D{el}, current_signal, all_ROIs, valid_ROIs);
    
    %% Show clusters, push on tree, show traces
    cluster_idx{el} = phate_figure(obj, Y_PHATE_3D{el}, mean([epsilon, suggested]), source_signal(tp,:), Fig_count, signal_indices);
    
    %for testing new PHATE (comment out if not using)
    %phate_figure(obj, New_PHATE_3D, mean([epsilon, suggested]), source_signal(tp,:), Fig_count, signal_indices);
    
    hold on;sgtitle(['Cluster for condition : ',strrep(conditions{el},'_','\_')])
    
    %% Rerun phate along time
%   which_phate = 3
%   loc = get_phate_on_events(Y_PHATE_3D{el}, which_phate, current_signal, 0.8, valid_ROIs, obj, tp) % half max phate
%     for which_phate = 1:9
%         loc = get_phate_on_events(Y_PHATE_3D{el}, which_phate, current_signal, 0.8, valid_ROIs, obj, tp) % half max phate
%     end

    %% Hierarchical clustering for reference    
%         figure(Fig_count + 1000);clf();
%         s1 = subplot(3,2,1); 
%         eva = evalclusters(Y_PHATE_3D{el},'linkage','silhouette','KList',1:100);
%         cluster_idx = clusterdata(Y_PHATE_3D{el},'Linkage', 'ward', 'MAXCLUST', eva.OptimalK);%, 'Criterion','distance' 'MAXCLUST', 40)
%         cluster_idx = clusterdata(Y_PHATE_3D{el},'Linkage', 'ward', 'MAXCLUST', 40);%, 'Criterion','distance' 'MAXCLUST', 40)
%         scatter3(Y_PHATE_3D{el}(:,1), Y_PHATE_3D{el}(:,2), Y_PHATE_3D{el}(:,3), 30, cluster_idx, 'filled'); hold on;
%         Z = linkage(Y_PHATE_3D{el},'ward');
%         s2 = subplot(3,2,2); dendrogram(Z,0,'ColorThreshold','default','Orientation','left');
%         s3 = subplot(3,2,3); 
%         values = split_values_per_voxel(cluster_idx, obj.ref.header.res_list(:,1), signal_indices);
%         obj.ref.plot_value_tree(values, '','','cluster tree (one value per voxel)','',s3,'curved',current_cmap);
%         s4 = subplot(3,2,4);cla();silhouette(Y_PHATE_3D{el},cluster_idx);hold on;axis fill
%         s5 = subplot(3,2,[5,6]);
%         for gp = unique(cluster_idx')
%             plot(nanmedian(source_signal(:,signal_indices(cluster_idx == gp)),2));hold on;
%         end

    Fig_count = Fig_count + 1    
end


%%
test_beh = 2;
is_matching     = false(1, N_Dim);
matching_thr    = 0.3;
for phate_dim = 1:N_Dim
    cc = corrcoef([Y_PHATE_3D{test_beh}(:,phate_dim),Y_PHATE_3D{1}]);cc = cc(1,:);cc(cc == 1) = NaN;
    figure();
    sub = subplot(2,2,1);
    obj.ref.plot_value_tree(Y_PHATE_3D{test_beh}(:,phate_dim), find(valid_ROIs),'',['phate #',num2str(phate_dim),' in current behaviour'],'',sub,'curved','viridis');
    sub = subplot(2,2,2);
    [m, best_match] = max(abs(cc));best_match = best_match - 1; % -1 because we added a ref in front of the cc matrix
    is_matching(phate_dim) = m > matching_thr;
    obj.ref.plot_value_tree(Y_PHATE_3D{1}(:,best_match), find(valid_ROIs),'',['phate #',num2str(phate_dim),' best match in ref is phate #',num2str(best_match),' with corr coeff == ',num2str(round(m, 2))],'',sub,'curved','viridis');
    subplot(2,2,[3,4]);
    %plot(cc)
    plot(Y_PHATE_3D{test_beh}(:,phate_dim));hold on; plot(Y_PHATE_3D{1}(:,best_match))
end

specific = find(~is_matching);
figure();
for el = 1:numel(specific)    
    sub = subplot(1,numel(specific),el);
    obj.ref.plot_value_tree(Y_PHATE_3D{test_beh}(:,specific(el)), find(valid_ROIs),'',['phate #',num2str(specific(el)),' is specific to this behaviour'],'',sub,'curved','viridis');
end


    
%     cluster_idx = dbscan(Y_PHATE_3D{el} , nanmedian(epsilon), 5);
%     f = figure();
%     subplot(1,2,1);
%     cluster_idx = dbscan(Y_PHATE_3D{el} , nanmedian(epsilon), 3);
%     scatter3(Y_PHATE_3D{el}(:,1), Y_PHATE_3D{el}(:,2), Y_PHATE_3D{el}(:,3), 30, cluster_idx, 'filled');
%     colormap(jet);
%     s = subplot(1,2,2);
%     obj.ref.plot_value_tree(cluster_idx, find(~all(isnan(signal),2)),'','','',s,'classic','jet');
%     title(num2str(e));
    

%     eva = evalclusters(Y_PHATE_3D{el},'kmeans','silhouette','KList',1:20);
% 	R = 0.01:0.01:5
%     ok = []
%     for v = 1:numel(R)
%         cluster_idx = dbscan(Y_PHATE_3D{el} , R(v), 5);
%         ok(v) = sum(cluster_idx == -1);
%     end

%     close all



