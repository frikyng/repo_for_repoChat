%% Build group by group Correlations
obj.prepare_binning({'depth',50});
obj.rescale_traces();
close all
for beh = {'peaks_groups','peaks_groups~BodyCam_L_whisker','groups_peaks_~encoder','groups_peaks_encoder'}
    beh = beh{1}
    t = obj.get_correlations(beh);
    saveas(figure(1018), ['tree_by_group_',beh,'.png'])
    saveas(figure(1008), ['crosscorr_by_group_',beh,'.png'])
end


%% Build ROI by ROI orrelations
obj.use_hd_data = false;
obj.prepare_binning({'depth',50});
obj.find_events();
obj.rescale_traces();
obj.set_median_traces();
obj.compute_similarity();
obj.fit_events();
close all
data2 = obj.rescaled_traces;   
med2 = nanmedian(data2,2);
for beh = {'peaks','peaks_~BodyCam_L_whisker','groups_peaks_~encoder','groups_peaks_encoder'}
    beh = beh{1}
    t = obj.get_tp_for_condition(beh);
    cross_corr  = corrcoef(data2(t,:) - med2(t));      
    to_keep = find(~all(isnan(cross_corr),1));
    cross_corr = cross_corr(to_keep, to_keep);
    eva = evalclusters(cross_corr,'linkage','silhouette','KList',1:30);
    cluster_idx = clusterdata(cross_corr,'Linkage', 'ward', 'MAXCLUST', eva.OptimalK);unique(cluster_idx)
    [a,b] = sort(cluster_idx);
    figure(1008);clf();imagesc(cross_corr(b,b))
    obj.ref.plot_value_tree(cluster_idx,to_keep,'','','',124,'','jet');
    saveas(figure(124), ['tree_LD_',beh,'_Nclust=',num2str(eva.OptimalK),'.png'])
    saveas(figure(1008), ['crosscorr_LD_',beh,'_Nclust=',num2str(eva.OptimalK),'.png'])
end

%% Build HD trees (you must detected event before hand from low D)
obj.use_hd_data = true;

% update excluded ROIs         
bad_ROI_list                    = find(any(isnan(obj.extracted_traces_conc),1));
bad_ROI_list(bad_ROI_list > obj.n_ROIs) = [];
signal_indices                  = true(1, obj.n_ROIs); %% ROIs or voxels, depending on the data source
signal_indices(bad_ROI_list)    = false;
signal_indices                  = find(signal_indices);

obj.rescale_traces();
data = obj.rescaled_traces;   
med = nanmedian(data,2);



% %% If you have a correlation matrix "cross_corr"
% [S,Q] = genlouvain(double(cross_corr),[],[],1); % this do the community detection
% [~,b] = sort(S);
% figure(1008);clf();imagesc(cross_corr(b,b))
% 
% %% plot HD data onto the tree
% HD_data = split_values_per_voxel(S, obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices)
% obj.ref.plot_value_tree(HD_data,'','','','',124,'','lines');
% R = unique(S);
% colorbar('Ticks',R);
% caxis([nanmin(R)-0.5, nanmax(R)+0.5])
% colormap(lines(numel(unique(S))));
% 
% figure(125);clf();
% for community = R'
%     plot(nanmean(obj.extracted_traces_conc(:,S == community),2)); hold on;                
% end
%             
for beh = {'peaks','peaks_~BodyCam_L_whisker','groups_peaks_~encoder','groups_peaks_encoder'}
    beh = beh{1}
    t = obj.get_tp_for_condition(beh);
    cross_corr  = corrcoef(data(t,:) - med(t));      
    to_keep = find(~all(isnan(cross_corr),1));
    cross_corr = cross_corr(to_keep, to_keep);
    eva = evalclusters(cross_corr,'linkage','silhouette','KList',1:30);
    cluster_idx = clusterdata(cross_corr,'Linkage', 'ward', 'MAXCLUST', eva.OptimalK);unique(cluster_idx)
    [a,b] = sort(cluster_idx);
    figure(1008);clf();imagesc(cross_corr(b,b))
    obj.ref.plot_value_tree(split_values_per_voxel(cluster_idx, obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), to_keep),'','','','',124,'','jet');
    saveas(figure(124), ['tree_HD_',beh,'_Nclust=',num2str(eva.OptimalK),'.png'])
    saveas(figure(1008), ['crosscorr_HD_',beh,'_Nclust=',num2str(eva.OptimalK),'.png'])
end

