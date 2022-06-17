%% Convert list of voxel values into a cell array of nROIs with their correct voxels values. Use this to plot a HD tree 
% obj.ref.plot_value_tree(values, '','','values per voxel','','','curved');


%% Example to run
% all_values = obj.extracted_traces_conc;
% for t = 1:size(all_values, 1)
%     v = nanmean(all_values(t, :),1);
%     formatted_values = split_values_per_voxel(v, obj.ref.header.res_list(:,1));
%     if t == 1
%         f = obj.ref.plot_value_tree(formatted_values, '','','cluster tree (one value per voxel)','','','curved');
%         caxis([1,9]);        
%         f  = f.Parent;
%     else
%         obj.ref.plot_value_tree(formatted_values, '','','cluster tree (one value per voxel)','',f,'curved');
%         drawnow
%     end
% end

function formatted_values = split_values_per_voxel(input_values, voxels_per_ROI, voxel_indexes)
    if (nargin < 3 || isempty(voxel_indexes)) && sum(voxels_per_ROI) == size(input_values, 2)
        voxel_indexes = 1:size(input_values, 2);
    elseif nargin < 3 || isempty(voxel_indexes)
        error('if the number of values don t match the number of voxel, you must indicate the voxel indexes!')
    end
        
    % typically, voxels_per_ROI = obj.ref.header.res_list(:,1)
    n_ROIs              = numel(voxels_per_ROI);
    formatted_values              = cell(1, n_ROIs);
    start_vxl_per_ROI   = [cumsum(voxels_per_ROI) - voxels_per_ROI(1) + 1; sum(voxels_per_ROI)+1];
    for ROI_idx = 1:n_ROIs
        formatted_values{ROI_idx}     = NaN(1, voxels_per_ROI(ROI_idx));   
        valid_ROIs          = ismember(voxel_indexes, start_vxl_per_ROI(ROI_idx):(start_vxl_per_ROI(ROI_idx+1)-1));
        v                   = input_values(valid_ROIs);
        valid_vxl_in_ROI    = ismember(start_vxl_per_ROI(ROI_idx):(start_vxl_per_ROI(ROI_idx+1)-1) , voxel_indexes);
        formatted_values{ROI_idx}(valid_vxl_in_ROI)     = v;
    end
end