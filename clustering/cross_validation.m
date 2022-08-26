%%%% Cross-validation of 'withheld' and 'test' data

%% Load OBJECT when already correct
% load('C:\Users\vanto\Documents\MATLAB\extracted_arboreal_scans 2\arboreal_scans_thin_mask.mat')
obj = arboreal_scan_experiment('');

%% Quick processing to have proper event detection. This is required for filtering based on activity
rescaling_method = 'peaks_trials';      % does this go here?
obj.prepare_binning({'depth',50});
obj.find_events();
obj.rescale_traces();
obj.set_median_traces();
obj.compute_similarity();
close all

%%
close all
clearvars -except obj

%% set various parameters number of PHATE loadings
N_Dim = 9;     % # of PHATE loadings
analysis_mode = 'space';
obj.filter_win = [10, 0];    % ensure no smoothing/no introduction of temporal correlation
conditions = {'active'};%,'~encoder_active'}
source_signal = obj.rescaled_traces-obj.global_median_raw;  % obj.extracted_traces_conc; obj.rescaled_traces - obj.global_median_raw;
   
%% number of shuffles/permutations
iter = 1; % 100, 500 or 1000? the latter 2 take a while, use 1 if just doing the control/withheld case

% pre-allocate for withheld group
Y_PHATE_3D_to_cross_corr_withheld = cell(1,iter);
var_PHATE_loadings_withheld = zeros(N_Dim,iter);
cumsum_PHATE_withheld = zeros(N_Dim,iter);
range_vec_all_withheld = zeros(N_Dim,iter);
cumsum_range_PHATE_loadings_withheld = zeros(N_Dim,iter);

% pre-allocate for test group
Y_PHATE_3D_to_cross_corr_test = cell(1,iter);
var_PHATE_loadings_test = zeros(N_Dim,iter);
cumsum_PHATE_test = zeros(N_Dim,iter);
range_vec_all_test = zeros(N_Dim,iter);
cumsum_range_PHATE_loadings_test = zeros(N_Dim,iter);


%% Don't keep any Inf values, if any (only happens in HD case, occasionally)
source_signal(isinf(source_signal)) = NaN;
Fig_count = 1000;

% Flag ROIs that have NaN values at one point as they may mess up later computations
bad_ROI_list                    = find(any(isnan(source_signal),1));
bad_ROI_list(bad_ROI_list > obj.n_ROIs) = [];
signal_indices                  = true(1, obj.n_ROIs); %% ROIs or voxels, depending on the data source
signal_indices(bad_ROI_list)    = false;
signal_indices                  = find(signal_indices);
    
%%
for iteration = 1:iter       
    %% Now, for each condition, classify and cluster the data
    for el = 1:numel(conditions)
        %% Get timpoints for current behavioural condition
        tp = obj.get_tp_for_condition(conditions{el});
        %     tp = diff(tp);
        %     tp(tp < 0) = 0;
        %     tp = [0, tp];
        
        %% Get signal
        current_signal = double(source_signal(:, 1:obj.n_ROIs));
        
        %% Filter out unrequired timepoints
        current_signal(~tp, :) = NaN;
        
        %% Filter out bad ROIs/voxels
        current_signal(:, bad_ROI_list) = NaN;
        
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
        
        %% Randomly sample points for cross-validation (e.g. control/withheld = %50 of data vs. test = %50 of data)
        current_signal = current_signal';
        %prepare to split the data
        length_data = size(current_signal);
        length_data = length_data(1,1);
        split_point = floor(length_data/2);                        % can change the fraction of the data split 70/30
        T_shuffled = randperm(length_data);
        first_half_T_shuffled = T_shuffled(1:split_point);
        num_col = size(current_signal);
        num_col = num_col(1,2);
        current_signal_for_PHATE_withheld = zeros(split_point,num_col);
        for ss = 1:num_col
            trace_to_shuffle = current_signal(:,ss);
            current_signal_for_PHATE_withheld(:,ss) = trace_to_shuffle(first_half_T_shuffled);
        end
        current_signal_for_PHATE_withheld = current_signal_for_PHATE_withheld';
        % use the remaining 50% of the data 
        second_half_T_shuffled = T_shuffled(split_point+1:end);
        if mod(length_data,2) == 0
            split_point = split_point;
        else split_point = split_point+1;       %% update this when splitting data in other ways (e.g. 70/30)
        end
        current_signal_for_PHATE_test = zeros(split_point,num_col);
        for ss = 1:num_col
            trace_to_shuffle = current_signal(:,ss);
            current_signal_for_PHATE_test(:,ss) = trace_to_shuffle(second_half_T_shuffled);
        end
        current_signal_for_PHATE_test = current_signal_for_PHATE_test';
                
         %% Randomly sample events (active time points) for cross-validation (e.g. control/withheld = %~50 of data vs. test = %~50 of data)
        current_signal = current_signal';
        %prepare to split the data
        length_data = size(current_signal);
        length_data = length_data(1,1);
        split_point = floor(length_data/2);                        % can change the fraction of the data split 70/30
        T_shuffled = randperm(length_data);
        first_half_T_shuffled = T_shuffled(1:split_point);
        num_col = size(current_signal);
        num_col = num_col(1,2);
        current_signal_for_PHATE_withheld = zeros(split_point,num_col);
        for ss = 1:num_col
            trace_to_shuffle = current_signal(:,ss);
            current_signal_for_PHATE_withheld(:,ss) = trace_to_shuffle(first_half_T_shuffled);
        end
        current_signal_for_PHATE_withheld = current_signal_for_PHATE_withheld';
        % use the remaining 50% of the data 
        second_half_T_shuffled = T_shuffled(split_point+1:end);
        if mod(length_data,2) == 0
            split_point = split_point;
        else split_point = split_point+1;       %% update this when splitting data in other ways (e.g. 70/30)
        end
        current_signal_for_PHATE_test = zeros(split_point,num_col);
        for ss = 1:num_col
            trace_to_shuffle = current_signal(:,ss);
            current_signal_for_PHATE_test(:,ss) = trace_to_shuffle(second_half_T_shuffled);
        end
        current_signal_for_PHATE_test = current_signal_for_PHATE_test';
        
       %% PHATE 3D
        Y_PHATE_3D_withheld = phate(current_signal_for_PHATE_withheld, 'ndim', N_Dim, 't', []);
        close(gcf);
%         figure(Fig_count + 3000);clf(); title('Phate first 3 dimensions scatter plot')
%         scatter3(Y_PHATE_3D_withheld(:,1), Y_PHATE_3D_withheld(:,2), Y_PHATE_3D_withheld(:,3), 30); hold on;
        
        Y_PHATE_3D_test = phate(current_signal_for_PHATE_test, 'ndim', N_Dim, 't', []);
        close(gcf);
%         figure(Fig_count + 3000);clf(); title('Phate first 3 dimensions scatter plot')
%         scatter3(Y_PHATE_3D_test(:,1), Y_PHATE_3D_test(:,2), Y_PHATE_3D_test(:,3), 30); hold on;
        
        
       %% fill the pre-allocations for withheld          
         % to cross-correlate
        Y_PHATE_3D_to_cross_corr_withheld{iteration} = Y_PHATE_3D_withheld;
        % determine var of each loading
        var_PHATE_loadings_withheld(:,iteration) = var(Y_PHATE_3D_withheld)';
        % determine cumsum of magnitude of each loading         
        cumsum_PHATE_withheld(:,iteration) = cumsum(mean(abs(Y_PHATE_3D_withheld),1))';
        % determine range of each loading
        range_vec = zeros(1,N_Dim);
        for jj = 1:N_Dim
            range_vec(jj) = range(Y_PHATE_3D_withheld(:,jj));
        end
        range_vec_all_withheld(:,iteration) = range_vec';
        % determine cumsum of range of each loading
        cumsum_range_PHATE_loadings_withheld(:,iteration) = cumsum(range(Y_PHATE_3D_withheld))';
        
        % fill the pre-allocations for test     
        % to cross-correlate
        Y_PHATE_3D_to_cross_corr_test{iteration} = Y_PHATE_3D_test;
        % determine var of each loading
        var_PHATE_loadings_test(:,iteration) = var(Y_PHATE_3D_test)';
        % determine cumsum of magnitude of each loading         
        cumsum_PHATE_test(:,iteration) = cumsum(mean(abs(Y_PHATE_3D_test),1))';
        % determine range of each loading
        range_vec = zeros(1,N_Dim);
        for jj = 1:N_Dim
            range_vec(jj) = range(Y_PHATE_3D_test(:,jj));
        end
        range_vec_all_test(:,iteration) = range_vec';
        % determine cumsum of range of each loading
        cumsum_range_PHATE_loadings_test(:,iteration) = cumsum(range(Y_PHATE_3D_test))';
                    
        %% Display Phates on tree for withheld
%         n_row = floor(sqrt(size(Y_PHATE_3D_withheld, 2))) + 1;
%         n_col = ceil(sqrt(size(Y_PHATE_3D_withheld, 2)));
%         figure(Fig_count + 2000);clf()
%         for dim = 1:size(Y_PHATE_3D_withheld, 2)
%             sub = subplot(n_row,n_col,dim);
%             if obj.use_hd_data
%                 obj.ref.plot_value_tree(split_values_per_voxel(Y_PHATE_3D_withheld(:,dim), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',sub,'curved','viridis');
%             else
%                 obj.ref.plot_value_tree(Y_PHATE_3D_withheld(:,dim), find(valid_ROIs),'',['phate #',num2str(dim),' Loadings (per ribbon)'],'',sub,'curved','viridis');
%             end
%         end
%         
%         figure();plot(Y_PHATE_3D_withheld)
%         figure();imagesc(Y_PHATE_3D_withheld);colorbar;%caxis([-1 1]);
        
         %% Display Phates on tree for test
%         n_row = floor(sqrt(size(Y_PHATE_3D_test, 2))) + 1;
%         n_col = ceil(sqrt(size(Y_PHATE_3D_test, 2)));
%         figure(Fig_count + 2000);clf()
%         for dim = 1:size(Y_PHATE_3D_test, 2)
%             sub = subplot(n_row,n_col,dim);
%             if obj.use_hd_data
%                 obj.ref.plot_value_tree(split_values_per_voxel(Y_PHATE_3D_test(:,dim), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',sub,'curved','viridis');
%             else
%                 obj.ref.plot_value_tree(Y_PHATE_3D_test(:,dim), find(valid_ROIs),'',['phate #',num2str(dim),' Loadings (per ribbon)'],'',sub,'curved','viridis');
%             end
%         end
%         
%         figure();plot(Y_PHATE_3D_test)
%         figure();imagesc(Y_PHATE_3D_test);colorbar;%caxis([-1 1]);
        
    end
end



%% plot the range of each loading in histograms for withheld & test groups

%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!

% withheld:

% if extracting data from an open figure, drag/drop the figre into Cmnd Line
% or use uiopen('C:\Users\...\filepath\filename',1), then run this script:
% Extract_data_from_open_figure

figure(44);subplot(4,3,1);histogram(range_vec_all_withheld(1,:),20); title('Distr of range for withheld data for Loading 1')
subplot(4,3,2);histogram(range_vec_all_withheld(2,:),20); title('Loading 2')
subplot(4,3,3);histogram(range_vec_all_withheld(3,:),20); title('Loading 3')
subplot(4,3,4);histogram(range_vec_all_withheld(4,:),20); title('Loading 4')
subplot(4,3,5);histogram(range_vec_all_withheld(5,:),20); title('Loading 5')
subplot(4,3,6);histogram(range_vec_all_withheld(6,:),20); title('Loading 6')
subplot(4,3,7);histogram(range_vec_all_withheld(7,:),20); title('Loading 7')
subplot(4,3,8);histogram(range_vec_all_withheld(8,:),20); title('Loading 8')
subplot(4,3,9);histogram(range_vec_all_withheld(9,:),20); title('Loading 9')
%subplot(4,3,10);histogram(range_vec_all_withheld(10,:),20); title('Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_0.33075s_100x_range_distr';
filepathc = append(filepathc,filenamec);
%saveas(figure(44),filepathc, 'png');
filepathd = append(filepathc, '.fig');
%savefig(figure(44),filepathd,'compact');

%test group :

figure(45);subplot(4,3,1);histogram(range_vec_all_test(1,:),20); title('Distr ranges test for Loading 1') % xlim([0, 6]); 
subplot(4,3,2);histogram(range_vec_all_test(2,:),20); title('Loading 2')
subplot(4,3,3);histogram(range_vec_all_test(3,:),20); title('Loading 3')
subplot(4,3,4);histogram(range_vec_all_test(4,:),20); title('Loading 4')
subplot(4,3,5);histogram(range_vec_all_test(5,:),20); title('Loading 5')
subplot(4,3,6);histogram(range_vec_all_test(6,:),20); title('Loading 6')
subplot(4,3,7);histogram(range_vec_all_test(7,:),20); title('Loading 7')
subplot(4,3,8);histogram(range_vec_all_test(8,:),20); title('Loading 8')
subplot(4,3,9);histogram(range_vec_all_test(9,:),20); title('Loading 9')
%subplot(4,3,10);histogram(range_vec_all_test(10,:),20); title('Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_0.33075s_100x_range_distr';
filepathc = append(filepathc,filenamec);
%saveas(figure(45),filepathc, 'png');
filepathd = append(filepathc, '.fig');
%savefig(figure(45),filepathd,'compact');



%% plot the summary data for withheld and test groups: individual (gray) and mean (red)

%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!

% withheld:

figure(71);subplot(2,2,1);plot(var_PHATE_loadings_withheld, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(var_PHATE_loadings_withheld,2),'color', [1 0 0], 'linewidth', 1.5); title('Var of each Loading withheld')
subplot(2,2,2);plot(cumsum_PHATE_withheld, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_PHATE_withheld,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each mean(abs(Loading)) withheld')
subplot(2,2,3);plot(range_vec_all_withheld, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(range_vec_all_withheld,2),'color', [1 0 0], 'linewidth', 1.5); title('Range for each Loading withheld')
subplot(2,2,4);plot(cumsum_range_PHATE_loadings_withheld, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_range_PHATE_loadings_withheld,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each range(Loading) withheld')

filepath = obj.source_folder;
filename = 'randperm_all_traces_same_100x_cumsum';
filepath = append(filepath,filename);
%saveas(figure(71),filepath, 'png');
filepatha = append(filepath, '.fig');
%savefig(figure(71),filepatha,'compact');

%test 

figure(72);subplot(2,2,1);plot(var_PHATE_loadings_test, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(var_PHATE_loadings_test,2),'color', [1 0 0], 'linewidth', 1.5); title('Var of each Loading_test')
subplot(2,2,2);plot(cumsum_PHATE_test, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_PHATE_test,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each mean(abs(Loading))_test')
subplot(2,2,3);plot(range_vec_all_test, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(range_vec_all_test,2),'color', [1 0 0], 'linewidth', 1.5); title('Range for each Loading_test')
subplot(2,2,4);plot(cumsum_range_PHATE_loadings_test, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_range_PHATE_loadings_test,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each range(Loading)_test')

filepath = obj.source_folder;
filename = 'randperm_all_traces_same_100x_cumsum';
filepath = append(filepath,filename);
%saveas(figure(72),filepath, 'png');
filepatha = append(filepath, '.fig');
%savefig(figure(72),filepatha,'compact');

%% prepare and run pwCC for each loading for the withheld group
col_val = size(Y_PHATE_3D_withheld);
col_val = col_val(1,1);
to_pwCC_matrix_L1_withheld = zeros(col_val,iter);
to_pwCC_matrix_L2_withheld = zeros(col_val,iter);
to_pwCC_matrix_L3_withheld = zeros(col_val,iter);
to_pwCC_matrix_L4_withheld = zeros(col_val,iter);
to_pwCC_matrix_L5_withheld = zeros(col_val,iter);
to_pwCC_matrix_L6_withheld = zeros(col_val,iter);
to_pwCC_matrix_L7_withheld = zeros(col_val,iter);
to_pwCC_matrix_L8_withheld = zeros(col_val,iter);
to_pwCC_matrix_L9_withheld = zeros(col_val,iter);
%to_pwCC_matrix_L10_withheld = zeros(col_val,iter);

for rr = 1:iter
    to_pwCC_matrix_L1_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,1)];
end
pwCCL1_withheld = corrcoef(to_pwCC_matrix_L1_withheld);

for rr = 1:iter
    to_pwCC_matrix_L2_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,2)];
end
pwCCL2_withheld = corrcoef(to_pwCC_matrix_L2_withheld);

for rr = 1:iter
    to_pwCC_matrix_L3_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,3)];
end
pwCCL3_withheld = corrcoef(to_pwCC_matrix_L3_withheld);

for rr = 1:iter
    to_pwCC_matrix_L4_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,4)];
end
pwCCL4_withheld = corrcoef(to_pwCC_matrix_L4_withheld);

for rr = 1:iter
    to_pwCC_matrix_L5_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,5)];
end
pwCCL5_withheld = corrcoef(to_pwCC_matrix_L5_withheld);

for rr = 1:iter
    to_pwCC_matrix_L6_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,6)];
end
pwCCL6_withheld = corrcoef(to_pwCC_matrix_L6_withheld);

for rr = 1:iter
    to_pwCC_matrix_L7_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,7)];
end
pwCCL7_withheld = corrcoef(to_pwCC_matrix_L7_withheld);

for rr = 1:iter
    to_pwCC_matrix_L8_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,8)];
end
pwCCL8_withheld = corrcoef(to_pwCC_matrix_L8_withheld);

for rr = 1:iter
    to_pwCC_matrix_L9_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,9)];
end
pwCCL9_withheld = corrcoef(to_pwCC_matrix_L9_withheld);

% for rr = 1:iter
%     to_pwCC_matrix_L10_withheld(:,rr) = [Y_PHATE_3D_to_cross_corr_withheld{rr}(:,10)];
% end
% pwCCL10_withheld = corrcoef(to_pwCC_matrix_L10_withheld);


% prepare and run pwCC for each loading for the test group
col_val = size(current_signal);
col_val = col_val(1,2);
to_pwCC_matrix_L1_test = zeros(col_val,iter);
to_pwCC_matrix_L2_test = zeros(col_val,iter);
to_pwCC_matrix_L3_test = zeros(col_val,iter);
to_pwCC_matrix_L4_test = zeros(col_val,iter);
to_pwCC_matrix_L5_test = zeros(col_val,iter);
to_pwCC_matrix_L6_test = zeros(col_val,iter);
to_pwCC_matrix_L7_test = zeros(col_val,iter);
to_pwCC_matrix_L8_test = zeros(col_val,iter);
to_pwCC_matrix_L9_test = zeros(col_val,iter);
%to_pwCC_matrix_L10_test = zeros(col_val,iter);

for rr = 1:iter
    to_pwCC_matrix_L1_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,1)];
end
pwCCL1_test = corrcoef(to_pwCC_matrix_L1_test);

for rr = 1:iter
    to_pwCC_matrix_L2_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,2)];
end
pwCCL2_test = corrcoef(to_pwCC_matrix_L2_test);

for rr = 1:iter
    to_pwCC_matrix_L3_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,3)];
end
pwCCL3_test = corrcoef(to_pwCC_matrix_L3_test);

for rr = 1:iter
    to_pwCC_matrix_L4_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,4)];
end
pwCCL4_test = corrcoef(to_pwCC_matrix_L4_test);

for rr = 1:iter
    to_pwCC_matrix_L5_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,5)];
end
pwCCL5_test = corrcoef(to_pwCC_matrix_L5_test);

for rr = 1:iter
    to_pwCC_matrix_L6_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,6)];
end
pwCCL6_test = corrcoef(to_pwCC_matrix_L6_test);

for rr = 1:iter
    to_pwCC_matrix_L7_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,7)];
end
pwCCL7_test = corrcoef(to_pwCC_matrix_L7_test);

for rr = 1:iter
    to_pwCC_matrix_L8_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,8)];
end
pwCCL8_test = corrcoef(to_pwCC_matrix_L8_test);

for rr = 1:iter
    to_pwCC_matrix_L9_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,9)];
end
pwCCL9_test = corrcoef(to_pwCC_matrix_L9_test);

% for rr = 1:iter
%     to_pwCC_matrix_L10_test(:,rr) = [Y_PHATE_3D_to_cross_corr_test{rr}(:,10)];
% end
% pwCCL10_test = corrcoef(to_pwCC_matrix_L10_test);




%% plot the pwCCs for each loading for the withheld and test data
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!

%withheld

figure(73);subplot(4,3,1);imagesc(pwCCL1_withheld);colorbar; caxis([-1 1]);title('pwCCs withheld for Loading 1')
subplot(4,3,2);imagesc(pwCCL2_withheld); colorbar; caxis([-1 1]); title('Loading 2')
subplot(4,3,3);imagesc(pwCCL3_withheld); colorbar; caxis([-1 1]); title('Loading 3')
subplot(4,3,4);imagesc(pwCCL4_withheld); colorbar; caxis([-1 1]); title('Loading 4')
subplot(4,3,5);imagesc(pwCCL5_withheld); colorbar; caxis([-1 1]); title('Loading 5')
subplot(4,3,6);imagesc(pwCCL6_withheld); colorbar; caxis([-1 1]); title('Loading 6')
subplot(4,3,7);imagesc(pwCCL7_withheld); colorbar; caxis([-1 1]); title('Loading 7')
subplot(4,3,8);imagesc(pwCCL8_withheld); colorbar; caxis([-1 1]); title('Loading 8')
subplot(4,3,9);imagesc(pwCCL9_withheld); colorbar; caxis([-1 1]); title('Loading 9')
%subplot(4,3,10);imagesc(pwCCL10_withheld); colorbar; caxis([-1 1]); title('Loading 10')

filepathb = obj.source_folder;
filenameb = 'blockperm_all_traces_0.33075s_100x_pwCC';
filepathb = append(filepathb,filenameb);
%saveas(figure(73),filepathb, 'png');
filepathc = append(filepathb, '.fig');
%savefig(figure(73),filepathc,'compact');

%get values of the pwCCs and put them in a col_vector for visualizing distributions
col_vect_pwCCL1_withheld = pwCCL1_withheld(~tril(ones(size(pwCCL1_withheld))));
col_vect_pwCCL2_withheld = pwCCL2_withheld(~tril(ones(size(pwCCL2_withheld))));
col_vect_pwCCL3_withheld = pwCCL3_withheld(~tril(ones(size(pwCCL3_withheld))));
col_vect_pwCCL4_withheld = pwCCL4_withheld(~tril(ones(size(pwCCL4_withheld))));
col_vect_pwCCL5_withheld = pwCCL5_withheld(~tril(ones(size(pwCCL5_withheld))));
col_vect_pwCCL6_withheld = pwCCL6_withheld(~tril(ones(size(pwCCL6_withheld))));
col_vect_pwCCL7_withheld = pwCCL7_withheld(~tril(ones(size(pwCCL7_withheld))));
col_vect_pwCCL8_withheld = pwCCL8_withheld(~tril(ones(size(pwCCL8_withheld))));
col_vect_pwCCL9_withheld = pwCCL9_withheld(~tril(ones(size(pwCCL9_withheld))));
%col_vect_pwCCL10_withheld = pwCCL10_withheld(~tril(ones(size(pwCCL10_withheld))));

figure(74);subplot(4,3,1);histogram(col_vect_pwCCL1_withheld,20); xlim([-1, 1.05]); title('Distr of pwCCs withheld for Loading 1')
subplot(4,3,2);histogram(col_vect_pwCCL2_withheld,20); xlim([-1, 1.05]); title('Loading 2')
subplot(4,3,3);histogram(col_vect_pwCCL3_withheld,20); xlim([-1, 1.05]); title('Loading 3')
subplot(4,3,4);histogram(col_vect_pwCCL4_withheld,20); xlim([-1, 1]); title('Loading 4')
subplot(4,3,5);histogram(col_vect_pwCCL5_withheld,20); xlim([-1, 1]); title('Loading 5')
subplot(4,3,6);histogram(col_vect_pwCCL6_withheld,20); xlim([-1, 1]); title('Loading 6')
subplot(4,3,7);histogram(col_vect_pwCCL7_withheld,20); xlim([-1, 1]); title('Loading 7')
subplot(4,3,8);histogram(col_vect_pwCCL8_withheld,20); xlim([-1, 1]); title('Loading 8')
subplot(4,3,9);histogram(col_vect_pwCCL9_withheld,20); xlim([-1, 1]); title('Loading 9')
%subplot(4,3,10);histogram(col_vect_pwCCL10_withheld,20); xlim([-1, 1]); title('Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_0.33075s_100x_distr';
filepathc = append(filepathc,filenamec);
%saveas(figure(74),filepathc, 'png');
filepathd = append(filepathc, '.fig');
%savefig(figure(74),filepathd,'compact');

%test
% plot the pwCCs for each loading
figure(75);subplot(4,3,1);imagesc(pwCCL1_test);colorbar; caxis([-1 1]);title('pwCCs_test for Loading 1')
subplot(4,3,2);imagesc(pwCCL2_test); colorbar; caxis([-1 1]); title('Loading 2')
subplot(4,3,3);imagesc(pwCCL3_test); colorbar; caxis([-1 1]); title('Loading 3')
subplot(4,3,4);imagesc(pwCCL4_test); colorbar; caxis([-1 1]); title('Loading 4')
subplot(4,3,5);imagesc(pwCCL5_test); colorbar; caxis([-1 1]); title('Loading 5')
subplot(4,3,6);imagesc(pwCCL6_test); colorbar; caxis([-1 1]); title('Loading 6')
subplot(4,3,7);imagesc(pwCCL7_test); colorbar; caxis([-1 1]); title('Loading 7')
subplot(4,3,8);imagesc(pwCCL8_test); colorbar; caxis([-1 1]); title('Loading 8')
subplot(4,3,9);imagesc(pwCCL9_test); colorbar; caxis([-1 1]); title('Loading 9')
%subplot(4,3,10);imagesc(pwCCL10_test); colorbar; caxis([-1 1]); title('Loading 10')

filepathb = obj.source_folder;
filenameb = 'blockperm_all_traces_0.33075s_100x_pwCC';
filepathb = append(filepathb,filenameb);
%saveas(figure(75),filepathb, 'png');
filepathc = append(filepathb, '.fig');
%savefig(figure(75),filepathc,'compact');

%get values of the pwCCs and put them in a col_vector for visualizing distributions
col_vect_pwCCL1_test = pwCCL1_test(~tril(ones(size(pwCCL1_test))));
col_vect_pwCCL2_test = pwCCL2_test(~tril(ones(size(pwCCL2_test))));
col_vect_pwCCL3_test = pwCCL3_test(~tril(ones(size(pwCCL3_test))));
col_vect_pwCCL4_test = pwCCL4_test(~tril(ones(size(pwCCL4_test))));
col_vect_pwCCL5_test = pwCCL5_test(~tril(ones(size(pwCCL5_test))));
col_vect_pwCCL6_test = pwCCL6_test(~tril(ones(size(pwCCL6_test))));
col_vect_pwCCL7_test = pwCCL7_test(~tril(ones(size(pwCCL7_test))));
col_vect_pwCCL8_test = pwCCL8_test(~tril(ones(size(pwCCL8_test))));
col_vect_pwCCL9_test = pwCCL9_test(~tril(ones(size(pwCCL9_test))));
%col_vect_pwCCL10_test = pwCCL10_test(~tril(ones(size(pwCCL10_test))));

figure(76);subplot(4,3,1);histogram(col_vect_pwCCL1_test,20); xlim([-1, 1.05]); title('Distr of pwCCs_test for Loading 1')
subplot(4,3,2);histogram(col_vect_pwCCL2_test,20); xlim([-1, 1.05]); title('Loading 2')
subplot(4,3,3);histogram(col_vect_pwCCL3_test,20); xlim([-1, 1.05]); title('oading 3')
subplot(4,3,4);histogram(col_vect_pwCCL4_test,20); xlim([-1, 1]); title('Loading 4')
subplot(4,3,5);histogram(col_vect_pwCCL5_test,20); xlim([-1, 1]); title('Loading 5')
subplot(4,3,6);histogram(col_vect_pwCCL6_test,20); xlim([-1, 1]); title('Loading 6')
subplot(4,3,7);histogram(col_vect_pwCCL7_test,20); xlim([-1, 1]); title('Loading 7')
subplot(4,3,8);histogram(col_vect_pwCCL8_test,20); xlim([-1, 1]); title('Loading 8')
subplot(4,3,9);histogram(col_vect_pwCCL9_test,20); xlim([-1, 1]); title('Loading 9')
%subplot(4,3,10);histogram(col_vect_pwCCL10_test,20); xlim([-1, 1]); title('Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_0.33075s_100x_distr';
filepathc = append(filepathc,filenamec);
%saveas(figure(76),filepathc, 'png');
filepathd = append(filepathc, '.fig');
%savefig(figure(76),filepathd,'compact');



%% generate DBSCAN clustering for withheld
% estimate epsilon for subsequent (h)DBScan clustering      
kD = pdist2(Y_PHATE_3D_withheld,Y_PHATE_3D_withheld,'euc','Smallest',5);
kd_sorted = sort(kD(:))';
slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
[~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
epsilon = kd_sorted(minloc);
%% estimate epsilon from 'knee' ## Uncomment second line to see how clustering evolves with epsilon
suggested = test_epsilon(obj, Y_PHATE_3D_withheld, current_signal, all_ROIs, valid_ROIs); %open test_epsilon and change rendering to true to see it
%epsilon = 1;
%suggested = 2;
%% Show clusters, push on tree, show traces
cluster_idx{el} = phate_figure(obj, Y_PHATE_3D_withheld, mean([epsilon, suggested]), source_signal(tp,:), Fig_count, signal_indices);
hold on;sgtitle(['Cluster for condition : ',strrep(conditions{el},'_','\_')])


%% generate DBSCAN clustering for test
% estimate epsilon for subsequent (h)DBScan clustering      
kD = pdist2(Y_PHATE_3D_test,Y_PHATE_3D_test,'euc','Smallest',5);
kd_sorted = sort(kD(:))';
slope = (kd_sorted(end) - kd_sorted(1)) / numel(kd_sorted);
[~, minloc] = min(kd_sorted - ((1:numel(kd_sorted)) * slope));
epsilon = kd_sorted(minloc);
%% estimate epsilon from 'knee' ## Uncomment second line to see how clustering evolves with epsilon
suggested = test_epsilon(obj, Y_PHATE_3D_test, current_signal, all_ROIs, valid_ROIs); %open test_epsilon and change rendering to true to see it
%epsilon = 1;
%suggested = 2;
%% Show clusters, push on tree, show traces
cluster_idx{el} = phate_figure(obj, Y_PHATE_3D_test, mean([epsilon, suggested]), source_signal(tp,:), Fig_count, signal_indices);
hold on;sgtitle(['Cluster for condition : ',strrep(conditions{el},'_','\_')])