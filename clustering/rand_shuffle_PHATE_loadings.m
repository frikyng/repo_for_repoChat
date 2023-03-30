%%%% Shuffle/permute time points for all ROIs (in user-defined ways below), then run PHATE

%% Load OBJECT when already correct
% load('C:\Users\vanto\Documents\MATLAB\extracted_arboreal_scans 2\arboreal_scans_thin_mask.mat')
obj = arboreal_scan_experiment('');

%% Quick processing to have proper event detection. This is required for filtering based on activity
rescaling_method = 'peaks_trials';
obj.prepare_binning({'distance',50});
obj.find_events();
obj.rescale_traces();
obj.set_median_traces();
obj.compute_similarity();
close all

%%
close all
clearvars -except obj

%% set various parameters number of PHATE loadings
N_Dim = 10;     % # of PHATE loadings
analysis_mode = 'space';
obj.time_smoothing = [0, 0];    % ensure no smoothing/no introduction of temporal correlation
conditions = {''};%,'~encoder_active'}
source_signal = obj.rescaled_traces;

%% number of shuffles/permutations
iter = 10; % 100, 500 or 1000? the latter 2 take a while, use 1 if just doing the control/test case

%pre-allocate 
Y_PHATE_3D_to_cross_corr = cell(1,iter);
var_PHATE_loadings = zeros(N_Dim,iter);
cumsum_PHATE = zeros(N_Dim,iter);
range_vec_all = zeros(N_Dim,iter);
cumsum_range_PHATE_loadings = zeros(N_Dim,iter);

%% Don't keep any Inf values, if any (only happens in HD case, occasionally)
source_signal(isinf(source_signal)) = NaN;

Fig_count = 1000;

% Flag ROIs that have NaN vaues at one point as they may mess up later computations
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
        
        %% Fully randomized points using same number of time points as original data (i.e. 'true baseline' lacking any structure)
%         current_signal_for_PHATE = randn(size(current_signal'))';
%         num_col = size(current_signal_for_PHATE);
%         num_col = num_col(1,2);     % this is used later in the code and is here for consistency with the other shuffles
 
        %% Randomize/shuffle all time points for each ROI separately (i.e. try to destroy all correlations)
        current_signal = current_signal';
        current_signal_for_PHATE = zeros(size(current_signal));
        num_col = size(current_signal_for_PHATE);
        num_col = num_col(1,2);
        for ii = 1:num_col
            trace_to_shuffle = current_signal(:,ii);
            current_signal_for_PHATE(:,ii) = trace_to_shuffle(randperm(length(trace_to_shuffle)));
        end
        current_signal_for_PHATE = current_signal_for_PHATE';
 
        %% Randomize/shuffle every ROI together in the same way
%           current_signal = current_signal';
%           current_signal_for_PHATE = zeros(size(current_signal));
%           num_col = size(current_signal_for_PHATE);
%           num_col = num_col(1,2);
%           shuffled_trace = randperm(length(current_signal));
%           for ss = 1:num_col
%               trace_to_shuffle = current_signal(:,ss);
%               current_signal_for_PHATE(:,ss) = trace_to_shuffle(shuffled_trace);
%           end
%           current_signal_for_PHATE = current_signal_for_PHATE';
          
        %% Block shuffle time points for each ROI separately using Alex's custom function: block_shuffle_time.m
        % first, get a vector/list of block shuffled timepoints
        % determine T = timepoints, acq_rate, and block_s to shuffle; est block_s from tau decay of median event width for our cell of interest 
        current_signal = current_signal';
        T = size(current_signal);
        T = T(1);
        acquisition_rate = obj.ref.header.points_per_s;
        list_event_widths = vertcat(2,obj.event.peak_width{:});
        median_event_width = median(list_event_widths);
        block_s = median_event_width*.37; % this est. just for starters, try shorter and longer block_s
        block_s = block_s*2; % note the divisor if present (should be 1 if want to use estimated tau)

        % check the shuffling procedure
        % T_shuffled = block_shuffle_time(T,acquisition_rate,block_s)
        % figure();plot(1:T);hold on; plot(T_shuffled)

        % randomize/shuffle all blocked time points for each ROI separately (i.e. try to destroy all correlations)
        current_signal_for_PHATE = zeros(size(current_signal));
        num_col = size(current_signal_for_PHATE);
        num_col = num_col(1,2);
        for gg = 1:num_col
            T_shuffled = block_shuffle_time(T,acquisition_rate,block_s);
            trace_to_shuffle = current_signal(:,gg);
            current_signal_for_PHATE(:,gg) = trace_to_shuffle(T_shuffled);
        end
        current_signal_for_PHATE = current_signal_for_PHATE';
        
        %% PHATE 3D
        Y_PHATE_3D = phate(current_signal_for_PHATE, 'ndim', N_Dim, 't', []);
        close(gcf);
%         figure(Fig_count + 3000);clf(); title('Phate first 3 dimensions scatter plot')
%         scatter3(Y_PHATE_3D(:,1), Y_PHATE_3D(:,2), Y_PHATE_3D(:,3), 30); hold on;
        
        
        %% Display Phates on tree
        nn_row = floor(sqrt(size(Y_PHATE_3D, 2))) + 1;
        nn_col = ceil(sqrt(size(Y_PHATE_3D, 2)));
        figure(Fig_count + 2000);clf()
        for dim = 1:size(Y_PHATE_3D, 2)
            sub = subplot(nn_row,nn_col,dim);
            if obj.use_hd_data
                obj.ref.plot_value_tree(split_values_per_voxel(Y_PHATE_3D(:,dim), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',sub,'curved','viridis');
            else
                obj.ref.plot_value_tree(Y_PHATE_3D(:,dim), find(valid_ROIs),'',['phate #',num2str(dim),' Loadings (per ribbon)'],'',sub,'curved','viridis');
            end
        end
        
        figure();plot(Y_PHATE_3D)
        figure();imagesc(Y_PHATE_3D);caxis([-1 1]);colorbar;
        
        % to cross-correlate 
        Y_PHATE_3D_to_cross_corr{iteration} = Y_PHATE_3D;
        % determine var of each loading
        var_PHATE_loadings(:,iteration) = var(Y_PHATE_3D)';
        % determine cumsum of magnitude of each loading        
        cumsum_PHATE(:,iteration) = cumsum(mean(abs(Y_PHATE_3D),1))';
        % determine range value of each loading
        range_vec = zeros(1,N_Dim);
        for jj = 1:N_Dim
            range_vec(jj) = range(Y_PHATE_3D(:,jj));
        end
        range_vec_all(:,iteration) = range_vec';
        % determine cumsum of range value of each loading
        cumsum_range_PHATE_loadings(:,iteration) = cumsum(range(Y_PHATE_3D))';
        
    end
end

%% plot the range of each loading in histograms
% if extracting data from an open figure, drag/drop the figre into Cmnd Line
% or use uiopen('C:\Users\...\filepath\filename',1), then run this script:
% Extract_data_from_open_figure

figure(88);subplot(4,3,1);histogram(range_vec_all(1,:),20); xlim([0, 8]); title('Distr of ranges for Loading 1')
subplot(4,3,2);histogram(range_vec_all(2,:),20); xlim([0, 8]); title('Distr of ranges for Loading 2')
subplot(4,3,3);histogram(range_vec_all(3,:),20); xlim([0, 8]); title('Distr of ranges for Loading 3')
subplot(4,3,4);histogram(range_vec_all(4,:),20); xlim([0, 6]); title('Distr of ranges for Loading 4')
subplot(4,3,5);histogram(range_vec_all(5,:),20); xlim([0, 6]); title('Distr of ranges for Loading 5')
subplot(4,3,6);histogram(range_vec_all(6,:),20); xlim([0, 6]); title('Distr of ranges for Loading 6')
subplot(4,3,7);histogram(range_vec_all(7,:),20); xlim([0, 6]); title('Distr of ranges for Loading 7')
subplot(4,3,8);histogram(range_vec_all(8,:),20); xlim([0, 6]); title('Distr of ranges for Loading 8')
subplot(4,3,9);histogram(range_vec_all(9,:),20); xlim([0, 6]); title('Distr of ranges for Loading 9')
subplot(4,3,10);histogram(range_vec_all(10,:),20); xlim([0, 6]); title('Distr of ranges for Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_2.646s_100x_range_distr';
filepathc = append(filepathc,filenamec);
saveas(figure(88),filepathc, 'png');
filepathd = append(filepathc, '.fig');
savefig(figure(88),filepathd,'compact');

%% plot the summary data: individual (gray) and mean (red)
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
figure(89);subplot(2,2,1);plot(var_PHATE_loadings, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(var_PHATE_loadings,2),'color', [1 0 0], 'linewidth', 1.5); title('Var of each Loading')
subplot(2,2,2);plot(cumsum_PHATE, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_PHATE,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each mean(abs(Loading))')
subplot(2,2,3);plot(range_vec_all, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(range_vec_all,2),'color', [1 0 0], 'linewidth', 1.5); title('Range for each Loading')
subplot(2,2,4);plot(cumsum_range_PHATE_loadings, 'color', [.5 .5 .5], 'linewidth', 1.0); hold on; plot(mean(cumsum_range_PHATE_loadings,2),'color', [1 0 0], 'linewidth', 1.5); title('Cum Sum of each range(Loading)')

filepath = obj.source_folder;
filename = 'blockperm_all_traces_2.646s_100x_cumsum';
filepath = append(filepath,filename);
saveas(figure(89),filepath, 'png');
filepatha = append(filepath, '.fig');
savefig(figure(89),filepatha,'compact');

%% prepare and run pwCC for each loading 
col_val = num_col;
to_pwCC_matrix_L1 = zeros(col_val,iter);
to_pwCC_matrix_L2 = zeros(col_val,iter);
to_pwCC_matrix_L3 = zeros(col_val,iter);
to_pwCC_matrix_L4 = zeros(col_val,iter);
to_pwCC_matrix_L5 = zeros(col_val,iter);
to_pwCC_matrix_L6 = zeros(col_val,iter);
to_pwCC_matrix_L7 = zeros(col_val,iter);
to_pwCC_matrix_L8 = zeros(col_val,iter);
to_pwCC_matrix_L9 = zeros(col_val,iter);
to_pwCC_matrix_L10 = zeros(col_val,iter);

for rr = 1:iter
    to_pwCC_matrix_L1(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,1)];
end

pwCCL1 = corrcoef(to_pwCC_matrix_L1);


for rr = 1:iter
    to_pwCC_matrix_L2(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,2)];
end

pwCCL2 = corrcoef(to_pwCC_matrix_L2);


for rr = 1:iter
    to_pwCC_matrix_L3(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,3)];
end

pwCCL3 = corrcoef(to_pwCC_matrix_L3);


for rr = 1:iter
    to_pwCC_matrix_L4(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,4)];
end

pwCCL4 = corrcoef(to_pwCC_matrix_L4);


for rr = 1:iter
    to_pwCC_matrix_L5(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,5)];
end

pwCCL5 = corrcoef(to_pwCC_matrix_L5);


for rr = 1:iter
    to_pwCC_matrix_L6(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,6)];
end

pwCCL6 = corrcoef(to_pwCC_matrix_L6);


for rr = 1:iter
    to_pwCC_matrix_L7(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,7)];
end

pwCCL7 = corrcoef(to_pwCC_matrix_L7);


for rr = 1:iter
    to_pwCC_matrix_L8(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,8)];
end

pwCCL8 = corrcoef(to_pwCC_matrix_L8);


for rr = 1:iter
    to_pwCC_matrix_L9(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,9)];
end

pwCCL9 = corrcoef(to_pwCC_matrix_L9);


for rr = 1:iter
    to_pwCC_matrix_L10(:,rr) = [Y_PHATE_3D_to_cross_corr{rr}(:,10)];
end

pwCCL10 = corrcoef(to_pwCC_matrix_L10);

%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!
%% CHANGE THE FILE NAME TO SAVE
%% !!!!!!!!!!!!!!!!!!!!!!!!!

% plot the pwCCs for each loading
figure(90);subplot(4,3,1);imagesc(pwCCL1);colorbar; caxis([-1 1]);title('pwCCs for Loading 1')
subplot(4,3,2);imagesc(pwCCL2); colorbar; caxis([-1 1]); title('pwCCs for Loading 2')
subplot(4,3,3);imagesc(pwCCL3); colorbar; caxis([-1 1]); title('pwCCs for Loading 3')
subplot(4,3,4);imagesc(pwCCL4); colorbar; caxis([-1 1]); title('pwCCs for Loading 4')
subplot(4,3,5);imagesc(pwCCL5); colorbar; caxis([-1 1]); title('pwCCs for Loading 5')
subplot(4,3,6);imagesc(pwCCL6); colorbar; caxis([-1 1]); title('pwCCs for Loading 6')
subplot(4,3,7);imagesc(pwCCL7); colorbar; caxis([-1 1]); title('pwCCs for Loading 7')
subplot(4,3,8);imagesc(pwCCL8); colorbar; caxis([-1 1]); title('pwCCs for Loading 8')
subplot(4,3,9);imagesc(pwCCL9); colorbar; caxis([-1 1]); title('pwCCs for Loading 9')
subplot(4,3,10);imagesc(pwCCL10); colorbar; caxis([-1 1]); title('pwCCs for Loading 10')

filepathb = obj.source_folder;
filenameb = 'blockperm_all_traces_2.646s_100x_pwCC';
filepathb = append(filepathb,filenameb);
saveas(figure(90),filepathb, 'png');
filepathc = append(filepathb, '.fig');
savefig(figure(90),filepathc,'compact');

%get values of the pwCCs and put them in a col_vector for visualizing distributions
col_vect_pwCCL1 = pwCCL1(~tril(ones(size(pwCCL1))));
col_vect_pwCCL2 = pwCCL2(~tril(ones(size(pwCCL2))));
col_vect_pwCCL3 = pwCCL3(~tril(ones(size(pwCCL3))));
col_vect_pwCCL4 = pwCCL4(~tril(ones(size(pwCCL4))));
col_vect_pwCCL5 = pwCCL5(~tril(ones(size(pwCCL5))));
col_vect_pwCCL6 = pwCCL6(~tril(ones(size(pwCCL6))));
col_vect_pwCCL7 = pwCCL7(~tril(ones(size(pwCCL7))));
col_vect_pwCCL8 = pwCCL8(~tril(ones(size(pwCCL8))));
col_vect_pwCCL9 = pwCCL9(~tril(ones(size(pwCCL9))));
col_vect_pwCCL10 = pwCCL10(~tril(ones(size(pwCCL10))));

figure(91);subplot(4,3,1);histogram(col_vect_pwCCL1,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 1')
subplot(4,3,2);histogram(col_vect_pwCCL2,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 2')
subplot(4,3,3);histogram(col_vect_pwCCL3,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 3')
subplot(4,3,4);histogram(col_vect_pwCCL4,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 4')
subplot(4,3,5);histogram(col_vect_pwCCL5,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 5')
subplot(4,3,6);histogram(col_vect_pwCCL6,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 6')
subplot(4,3,7);histogram(col_vect_pwCCL7,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 7')
subplot(4,3,8);histogram(col_vect_pwCCL8,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 8')
subplot(4,3,9);histogram(col_vect_pwCCL9,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 9')
subplot(4,3,10);histogram(col_vect_pwCCL10,20); xlim([-1, 1.05]); title('Distr of pwCCs for Loading 10')

filepathc = obj.source_folder;
filenamec = 'blockperm_all_traces_2.646s_100x_distr';
filepathc = append(filepathc,filenamec);
saveas(figure(91),filepathc, 'png');
filepathd = append(filepathc, '.fig');
savefig(figure(91),filepathd,'compact');