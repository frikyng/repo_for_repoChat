function partition = block_shuffle(timepoints, shuffle_window_pts, total_tp, holdout)

    shuffle_window_pts = abs(shuffle_window_pts);

	%% Cut data into windows
    win_start = 1:shuffle_window_pts:total_tp;
    windows  = {};
    for t = win_start
         current_win = t:t+shuffle_window_pts-1;
         windows{end+1} = timepoints(ismember(timepoints, current_win));      
    end
    windows = windows(~cellfun(@(x) isempty(x), windows));
    
    %% Shuffle windows
    windows = windows(randperm(numel(windows)));
    
    %% Extract held out blocks
    n_pts_target = ceil(numel(timepoints) * holdout);
    n_cells_to_use = find(cumsum(cellfun(@(x) numel(x), windows)) < n_pts_target, 1, 'last')+1;
    heldout_batch = windows(1:n_cells_to_use);
    training_batch = windows(n_cells_to_use+1:end);

    partition = {};
    partition.x_train = [training_batch{:}];
    partition.x_test  = [heldout_batch{:}];
    
    [~, partition.x_train] = ismember(partition.x_train, timepoints);
    [~, partition.x_test] = ismember(partition.x_test, timepoints);
%     %% Rearrange order
%     timepoints = [windows{:}]; 
%     data = data(:, timepoints);
%     processed_behaviours = processed_behaviours(:,timepoints);
%     

end