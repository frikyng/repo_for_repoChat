load('C:\Users\THE BEASTWO\Documents\MATLAB\curation_code\test_predicition.mat'); % in case you cleared everything

%% go through all experiments (in reverse)
for current_cell_idx = numel(fold):-1:1
    
    %current_cell_idx = 79 % change manually to just run one cell
    
    %% Get the output of the classifier for the current cell
    A = out{current_cell_idx};
    
    %% Plot classifier performances for every behaviour
    for beh_idx = 1:numel(A.peak_tp)
        plot_prediction(A.calcium, A.bin_beh{beh_idx}, A.peak_tp{beh_idx}, A.train_range{beh_idx}, A.prediction{beh_idx}, A.full_beh{beh_idx}, A.beh_type{beh_idx}, A.score{beh_idx}, beh_idx);
        
        %% Comment out the pause to go faster
        pause(0)
    end
    
    %% Rendering tweaking. you can change the pause etc...
    arrangefigures;
    pause(5);
    close all
end
