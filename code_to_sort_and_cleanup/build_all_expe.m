top_export_folder = 'C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans_2\extracted_arboreal_scans';
fold = dir([top_export_folder,'/*-*-*_exp_*']);
fold = fold([fold.isdir]);
errors = {};

%to_redo = numel(fold)-7:-1:1%find(~cellfun(@isempty, errors));
for idx = redo %;%(numel(fold)):-1:1% PB  5     7    11    13    16
    try
%         expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
%%         if contains(errors{idx}.message, 'value')
%             [~,~,expe_f] = parse_paths(expe.ref.data_folder, true);
%             meta_batch_process_ribbon_scan(expe_f,'', '','C:\Users\THE BEASTWO\Documents\MATLAB\arboreal_scans');
%         end
        expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
        expe.auto_save_analysis = 1;
        expe.auto_save_figures = 0;
        expe.process({'depth',50},[ceil(1/nanmedian(expe.timescale.sr)), 0]); %% Add grouping settings here
       % errors{idx} = [];
    catch err
        errors{idx} = err;
        1
    end
    close all
end

% Columns 1 through 26
% 
%      1     2     3     4     5     6     9    12    22    25    26    27    28    33    37    39    41    47    52    55    56    59    60    64    68    70
% 
%   Columns 27 through 35
% 
%     73    74    79    81    84    85    86    90    94

%%% 'C:/Users/THE
%%% BEASTWO/Documents/MATLAB/arboreal_scans/extracted_arboreal_scans/2019-09-12_exp_7/'
%%% has an issue (it's the duplicated expe)