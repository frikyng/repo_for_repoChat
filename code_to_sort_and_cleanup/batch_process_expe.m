%% Set some analysis options
opengl software
save_data   = true;
demo        = 0;                % 0 / false for no demo ; 1 for fitting demo ; 2 for fitting and scaling demo
filter      = '';
condition   = {'distance', 100};
    
%% Define the folders to use
%% source_folder points indicates where your ribbon scan experiments are 
% setting_file_path = 'D:\Curated Data\settings.txt';
% source_folder = 'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER\Ribbon Scan paper\code for paper\test7\meta\';
export_folder = 'C:\Users\vanto\Desktop\new_meta_nov_2020'

if exist('dataset','var') && isa(dataset, 'arboral_scan_meta_analysis')
    
elseif isfile([export_folder, '/summary_full.mat'])
    %% Resume analysis
    load([export_folder, '/summary_full.mat']);
else
    %% Start new analysis
    dataset                 = arboreal_scan_dataset(source_folder, export_folder);
    setting_file_path       = SETTINGS_FILE                     %'D:\Curated Data\settings.txt';
    source_folder           = [EXPORT_PATH, '\meta']            %'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER\Ribbon Scan paper\code for paper\26-10-2020\meta\';
    export_folder           = [pwd, '\new_meta_nov_2020'];
    export_folder           = parse_paths(export_folder);
end
%data_folders_per_exp    = dataset.filter_expe_subset();

%data_folders_per_exp    = data_folders_per_exp(dataset.need_update)

failed_analysis         = {};
failed_factoran         = {};
whitebg('w')

%% Now, do the analysis expe-by-expe
for expe = 1:numel(data_folders_per_exp) % review expe 23 for fitting warning, saturated peak scaling and one group having an overestimated baseline
    dataset.process(expe, data_folders_per_exp, condition, '', true)
end