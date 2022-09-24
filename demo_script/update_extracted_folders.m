%% If you removed folders ater extraction in your original TOP FOLDER, clean the extracted folders from here

TOP_FOLDER = 'D:\Curated Data\' % Top folder used for meta_batch_process_ribbon_scan
extracted_folder= 'C:\Users\vanto\Documents\MATLAB\RIBBON_SCAN_PAPER\Ribbon Scan paper\code for paper\26-10-2020\meta\' % export folder of meta_batch_process_ribbon_scan


[~, day_folders] = list_day_folders(TOP_FOLDER, false);

tic;

all_tags = {};
for day_idx = numel(day_folders):-1:1
    [~, expe_folders] = list_expe_folders(day_folders{day_idx}, false);
    for expe_idx = 1:numel(expe_folders)
        [~, data_folders] = list_data_folders(expe_folders{expe_idx}, false);
        for data_idx = 1:numel(data_folders)
            data_folder = data_folders{data_idx};
            [~,~,~,~,~,~,all_tags{end+1},~] = parse_paths(data_folder);            
        end
    end
end


extracted = dir([extracted_folder,'**\*-*-*_exp_*_*-*-*']);
extracted = extracted([extracted.isdir]);
extracted = arrayfun(@(x) [x.folder,'\',x.name], extracted, 'UniformOutput', false);

%% Folders to delete
to_delete = extracted(~(cellfun(@(x) any(contains(x, all_tags)), extracted')))

% for el = 1:numel(to_delete)
%     rmdir(to_delete{el},'s');
% end