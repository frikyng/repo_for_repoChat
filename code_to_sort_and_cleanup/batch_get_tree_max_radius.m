function out = batch_get_tree_max_radius(top_folder, settings_file_path)
    if nargin < 2 || isempty(settings_file_path) 
        [~,~,~,~,top_folder] = parse_paths(top_folder, true);
        settings_file_path = [top_folder, 'settings.txt'];
    end
    
    %% Load batch processing settings
    batch_params = get_paths_from_settings(settings_file_path); 

    %% Now generate the figures
    radius = batch_process_helper(  top_folder          ,...
                                    {'expe_folder'}       ,...
                                    @(x,batch_params) get_tree_max_radius(x, batch_params), batch_params);
                                
                                
    [~, exp] = list_expe_folders(top_folder);  
    out = cellfun(@(x, y) [x, y], exp',radius,'UniformOutput',false);
end


% 
% TOP_FOLDER = 'D:\Curated Data\';
% 
% %% Get radius for all experiments listed in settings.txt
% radius = batch_get_tree_max_radius(TOP_FOLDER, [TOP_FOLDER,'\settings.txt']);
% %radius = cellfun(@(x) x{1}, radius);
% 
% % figure();hist(radius,0:20:400)
% % xlim([0,400])
% % xlabel('Diameter (um)')
% % ylabel('# experiments')
% 
% %% 
% fill_specific_spreadsheet_column(TOP_FOLDER, 'Curated log 20-10-2020.xlsx', radius, false, 'apicaltuftradialdiameterms')
% 
% 
% 
% TOP_FOLDER = 'D:\Curated Data\';
% 
% %% Get radius for all experiments listed in settings.txt
% radius = batch_get_tree_max_radius(TOP_FOLDER, [TOP_FOLDER,'\settings.txt']);
% 
% %% Update the "lab book" spreasheet with the missing data
% fill_specific_spreadsheet_column(TOP_FOLDER, 'Curated log 20-10-2020.xlsx', radius, false, 'apicaltuftradialdiameterms')
