extract_folder      = 'D:\SuperBackup Laptop TO DELETE\LaCie\Curated Data\';
raw_export_folder   = 'C:\Users\Antoine.Valera\MATLAB\newextraction';
L23_go_extracted    = {'2019-11-05_exp_3','2019-10-03_exp_5','2019-10-01_exp_1','2019-09-12_exp_2','2019-09-07_exp_7','2019-08-02_exp_1'}; % use shortnames
L23_go_original     = cellfun(@(x) parse_paths([extract_folder, '\', strrep(x, '_exp','\experiment')]), L23_go_extracted, 'UniformOutput', false);


processed_export_folder     = [raw_export_folder, '_raw_zscored'];

use_existing                = true;
errors = {}
for cell = L23_go_extracted
    expe_folder = parse_paths([processed_export_folder,'/',  cell{1},'/']);
    if isfile([expe_folder, cell{1}, '.mat'])
        expe_folder = [expe_folder, cell{1}, '.mat'];
    end
%     try
    Demo_FINAL_SERIES_OF_ANALYSIS();
%     end
    close all
end