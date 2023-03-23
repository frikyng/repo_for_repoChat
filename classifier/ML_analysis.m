classdef ML_analysis < handle
    %PAPER_FIGURES Assuming you extracted data using the demo script, use
    %   this class to build the figures for the paper
    
    properties
        processed_export_folder = '';
        cell_tag                = '';
        expe_folder             = '';
        settings                = [];
    end
    
    methods
        function obj = ML_analysis(processed_export_folder, cell_tag)
            obj.processed_export_folder = processed_export_folder;
            obj.cell_tag                = cell_tag;

            obj.expe_folder             = parse_paths([obj.processed_export_folder,'/',  obj.cell_tag,'/']);
            
            use_existing                = true;
            assignin('base','use_existing',use_existing);
            
            obj.set_max_score();
        end
        
        function set.cell_tag(obj, cell_tag)
            if iscell(cell_tag)
                cell_tag = cell_tag{1};
            end
            obj.cell_tag = cell_tag;
        end
        
        function set.expe_folder(obj, expe_folder)
            if isfile([expe_folder,'/', obj.cell_tag, '.mat'])
                expe_folder = [expe_folder,'/', obj.cell_tag, '.mat'];
            end
            obj.expe_folder = parse_paths(expe_folder);
        end
        
        function obj = set_max_score(obj)
            path            = obj.get_file_path_from_title({'_lag_All_ROIs'});
            load(path);
            [max_score, max_lag] = max(smoothdata(cell2mat(R),2,'movmean',2)');
            x               = lag_list/sr;
            max_lag         = x(max_lag);
            obj.settings.max_score   = max_score;
            obj.settings.max_lag     = max_lag;
        end
        
        function show_shuffle(obj, title, behaviour_filter)
            if nargin < 2 || isempty(title)
                title = {'Prediction on all behaviours',' using events amplitude modulations on all ROIs'};
            end
            if nargin < 3 || isempty(behaviour_filter)
                behaviour_filter = '';
            end
            path = obj.get_file_path_from_title(title);
            load(path);
            [~,~, fig_handle]           = bar_chart(results,behaviour_filter,'','','','','',obj.settings);
            save_myfig(fig_handle,'temp_chart','png')
        end
        
        function show_selected_conditions(obj, behaviour_filter)
            if nargin < 2 || isempty(behaviour_filter)
                behaviour_filter = '';
            elseif ~iscell(behaviour_filter)
                behaviour_filter = {behaviour_filter};
            end
            cell_average = load(obj.get_file_path_from_title({'Prediction on all behaviours',' using events amplitude modulations on perisomatic ROIs (averaged)'}));
            cell_average = cell_average.results;
            random_sets = load(obj.get_file_path_from_title({'_ML_individual_clusters'}));
            random_sets = random_sets.results{2};
            all_segments = load(obj.get_file_path_from_title({'Prediction on all behaviours',' using events amplitude modulations on all ROIs'}));
            all_segments = all_segments.results;
            N_phates = 3;
            high_phates = load(obj.get_file_path_from_title({'Prediction on all behaviours';[' using events amplitude modulations on the ',num2str(N_phates),' Highest Phates'];' (* ROIs)'}));
            high_phates = high_phates.results;
            
            results = {cell_average, random_sets, high_phates, all_segments};
            
            if isempty(behaviour_filter)
                behaviour_filter = results{1}{1}.beh_type(~contains(results{1}{1}.beh_type, 'shuffle'));
            end
            for beh = behaviour_filter
                [~,~, fig_handle]           = bar_chart(results, beh{1},'','',{'Cell Average','Random Sets','High Phates','All ROIs'},'','',obj.settings);
                save_myfig(fig_handle,'temp_chart','png')
            end
        end
        
        
        function path = get_file_path_from_title(obj,title)
            if iscolumn(title)
                title = title';
            end
            fname   = [title,{'.mat'}];
            fname   = [fname{:}];
            folder  = parse_paths(get_nth_parent_folder(obj.expe_folder, 0));  
            path    = [folder, fname];
            if contains(fname, '*')
                fname   = dir(path);
                path    = [folder, fname.name];
            end
        end
    end
end

