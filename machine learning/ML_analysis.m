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
            obj.settings.score_metrics = 'explained_variance'
        end
        
        function [meanvalue, sem_values, stats_results, values] = show_result(obj, condition, behaviour_filter)
            %% Display one or several behaviour and its/theirs shuffled control(s)
            if nargin < 2 || isempty(condition)
                condition = {'Prediction on all behaviours',' using events modulation on all ROIs'}; % condition is the file title
            end
            if nargin < 3 || isempty(behaviour_filter)
                behaviour_filter = '';
            end
            path = obj.get_file_path_from_title(condition);
            load(path);
            settings = obj.settings;
            settings.split_shuffle = false;
            [meanvalue, sem_values, fig_handle, stats_results, values]      = bar_chart(results,behaviour_filter,'','','','','',settings);
            save_myfig(fig_handle,'temp_chart','png')
        end
        
        function temporal_shuffle(obj, behaviour_filter)
            if nargin < 2 || isempty(behaviour_filter)
                behaviour_filter = '';
            elseif ~iscell(behaviour_filter)
                behaviour_filter = {behaviour_filter};
            end
            
            results     = {};            
            all_titles  = {
                'Prediction on all behaviours using events magnitudes on all ROIs',...
                'Prediction on all behaviours using events magnitude on all ROIs, but Behaviour and Ca2+ events are shuffled',...
                'Prediction on all behaviours using events modulation on all ROIs',...
                'Prediction on all behaviours using events modulation on all ROIs, but Behaviour and Ca2+ events are shuffled'};
            all_ticks   = {'Magnitude','Temporal Shuffle on Magnitude','Modulation','Temporal Shuffle on Modulation'};
            settings = obj.settings;
            settings.use_shuffle = false;
            x_ticks     = {};
            
            for cond = 1:numel(all_titles)
                new = load(obj.get_file_path_from_title(all_titles{cond}));
                if iscell(new.results{1})
                    new = {new.results{3}}; % {horzcat(new.results{:})}; % merge all groups
                else
                    new = {new.results};
                end
                newtick = all_ticks{cond};
            	results = [results, new];
                x_ticks  = [x_ticks, newtick];
            end
            if isempty(behaviour_filter)
                behaviour_filter = results{1}{1}.beh_type(~contains(results{1}{1}.beh_type, 'shuffle'));
            end            
            for beh = behaviour_filter
                [~,~, fig_handle]           = bar_chart(results, beh{1},'','',x_ticks,'','',settings);
                save_myfig(fig_handle,'temp_chart','png')
            end
        end
        
        function spatial_shuffle(obj, behaviour_filter)
            if nargin < 2 || isempty(behaviour_filter)
                behaviour_filter = '';
            elseif ~iscell(behaviour_filter)
                behaviour_filter = {behaviour_filter};
            end
            
            results     = {};            
            all_titles  = {
                'Prediction on all behaviours using events magnitudes on all ROIs',...
                'Prediction on all behaviours using events magnitude on all ROIs but ROIs are shuffled for each event',...
                'Prediction on all behaviours using events modulation on all ROIs',...
                'Prediction on all behaviours using events modulation on all ROIs but ROIs are shuffled for each event'};
            all_ticks   = {'Magnitude','Spatial Shuffle on Magnitude','Modulation','Spatial Shuffle on Modulation'};
            settings = obj.settings;
            settings.use_shuffle = false;
            x_ticks     = {};
            
            for cond = 1:numel(all_titles)
                new = load(obj.get_file_path_from_title(all_titles{cond}));
                if iscell(new.results{1})
                    new = {new.results{3}}; % {horzcat(new.results{:})}; % merge all groups
                else
                    new = {new.results};
                end
                newtick = all_ticks{cond};
            	results = [results, new];
                x_ticks  = [x_ticks, newtick];
            end
            if isempty(behaviour_filter)
                behaviour_filter = results{1}{1}.beh_type(~contains(results{1}{1}.beh_type, 'shuffle'));
            end            
            for beh = behaviour_filter
                [~,~, fig_handle]           = bar_chart(results, beh{1},'','',x_ticks,'','',settings);
                save_myfig(fig_handle,'temp_chart','png')
            end
        end
        
        function show_magnitude_and_modulation(obj, behaviour_filter)
            if nargin < 2 || isempty(behaviour_filter)
                behaviour_filter = '';
            elseif ~iscell(behaviour_filter)
                behaviour_filter = {behaviour_filter};
            end
            
            results = {};            
            all_titles  = { {'Prediction on all behaviours',' using events magnitudes on all ROIs'},...
                            {'Prediction on all behaviours',' using events modulation on all ROIs'}};
            all_ticks = {'Magnitude', 'Modulation'};
            x_ticks  = {};
            
            for cond = 1:numel(all_titles)
                new = load(obj.get_file_path_from_title(all_titles{cond}));
                if iscell(new.results{1})
                    new = {new.results{3}}; % {horzcat(new.results{:})}; % merge all groups
                else
                    new = {new.results};
                end
                newtick = all_ticks{cond};
            	results = [results, new];
                x_ticks  = [x_ticks, newtick];
            end
            if isempty(behaviour_filter)
                behaviour_filter = results{1}{1}.beh_type(~contains(results{1}{1}.beh_type, 'shuffle'));
            end
            for beh = behaviour_filter
                [~,~, fig_handle]           = bar_chart(results, beh{1},'','',x_ticks,'','',obj.settings);
                save_myfig(fig_handle,'temp_chart','png')
            end
        end
        
        function [meanvalue, sem_values, stats_results, values] = show_selected_conditions(obj, behaviour_filter, comparaison)
            if nargin < 2 || isempty(behaviour_filter)
                behaviour_filter = '';
            elseif ~iscell(behaviour_filter)
                behaviour_filter = {behaviour_filter};
            end
            if nargin < 3 || isempty(comparaison)
                comparaison = 0;                
            end
            
            all_results = {};         
            if ~comparaison
                all_titles  = {{'Prediction on all behaviours',' using events magnitude on perisomatic ROIs (averaged)'},... %{'Prediction on all behaviours',' using events modulation on perisomatic ROIs (averaged)'}
                               {'Prediction on all behaviours',' using events magnitude on all ROIs (averaged)'},...
                               {'Prediction on all behaviours',' using events modulation on the individual highest phate',[' of group ', num2str(1)]},...
                               {'Prediction on all behaviours';[' using events modulation on the 3 Highest Phates'];' (* ROIs)'},...
                               {'Prediction on all behaviours',' using events modulation on all ROIs'}};
                           all_ticks = {'Soma Average','Cell Average','Highest Phate','3 High Phates','All ROIs'};
            elseif comparaison == 1
                all_titles  = {{'Prediction on all behaviours',' using events magnitude on perisomatic ROIs (averaged)'},... % {'Prediction on all behaviours',' using events modulation on perisomatic ROIs (averaged)'}
                               {'Prediction on all behaviours',' using events magnitude on all ROIs (averaged)'},...
                               {'Prediction on all behaviours',' using events modulation on the individual highest phate',[' of group ?']},...
                               {'Prediction on all behaviours';[' using events modulation on the 3 Highest Phates'];' (* ROIs)'},...
                               {'Prediction on all behaviours',' using events modulation on all ROIs'}};
                           all_ticks = {'Soma Average','Cell Average','Highest Phate','3 High Phates','All ROIs'};
            elseif  comparaison == 2
                all_titles  = {{'Prediction on all behaviours',' using events magnitude on perisomatic ROIs (averaged)'},... % {'Prediction on all behaviours',' using events modulation on perisomatic ROIs (averaged)'}
                               {'Prediction on all behaviours',' using events modulation on the individual highest phate',[' of group ?']}};
                all_ticks = {'Soma Average','Highest Phate'};
            end
            
            x_ticks  = {};
             %                  {'Prediction on all behaviours';[' using events modulation on the highest beta coefs trained on ',behaviour_filter{:}]},...
            for cond = 1:numel(all_titles)
                if any(contains(all_titles{cond}, '?'))
                    %% When not specify a subgroup, go through all groups and find the highest score
                    temp_name = [[all_titles{cond}{:}]];
                    path = obj.get_file_path_from_title(temp_name);
                    [~, order] = sort(cellfun(@(x) str2double(regexp(x,  '\d+(?=\D*$)', 'match', 'once')), path));
                    path = path(order);
                    v = [];
                    for el = 1:numel(path)
                        [~, title] = fileparts(path{el});
                        load(obj.get_file_path_from_title(title));
                        [v(el), ~, ~, ~, ~]      = bar_chart(results,behaviour_filter,'','','',false,'',obj.settings);
                    end
                    [~, loc] = max(v);
                    [~, all_titles{cond}] = fileparts(path{loc});
                end
                new = load(obj.get_file_path_from_title(all_titles{cond}));
                if iscell(new.results{1})
                    new = {new.results{3}}; % {horzcat(new.results{:})}; % merge all groups
                else
                    new = {new.results};
                end
                newtick = all_ticks{cond};
            	all_results = [all_results, new];
                x_ticks  = [x_ticks, newtick];
            end
            if isempty(behaviour_filter)
                behaviour_filter = all_results{1}{1}.beh_type(~contains(all_results{1}{1}.beh_type, 'shuffle'));
            end
            values = {};
            for beh_idx = 1:numel(behaviour_filter)
                if ~nargout > 1
                    [~,~, fig_handle]           = bar_chart(all_results, behaviour_filter{beh_idx},'','',x_ticks,'','',obj.settings);
                    save_myfig(fig_handle,'temp_chart','png');
                else
                    [meanvalue{beh_idx}, sem_values{beh_idx}, ~, stats_results{beh_idx}, values{beh_idx}]           = bar_chart(all_results, behaviour_filter{beh_idx},'','',x_ticks,'','',obj.settings);
                end
            end
        end
        
       
        
%         function show_selected_conditions(obj, behaviour_filter)
%             if nargin < 2 || isempty(behaviour_filter)
%                 behaviour_filter = '';
%             elseif ~iscell(behaviour_filter)
%                 behaviour_filter = {behaviour_filter};
%             end
%             
%             subset = 1:6;
%             results = {};            
%             all_titles  = {{'Prediction on all behaviours',' using events modulation on perisomatic ROIs (averaged)'},...
%                            {'Prediction on all behaviours',' using events modulation on perisomatic ROIs'},...
%                            {'_ML_individual_clusters'},...
%                            {'Prediction on all behaviours';[' using events modulation on the 3 Highest Phates'];' (* ROIs)'},...
%                            {'Prediction on all behaviours';[' using events modulation on the 9 Highest Phates'];' (* ROIs)'},...
%                            {'Prediction on all behaviours',' using events modulation on all ROIs'}};
%             all_ticks = {'Cell Average','Perisomatic Segments','Random Sets','3 High Phates','9 High Phates','All ROIs'};
%             x_ticks  = {};
%             
%             for cond = subset
%                 new = load(obj.get_file_path_from_title(all_titles{cond}));
%                 if iscell(new.results{1})
%                     new = {new.results{3}}; % {horzcat(new.results{:})}; % merge all groups
%                 else
%                     new = {new.results};
%                 end
%                 newtick = all_ticks{cond};
%             	results = [results, new];
%                 x_ticks  = [x_ticks, newtick];
%             end
%             if isempty(behaviour_filter)
%                 behaviour_filter = results{1}{1}.beh_type(~contains(results{1}{1}.beh_type, 'shuffle'));
%             end
%             for beh = behaviour_filter
%                 [~,~, fig_handle]           = bar_chart(results, beh{1},'','',x_ticks,'','',obj.settings);
%                 save_myfig(fig_handle,'temp_chart','png')
%             end
%         end
        
        
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
            elseif  contains(fname, '?')
                path = strrep(path, '?', '*')
                path = sort(arrayfun(@(x) fullfile(folder, x.name), dir(path), 'UniformOutput', false));
            end
        end
    end
end

