classdef arboral_scan_dataset < handle
    %% This class handles meta analysis or arboreal scans
    % see meta_analyze_arboreal_scans for examples
    properties
        source_folder
        experiments
        
        export_folder
        expe_list      % List of experiments analyzed
        rec_list
        general_info   % 
        extracted_traces
        binned_data
        timescale
        event_fitting
        crosscorr      % Data relate to correlation study between all_traces
        dimensionality % Data related to dimensionality analysis of all_traces
        external_variables
        demo            = 0; %
        filter_win      = 0;
        filter_type     = 'gaussian';
        current_expe    = '';
        need_update     = [];
        rendering       = false;
    end
    
    methods
        function obj = arboral_scan_dataset(source_folder)
            obj.source_folder = source_folder;
            
            %% Find arboreal_scan_experiments
            fold = dir([source_folder,'/*-*-*_exp_*']);
            fold = fold([fold.isdir]);            
            for idx = 1:numel(fold)
                files = dir([fold(idx).folder,'/',fold(idx).name,'/*.mat']);
                if numel(files) == 0
                    warning([fold(idx).folder,'/',fold(idx).name, ' has not been extracted yet. A simple experiment (no processing) will be generated for you. Please wait'])
                    expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
                    expe.save(true);
                    files = dir([fold(idx).folder,'/',fold(idx).name,'/*.mat']);
                end
                
                for el = 1:numel(files)
                    m = matfile([files(el).folder,'/',files(el).name]);
                    if isprop(m,'obj') && isa(m.obj, 'arboreal_scan_experiment')
                        if isempty(obj.experiments)
                            obj.experiments = m.obj;
                        else
                            %% SQUEEZE IDENTIFY HERE
                            obj.experiments(idx) = m.obj;
                        end
                        break
                    end
                end
            end       
        end
        

%         
%         function expe_list = get.expe_list(obj)
%             names                   = vertcat(obj.rec_list{1}.name);
%             names                   = names(:,1:end-9);
%             names                   = cellstr(names);
%             [experiments_available, ~, idx2] = unique(names(cellfun('isclass', names, 'char')));           
%             expe_list               = experiments_available';
%         end
% 
%         function data_folders_per_exp = filter_expe_subset(obj, filter)
%             %% QQ UNTESTED AND UNFINISHED
%             if nargin < 2 || isempty(filter)
%                 filter = '';
%             end
%             
%             %% Run analysis on a subset if you specified a filter
%             if isempty(filter) 
%                 filter = [vertcat(obj.expe_list{:})];
%             end
% 
%             %% Regroup recordings per experiment
%             data_folders_per_exp = {};
%             count = 1;
%             for idx = 1:numel(obj.expe_list)    
%                 if ismember(obj.expe_list(idx), filter)
%                     data_folders_per_exp{count} = obj.rec_list{1}(obj.rec_list{2} == idx);
%                     count = count+1;
%                 end
%             end            
%         end
%         
%         function expe = add_experiment(obj, path)  
%             %% Identify if the experiment is already here to avoid duplicate
%             expe = [];
%             for el = 1:numel(obj.general_info)
%                 if isstruct(obj.general_info{el}) && isfield(obj.general_info{el}, 'short_path') && any(strcmp(obj.general_info{el}.short_path, path(1).name(:,1:end-9)))
%                     expe = el;
%                     break
%                 end                    
%             end
%             if isempty(expe)
%                 expe = numel(obj.general_info) + 1;
%             end
%             
%             %% Prepare the fields for the indicated experiment
%             obj.general_info{expe}      = [];
%             obj.extracted_traces{expe}  = [];
%             obj.binned_data{expe}       = [];
%             obj.timescale{expe}         = [];
%             obj.event_fitting{expe}     = [];
%             obj.crosscorr{expe}         = [];
%             obj.dimensionality{expe}    = [];
%             obj.external_variables{expe}= [];
%             %obj.expe_list{expe}         = [];
% 
%             obj.general_info{expe}.original_path = path;
%             obj.general_info{expe}.short_path = path(1).name(:,1:end-9);
%             
%             try
%                 obj.cleanup_expe();
%             catch
%                 obj.cleanup_expe();
%             end
%             expe = find(cellfun(@(x) contains(x.short_path, path(1).name(:,1:end-9)), obj.general_info));
%             obj.current_expe = expe;
%         end
%         
% 
%         
%         function cleanup_expe(obj, deep_cleanup, updated_path)
%             if nargin < 2 || isempty(deep_cleanup)
%                 deep_cleanup = false; source = '';
%             end
%             %% Remove empty location and sort by expe name
%             to_keep = find(~cellfun(@isempty, obj.general_info));
%             
%             obj.general_info = obj.general_info(to_keep);
%             obj.extracted_traces = obj.extracted_traces(to_keep);
%             obj.binned_data = obj.binned_data(to_keep);
%             obj.timescale = obj.timescale(to_keep);
%             obj.event_fitting = obj.event_fitting(to_keep);
%             obj.crosscorr = obj.crosscorr(to_keep);
%             obj.dimensionality = obj.dimensionality(to_keep);
%             obj.expe_list = obj.expe_list(to_keep);
%             obj.need_update = obj.need_update(to_keep);
%             obj.external_variables  = obj.external_variables(to_keep);
%             
%             [~,order] = sort(cellfun(@(x) x.short_path, obj.general_info, 'UniformOutput', false));
%             
%             obj.expe_list = obj.expe_list(order);
%             obj.general_info = obj.general_info(order);
%             obj.extracted_traces = obj.extracted_traces(order);
%             obj.binned_data = obj.binned_data(order);
%             obj.timescale = obj.timescale(order);
%             obj.event_fitting = obj.event_fitting(order);
%             obj.crosscorr = obj.crosscorr(order);
%             obj.dimensionality = obj.dimensionality(order);
%             obj.need_update = obj.need_update(order);
%             obj.external_variables  = obj.external_variables(order);
%             
%             if strcmp(deep_cleanup, 'extracted') && isfolder(obj.source_folder)
%                 source = obj.source_folder;
%                 if nargin > 2 && ~isempty(updated_path)
%                     source = updated_path;
%                 end
%             elseif strcmp(deep_cleanup, 'original') && isfolder(obj.source_folder) && ~isempty(updated_path)
%                 source = obj.source_folder;
%             elseif deep_cleanup
%                 source = '';
%                 warning('deep_cleanup must be "original", or "extracetd". If the extracted path changed or if you want to use the original data for cleaning, indicate a valid path in "updated_path" input')
%             end
%             
%             if ~isempty(source)
%                 for exp = 1:numel(obj.expe_list)
%                     if strcmp(deep_cleanup, 'original') 
%                     	to_remove = cellfun(@(x) ~isdir(x), obj.general_info{exp}.sources);
%                     elseif strcmp(deep_cleanup, 'extracted') 
%                         to_remove = cellfun(@(x) ~isdir(x), arrayfun(@(y) [y.folder,'\',y.name], obj.general_info{exp}.original_path, 'UniformOutput', false))';
%                     end
%                     if any(to_remove)
%                         obj.need_update(exp) = true;
%                         
%                         %% Remove unnecessary path
%                         all_names = {obj.general_info{exp}.original_path(to_remove).name};
%                         to_remove_too = find(cellfun(@(x) any(strcmp(x, all_names)), {obj.rec_list{1}.name}));
%                         obj.rec_list{1}(to_remove_too) = [];obj.rec_list{2}(to_remove_too) = [];
%                         obj.general_info{exp}.sources = obj.general_info{exp}.sources(~to_remove);
%                         obj.general_info{exp}.original_path = obj.general_info{exp}.original_path(~to_remove); 
%                     end                    
%                 end
%             end
%         end
%         
        
  
%         
%         function failed = update_external_metrics(obj, expe)
%             if nargin < 2 || isempty(expe)
%                 expe = 1:numel(obj.expe_list);
%             end
%             
%             failed = {};
% %             for el = 1:numel(obj.general_info)
% %                 current_sources = obj.general_info{el}.sources;
% %                 current_external_variables = {};
% %                 for rec = 1:numel(current_sources)
% %                     fprintf([current_sources{rec}, '\n']);
% %                     current_external_variables{rec} = [];                 
% %                 end  
% %                 obj.external_variables{el} = current_external_variables;
% %             end        
% %             
%             failed = {};
%             for el = expe
%                 current_sources = obj.general_info{el}.sources;
%                 current_external_variables = {};
%                 for rec = 1:numel(current_sources)
%                     fprintf([current_sources{rec}, '\n']);
%                     try
%                     parameters = analysis_params('source',current_sources{rec},'data_type','raw');
%                     [parameters, ~, ~, ~] = load_generic_scan(parameters);
%                     parameters.external_var = structfun(@(x) struct('time',x.time(:),'data',x.speed(:)), parameters.external_var,'UniformOutput', false);
%                     current_external_variables{rec} = parameters.external_var;
%                     catch err
%                         current_external_variables{rec} = [];
%                         failed{end + 1} = [{err}, current_sources{rec}];
%                     end                    
%                 end  
%                 obj.external_variables{el} = current_external_variables;
%             end            
%         end
% 
%         
%         function expe = identify(obj, filter)
%             expe = find(cellfun(@(x) contains(x, filter), obj.expe_list));
%             fprintf(['Required experiment is # ', num2str(expe),'\n'])
%             
%         end

        function mean_cc(obj)            
            CCs = {};
            for expe = 1:numel(obj.experiments)
                CCs{expe} = obj.experiments(expe).crosscorr;                
            end
            
            [~, biggest] = max(cellfun(@numel, CCs));
            
            CCs = cat(3, CCs{:});
            CCs = nanmean(CCs, 3);
            figure(1008);cla();imagesc(CCs); hold on;set(gcf,'Color','w');
            caxis([0,1]); hold on;
            xticks([1:size(CCs, 1)])
            yticks([1:size(CCs, 1)])
            colorbar; hold on;
            plot([1.5,1.5],[0.5,size(CCs, 1)+0.5],'k-');xticklabels(['whole cell',obj.experiments(biggest).binned_data.bin_legend]);xtickangle(45);
            plot([0.5,size(obj.crosscorr, 1)+0.5],[1.5,1.5],'k-');yticklabels(['whole cell',obj.experiments(biggest).binned_data.bin_legend]);
            title('peak amplitude correlation per subgroup');
            arrangefigures([1,2]); 
        end

        function plot_gallery(obj, mode, range, proj_axis)
            if nargin < 2 || isempty(mode)
                mode = 'corr';
            end
            if nargin < 3 || isempty(range)
                range = 1:numel(obj.experiments);%find(~cellfun(@isempty, obj.dimensionality));
            end
            if nargin < 4 || isempty(proj_axis)
                proj_axis = 'xz';
            end
            
            %% Now generate the figures
            [all_trees,  all_soma, all_values]  = deal({});
            for expe = range
                if strcmp(mode, 'corr')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).plot_corr_tree();
                elseif strcmp(mode, 'dim')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).obj.plot_dim_tree(0); 
                elseif strcmp(mode, 'dist')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).obj.plot_dist_tree(0); 
                elseif strcmp(mode, 'order')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).obj.plot_ord_tree(0); 
                end
                close all;
            end

            
            L23 = 100 * double(~strcmp(proj_axis, 'xy'));
            L5a = 444 * double(~strcmp(proj_axis, 'xy'));
            L5b = 644 * double(~strcmp(proj_axis, 'xy'));       
            figure(); axis equal;set(gcf,'Color','w');hold on;
            x_offset = 0; y_offset = 0; max_fig_w = 5000; y_max_in_row = 0;
            for expe = 1:numel(all_values)
                if ~isempty(all_values{expe})
                    if x_offset > (max_fig_w-500) || expe == 1
                        x_offset = 0;
                        if expe > 1
                            y_offset     = y_offset + y_max_in_row + 100;
                            y_max_in_row = 0;
                        end
                        plot([0,max_fig_w],[y_offset,y_offset],'k-'); hold on;
                        plot([0,max_fig_w],[y_offset+L23,y_offset+L23],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
                        plot([0,max_fig_w],[y_offset+L5a,y_offset+L5a],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
                        plot([0,max_fig_w],[y_offset+L5b,y_offset+L5b],'Color',[0.3,0.3,0.3],'LineStyle','--'); hold on;
                    end
                    if strcmpi(proj_axis, 'xz')
                        req_x = 1; req_y = 3; req_z = 2;
                    elseif strcmpi(proj_axis, 'yz')
                        req_x = 2; req_y = 3; req_z = 1;
                    elseif strcmpi(proj_axis, 'xy')
                        req_x = 1; req_y = 2; req_z = 3;
                    else
                        error('proj axis must be xz, yz or xy')                        
                    end
                    shifted_x_coor = all_values{expe}{req_x};
                    shifted_z_coor = all_values{expe}{req_y};
                    
                    zero_shifted_coor = -min(shifted_x_coor) + x_offset;
                    shifted_x_coor = shifted_x_coor + zero_shifted_coor;
                    X = shifted_x_coor;
                    Y = shifted_z_coor+y_offset;
                    Z = all_values{expe}{req_z};
                    S = all_values{expe}{4};
                    HP = surface('XData'    ,[X,X],...
                                 'YData'    ,[Y,Y],...
                                 'ZData'    ,[Z,Z],...
                                 'CData'    ,[S,S],...
                                 'facecol'  ,'no',...
                                 'edgecol'  ,'interp',...
                                 'linew'    ,1); hold on;
                    scatter(all_soma{expe}(req_x) + zero_shifted_coor,all_soma{expe}(req_y)+y_offset,50,'k','filled'); hold on;
                    x_offset = max(shifted_x_coor);
                    y_max_in_row = max([y_max_in_row, max(Y-y_offset), all_soma{expe}(req_y), L5b]);
                end
            end
            ylim([0, y_offset + y_max_in_row]);
            xlim([0,max_fig_w]); hold on;set(gca, 'Ydir', 'reverse');colorbar                                                       
        end   

%         function extracted_traces = get.extracted_traces(obj)
%             extracted_traces = obj.extracted_traces;
%             if ~isempty(obj.current_expe) && ~isempty(extracted_traces{obj.current_expe}) && any(obj.filter_win)
%                 temp = extracted_traces{obj.current_expe};
%                 extracted_traces{obj.current_expe} = cellfun(@(x) smoothdata(x, 1, obj.filter_type, obj.filter_win), temp, 'UniformOutput', false);
%             end
%         end
    end
end

