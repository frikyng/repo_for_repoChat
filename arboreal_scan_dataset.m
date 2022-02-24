classdef arboreal_scan_dataset < handle
    %% This class handles meta analysis or arboreal scans
    % see meta_analyze_arboreal_scans for examples
    properties
        source_folder
        experiments;
        valid_expe = [];
    end
    
    methods
        function obj = arboreal_scan_dataset(source_folder)
            if isempty(obj.experiments)
                obj.experiments = arboreal_scan_experiment();% initial empty experiment
            end
            obj.source_folder = source_folder;
            
            %% Find arboreal_scan_experiments
            fold = dir([source_folder,'/*-*-*_exp_*']);
            fold = fold([fold.isdir]);            
            for idx = 1:(numel(fold))%1:numel(fold)
                files = dir([fold(idx).folder,'/',fold(idx).name,'/*.mat']);
%                 if numel(files) == 0
%                     warning([fold(idx).folder,'/',fold(idx).name, ' has not been extracted yet. A simple experiment (no processing) will be generated for you. Please wait'])
%                     expe = arboreal_scan_experiment([fold(idx).folder,'/',fold(idx).name]);
%                     expe.save(true);
%                     files = dir([fold(idx).folder,'/',fold(idx).name,'/*.mat']);
%                 end
                
               % for el = 1:numel(files)
               fname = [files(1).folder,'/',files(1).name];
               if strcmp(obj.fast_class_check(fname), 'arboreal_scan_experiment')
                   fprintf(['loading ',files(1).name,' ... please wait\n'])
                   m = load(fname);  
%                    test = m.obj;
%                    close all
%                    test.find_events
%                    test.rescale_traces      
%                    %% CHECK 3,9,15    for low gain,
%                    %% CHeck 10 for low offset          
%                    test.plot_detected_events 
                   obj.experiments(idx) = m.obj;
               end
               
                    
                  % if isfield(m,'obj') && isa(m.obj, 'arboreal_scan_experiment')
%                         if isempty(obj.experiments)
%                             obj.experiments = m.obj;
%                         else
%                             %% SQUEEZE IDENTIFY HERE
%                             try
                            
%                             catch
%                                 1
%                             end
                  %  end
                        %break
                  %  end
               % end
            end       
        end
        
        function class_name = fast_class_check(obj,fname)
            class_name = '';
            fid = H5F.open(fname); % Use low level H5F builtin to open
            try 
                 obj_id = H5O.open(fid,'/obj','H5P_DEFAULT');
                 attr_id = H5A.open(obj_id,'MATLAB_class');
                 class_name = [H5A.read(attr_id)'];

                 H5O.close(obj_id);
                 H5A.close(attr_id);
            end
            H5F.close(fid);
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
%             for el = 1:numel(obj.rescaling_info)
%                 if isstruct(obj.rescaling_info{el}) && isfield(obj.rescaling_info{el}, 'short_path') && any(strcmp(obj.rescaling_info{el}.short_path, path(1).name(:,1:end-9)))
%                     expe = el;
%                     break
%                 end                    
%             end
%             if isempty(expe)
%                 expe = numel(obj.rescaling_info) + 1;
%             end
%             
%             %% Prepare the fields for the indicated experiment
%             obj.rescaling_info{expe}      = [];
%             obj.extracted_traces{expe}  = [];
%             obj.binned_data{expe}       = [];
%             obj.timescale{expe}         = [];
%             obj.event.fitting{expe}     = [];
%             obj.crosscorr{expe}         = [];
%             obj.dimensionality{expe}    = [];
%             obj.external_variables{expe}= [];
%             %obj.expe_list{expe}         = [];
% 
%             obj.rescaling_info{expe}.original_path = path;
%             obj.rescaling_info{expe}.short_path = path(1).name(:,1:end-9);
%             
%             try
%                 obj.cleanup_expe();
%             catch
%                 obj.cleanup_expe();
%             end
%             expe = find(cellfun(@(x) contains(x.short_path, path(1).name(:,1:end-9)), obj.rescaling_info));
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
%             to_keep = find(~cellfun(@isempty, obj.rescaling_info));
%             
%             obj.rescaling_info = obj.rescaling_info(to_keep);
%             obj.extracted_traces = obj.extracted_traces(to_keep);
%             obj.binned_data = obj.binned_data(to_keep);
%             obj.timescale = obj.timescale(to_keep);
%             obj.event.fitting = obj.event.fitting(to_keep);
%             obj.crosscorr = obj.crosscorr(to_keep);
%             obj.dimensionality = obj.dimensionality(to_keep);
%             obj.expe_list = obj.expe_list(to_keep);
%             obj.need_update = obj.need_update(to_keep);
%             obj.external_variables  = obj.external_variables(to_keep);
%             
%             [~,order] = sort(cellfun(@(x) x.short_path, obj.rescaling_info, 'UniformOutput', false));
%             
%             obj.expe_list = obj.expe_list(order);
%             obj.rescaling_info = obj.rescaling_info(order);
%             obj.extracted_traces = obj.extracted_traces(order);
%             obj.binned_data = obj.binned_data(order);
%             obj.timescale = obj.timescale(order);
%             obj.event.fitting = obj.event.fitting(order);
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
%                     	to_remove = cellfun(@(x) ~isdir(x), obj.rescaling_info{exp}.sources);
%                     elseif strcmp(deep_cleanup, 'extracted') 
%                         to_remove = cellfun(@(x) ~isdir(x), arrayfun(@(y) [y.folder,'\',y.name], obj.rescaling_info{exp}.original_path, 'UniformOutput', false))';
%                     end
%                     if any(to_remove)
%                         obj.need_update(exp) = true;
%                         
%                         %% Remove unnecessary path
%                         all_names = {obj.rescaling_info{exp}.original_path(to_remove).name};
%                         to_remove_too = find(cellfun(@(x) any(strcmp(x, all_names)), {obj.rec_list{1}.name}));
%                         obj.rec_list{1}(to_remove_too) = [];obj.rec_list{2}(to_remove_too) = [];
%                         obj.rescaling_info{exp}.sources = obj.rescaling_info{exp}.sources(~to_remove);
%                         obj.rescaling_info{exp}.original_path = obj.rescaling_info{exp}.original_path(~to_remove); 
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
% %             for el = 1:numel(obj.rescaling_info)
% %                 current_sources = obj.rescaling_info{el}.sources;
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
%                 current_sources = obj.rescaling_info{el}.sources;
%                 current_external_variables = {};
%                 for rec = 1:numel(current_sources)
%                     fprintf([current_sources{rec}, '\n']);
%                     try
%                     parameters = analysis_params('source',current_sources{rec},'data_type','raw');
%                     [parameters, ~, ~, ~] = load_generic_scan(parameters);
%                     parameters.external_var = structfun(@(x) struct('time',x.time(:),'data',x.value(:)), parameters.external_var,'UniformOutput', false);
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







        function expe = identify(obj, filter)
            expe = find(arrayfun(@(x) contains(x.ref.data_folder, filter), obj.experiments));
            fprintf(['Required experiment is # ', num2str(expe),'\n'])            
        end

        function CCs = mean_cc(obj, depth_range, corr_type)  % if that's the overall mean, we can keep it   
            if nargin < 2 || isempty(depth_range)
                min_depth   = -inf;
                max_depth   = inf;
            else
                min_depth   = depth_range(1);
                max_depth   = depth_range(2);
            end
            if nargin < 3 || isempty(corr_type)
                corr_type   = 'peaks_groups';
            end
            
            %% Get CCs
            CCs_ori = {};
            for expe = 1:numel(obj.experiments)
                obj.experiments(expe).cc_mode = corr_type;
                CCs_ori{expe} = obj.experiments(expe).crosscorr;                
            end            
            
            %% Filter CCs
            subset          = obj.valid_expe(find(arrayfun(@(x) x.ref.soma_location(3) > min_depth, obj.experiments(obj.valid_expe)) & arrayfun(@(x) x.ref.soma_location(3) < max_depth, obj.experiments(obj.valid_expe))));
            CCs_ori         = CCs_ori(subset);
            
            %% Make sure there is no empty matrix
            CCs_ori = CCs_ori(~cellfun(@isempty, CCs_ori));
            
            %% Full cell only
           % CCs_ori = CCs_ori(~cellfun(@(x) any(isnan(x(:))), CCs_ori));
            
%             %% Big cell only
            % CCs_ori = CCs_ori(cellfun(@(x) any(isnan(x(:))), CCs_ori));
            
            %% Create output matrix and insert values
            if ~isempty(CCs_ori)
                [~, biggest]    = max(cellfun(@numel, CCs_ori)); 
                biggest_idx     = subset(biggest); % use real index
                CCs             = repmat({NaN(size(CCs_ori{biggest}))},1,numel(CCs_ori));
                for expe = 1:numel(CCs)
                    new = CCs_ori{expe};
                    CCs{expe}(1:size(new, 1),1:size(new, 2)) = new;
                end
                CCs = nanmean(cat(3, CCs{:}), 3);

                %% Plot result
                figure(10080);clf();imagesc(CCs); hold on;
                set(gcf,'Color','w');axis image;title('peak amplitude correlation per subgroup');caxis([0,1]); hold on;
                colorbar; hold on;
                plot([1.5,1.5],[0.5,size(CCs, 1)+0.5],'k-');plot([0.5,size(CCs, 1)+0.5],[1.5,1.5],'k-');
                xticks(1:size(CCs, 1));xticklabels(['whole cell',obj.experiments(biggest_idx).binned_data.bin_legend]);xtickangle(45); 
                yticks(1:size(CCs, 1));yticklabels(['whole cell',obj.experiments(biggest_idx).binned_data.bin_legend]);  

                test = CCs; test(test == 1) = NaN;
                figure(10081);cla();bar(nanmean(test))
                nanmean(test(:,1))
                arrangefigures([1,2]); 
            else
                CCs = [];
            end
        end
        
        function valid_expe = get.valid_expe(obj)
            valid_expe = find(~cellfun(@isempty, {obj.experiments.source_folder}));
        end

        function plot_gallery(obj, mode, range, proj_axis)
            if nargin < 2 || isempty(mode)
                mode = 'corr';
            end
            if nargin < 3 || isempty(range)
                range = obj.valid_expe;
            end
            if nargin < 4 || isempty(proj_axis)
                proj_axis = 'xz';
            end
            
            %% Now generate the figures
            [all_trees,  all_soma, all_values]  = deal({});
            %range(2) = []
            for expe = range
                if strcmp(mode, 'corr')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).plot_corr_tree();
                elseif strcmp(mode, 'dim')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).plot_dim_tree(0); 
                elseif strcmp(mode, 'dist')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).plot_dist_tree(0); 
                elseif strcmp(mode, 'order')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.experiments(expe).plot_ord_tree(0); 
                end
                close all;
            end

            
            L23 = 100 * double(~strcmp(proj_axis, 'xy'));
            L5a = 444 * double(~strcmp(proj_axis, 'xy'));
            L5b = 644 * double(~strcmp(proj_axis, 'xy'));       
            figure(); axis equal;set(gcf,'Color','w');hold on; % FF
            x_offset = 0; y_offset = 0; max_fig_w = 5000; y_max_in_row = 0;
            for expe = range
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
                    HP = surface('XData'    ,[X;X],...
                                 'YData'    ,[Y;Y],...
                                 'ZData'    ,[Z;Z],...
                                 'CData'    ,[S;S],...
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
        
        function plot_group_correlation_data(obj, filter)
                    % Plot Group pearson  

        %     if ~isempty(filter)
        %         filter = cellfun(@(x) size(x, 2), [results.cumsum]) < 9;
        %         results = structfun(@(x) x(filter), results, 'UniformOutput', false);
        %     end       
        %             
        
            [biggest_gp,idx]      = max(cellfun(@(x) size(x,2) ,{obj.experiments.crosscorr}));
            CCs             = {obj.experiments.crosscorr};
            valid           = ~cellfun(@isempty, {obj.experiments.crosscorr}) & ~cellfun(@(x) any(isnan(x(:))), CCs);
            
            CCs             = CCs(valid);
            %obj.cumsum      = obj.cumsum(~cellfun(@isempty, obj.cumsum));
            peaks           = {obj.experiments.event};
            peaks           = peaks(valid);
            peaks           = cellfun(@(x) x.fitting.post_correction_peaks, peaks, 'UniformOutput', false);
            
             all_cc = cellfun(@(x) [x(1:size(x,1),1:size(x,2)), NaN(size(x,1),biggest_gp-size(x,2))], CCs, 'UniformOutput', false);
             all_cc = cellfun(@(x) [x(1:size(x,1),1:size(x,2)); NaN(biggest_gp-size(x,1),size(x,2))], all_cc, 'UniformOutput', false);
             all_cc = cat(3,all_cc{:});  
             
%             figure();imagesc(nanmean(all_cc,3));axis equal;colorbar
%             figure();imagesc(sum(~isnan(all_cc),3));axis equal;colorbar

            first_gp_cc = cellfun(@(x) x(2,:), CCs, 'UniformOutput', false);
            first_gp_cc = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,biggest_gp-size(x,2))], first_gp_cc, 'UniformOutput', false);
            first_gp_cc = vertcat(first_gp_cc{:});
            figure();bar(nanmean(first_gp_cc));hold on; errorbar(nanmean(first_gp_cc), nanstd(first_gp_cc)./sqrt(sum(~isnan(first_gp_cc))), 'k') ;
            xticklabels(['whole cell',obj.experiments(idx).binned_data.bin_legend]);xtickangle(45);

            % Std per event (need behaviour) 
            %test = cellfun(@(x) nanstd(x'), results.peaks, 'UniformOutput', false); % std per event

            % Variability per group (from mean)
            test = cellfun(@(x) nanvar(x)./nanmean(x), peaks, 'UniformOutput', false); % variance from mean, per subgroup
            test = cellfun(@(x) [x(1,1:size(x,2)), NaN(1,biggest_gp-size(x,2))], test, 'UniformOutput', false);
            test = vertcat(test{:});
            figure();bar(nanmean(test));hold on; errorbar(nanmean(test), nanstd(test)./sqrt(sum(~isnan(test))), 'k') 
            xticklabels(['whole cell',obj.experiments(idx).binned_data.bin_legend]);xtickangle(45);
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

