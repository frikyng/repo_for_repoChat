classdef arboral_scan_meta_analysis < handle
    %% This class handles meta analysis or arboreal scans
    % see meta_analyze_arboreal_scans for examples
    properties
        source_folder
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
        function obj = arboral_scan_meta_analysis(source_folder, export_folder)
            if nargin < 2 || isempty(export_folder)
                export_folder = [export_folder, '/meta/'];
            end
            
            %% Fix paths
            obj.source_folder       = parse_paths(source_folder);
            obj.export_folder       = parse_paths(export_folder);

            %% Create export folder if necessary
            if ~isdir(obj.export_folder)
                mkdir(obj.export_folder);
            end

            %% Get all analysed data folder in the source folder, for the type specified
            all_recordings          = dir([obj.source_folder,'/**/*-*-*_exp_*_*-*-*']);
            all_recordings          = all_recordings([all_recordings(:).isdir]);
            all_recordings          = all_recordings(~(arrayfun(@(x) strcmp(x.name, '.'), all_recordings) | arrayfun(@(x) strcmp(x.name, '..'), all_recordings)));
            names                   = vertcat(all_recordings.name);
            names                   = names(:,1:end-9);
            names                   = cellstr(names);
            [~, ~, idx2] = unique(names(cellfun('isclass', names, 'char')));  
            obj.rec_list            = {all_recordings, idx2};
            
            obj.general_info        = cell(size(obj.expe_list));
            obj.extracted_traces    = cell(size(obj.expe_list));
            obj.binned_data         = cell(size(obj.expe_list));
            obj.timescale           = cell(size(obj.expe_list));
            obj.event_fitting       = cell(size(obj.expe_list));
            obj.crosscorr           = cell(size(obj.expe_list));
            obj.dimensionality      = cell(size(obj.expe_list));
            obj.external_variables  = cell(size(obj.expe_list));
            obj.need_update         = true(size(obj.expe_list));            
        end
        
        function expe_list = get.expe_list(obj)
            names                   = vertcat(obj.rec_list{1}.name);
            names                   = names(:,1:end-9);
            names                   = cellstr(names);
            [experiments_available, ~, idx2] = unique(names(cellfun('isclass', names, 'char')));           
            expe_list               = experiments_available';
        end

        function data_folders_per_exp = filter_expe_subset(obj, filter)
            %% QQ UNTESTED AND UNFINISHED
            if nargin < 2 || isempty(filter)
                filter = '';
            end
            
            %% Run analysis on a subset if you specified a filter
            if isempty(filter) 
                filter = [vertcat(obj.expe_list{:})];
            end

            %% Regroup recordings per experiment
            data_folders_per_exp = {};
            count = 1;
            for idx = 1:numel(obj.expe_list)    
                if ismember(obj.expe_list(idx), filter)
                    data_folders_per_exp{count} = obj.rec_list{1}(obj.rec_list{2} == idx);
                    count = count+1;
                end
            end            
        end
        
        function expe = add_experiment(obj, path)  
            %% Identify if the experiment is already here to avoid duplicate
            expe = [];
            for el = 1:numel(obj.general_info)
                if isstruct(obj.general_info{el}) && isfield(obj.general_info{el}, 'short_path') && any(strcmp(obj.general_info{el}.short_path, path(1).name(:,1:end-9)))
                    expe = el;
                    break
                end                    
            end
            if isempty(expe)
                expe = numel(obj.general_info) + 1;
            end
            
            %% Prepare the fields for the indicated experiment
            obj.general_info{expe}      = [];
            obj.extracted_traces{expe}  = [];
            obj.binned_data{expe}       = [];
            obj.timescale{expe}         = [];
            obj.event_fitting{expe}     = [];
            obj.crosscorr{expe}         = [];
            obj.dimensionality{expe}    = [];
            obj.external_variables{expe}= [];
            %obj.expe_list{expe}         = [];

            obj.general_info{expe}.original_path = path;
            obj.general_info{expe}.short_path = path(1).name(:,1:end-9);
            
            try
                obj.cleanup_expe();
            catch
                    obj.cleanup_expe();
            end
            expe = find(cellfun(@(x) contains(x.short_path, path(1).name(:,1:end-9)), obj.general_info));
            obj.current_expe = expe;
        end
        
        function process(obj, expe, data_folders_per_exp, condition, filter_win, rendering)
            if nargin < 5 || isempty(filter_win)
                filter_win = [100,0];                
            end            
            if nargin >= 6 && ~isempty(rendering)
                obj.rendering = rendering;
            end
            demo = 0
            save_data = false
            
            obj.filter_win = filter_win;
            original_expe_idx = expe; % for the saving part 
            fprintf(['Processing experiment #',num2str(expe),'\n'])
            
            %% Prepare fields some useful information
            expe = obj.add_experiment(data_folders_per_exp{original_expe_idx}); 

            %% Load and concatenate traces for the selected experiment
            obj.load_extracted_data(expe);   % this also sets the current expe #
% 
%             %% Check for consistent number of ROIs
%             quality_control(1, obj.extracted_traces{expe})
%             %quality_control(2, cell2mat(all_sr))

            %% Prepare binning of ROIs based on specific grouping condition
            [bins, metrics, bin_legend] = obj.prepare_binning(condition);

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces();

            %% Create median trace per bins
            obj.set_median_traces()

            %% Get summary covariance plots
            obj.similarity_plot();            

            %% Correct for decay to avoid overstimating peak amplitude
            if isempty(obj.event_fitting{expe})
                obj.event_fitting{expe} = detect_and_fit_events(obj.binned_data{expe}.median_traces, obj.timescale{expe}.global_timescale, demo, obj.binned_data{expe}.bin_legend);arrangefigures([1,2]);
            end

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.plot_events_distribution();

            %% Study bAP heterogeneity
            obj.assess_variability()

            %% Check how correlated are the peaks between different parts of the tree
            if isempty(obj.crosscorr{expe})
                obj.crosscorr{expe} = corrcoef([nanmean(obj.event_fitting{expe}.post_correction_peaks, 2), obj.event_fitting{expe}.post_correction_peaks]);
            end
            obj.plot_cc(); 

            %% Plot a tree correlation metric     
            %ROIs_list  = unique([all_ROI_ids_per_bin{:}]);
            %[fixed_tree, soma_location, ROIs_list, listing, batch_params] = rebuild_tree(data_folders_per_exp{original_expe_idx}(1), setting_file_path, true);

            %% Assign value of group to these ROIs
            obj.plot_corr_tree(); 

            cross_validate = false;
            obj.get_dimensionality(cross_validate);

            %% Plot weight-tree for each component
            for comp = 1:5
                obj.plot_dim_tree(comp);
            end

            %% Plot a map of tree weights by ROI number for each component
            obj.plot_weight_map();

            %% Plot strongest componenet per ROI
            obj.plot_strongest_comp_tree();

            %% Optionally, if external variables need an update
            %obj.update_external_metrics(60)

            %% Save figure
            if save_data
                obj.save_figures()
            end

            %catch
            %   failed{expe} = data_folders_per_exp{original_expe_idx}; 
            %end
            %close all
            figure(999);plot(obj.binned_data{expe}.median_traces - nanmean(obj.binned_data{expe}.median_traces,2));set(0,'DefaultAxesColorOrder',magma(size(obj.binned_data{expe}.median_traces,2)))
        end
        
        function cleanup_expe(obj, deep_cleanup, updated_path)
            if nargin < 2 || isempty(deep_cleanup)
                deep_cleanup = false; source = '';
            end
            %% Remove empty location and sort by expe name
            to_keep = find(~cellfun(@isempty, obj.general_info));
            
            obj.general_info = obj.general_info(to_keep);
            obj.extracted_traces = obj.extracted_traces(to_keep);
            obj.binned_data = obj.binned_data(to_keep);
            obj.timescale = obj.timescale(to_keep);
            obj.event_fitting = obj.event_fitting(to_keep);
            obj.crosscorr = obj.crosscorr(to_keep);
            obj.dimensionality = obj.dimensionality(to_keep);
            obj.expe_list = obj.expe_list(to_keep);
            obj.need_update = obj.need_update(to_keep);
            obj.external_variables  = obj.external_variables(to_keep);
            
            [~,order] = sort(cellfun(@(x) x.short_path, obj.general_info, 'UniformOutput', false));
            
            obj.expe_list = obj.expe_list(order);
            obj.general_info = obj.general_info(order);
            obj.extracted_traces = obj.extracted_traces(order);
            obj.binned_data = obj.binned_data(order);
            obj.timescale = obj.timescale(order);
            obj.event_fitting = obj.event_fitting(order);
            obj.crosscorr = obj.crosscorr(order);
            obj.dimensionality = obj.dimensionality(order);
            obj.need_update = obj.need_update(order);
            obj.external_variables  = obj.external_variables(order);
            
            if strcmp(deep_cleanup, 'extracted') && isfolder(obj.source_folder)
                source = obj.source_folder;
                if nargin > 2 && ~isempty(updated_path)
                    source = updated_path;
                end
            elseif strcmp(deep_cleanup, 'original') && isfolder(obj.source_folder) && ~isempty(updated_path)
                source = obj.source_folder;
            elseif deep_cleanup
                source = '';
                warning('deep_cleanup must be "original", or "extracetd". If the extracted path changed or if you want to use the original data for cleaning, indicate a valid path in "updated_path" input')
            end
            
            if ~isempty(source)
                for exp = 1:numel(obj.expe_list)
                    if strcmp(deep_cleanup, 'original') 
                    	to_remove = cellfun(@(x) ~isdir(x), obj.general_info{exp}.sources);
                    elseif strcmp(deep_cleanup, 'extracted') 
                        to_remove = cellfun(@(x) ~isdir(x), arrayfun(@(y) [y.folder,'\',y.name], obj.general_info{exp}.original_path, 'UniformOutput', false))';
                    end
                    if any(to_remove)
                        obj.need_update(exp) = true;
                        
                        %% Remove unnecessary path
                        all_names = {obj.general_info{exp}.original_path(to_remove).name};
                        to_remove_too = find(cellfun(@(x) any(strcmp(x, all_names)), {obj.rec_list{1}.name}));
                        obj.rec_list{1}(to_remove_too) = [];obj.rec_list{2}(to_remove_too) = [];
                        obj.general_info{exp}.sources = obj.general_info{exp}.sources(~to_remove);
                        obj.general_info{exp}.original_path = obj.general_info{exp}.original_path(~to_remove); 
                    end                    
                end
            end
        end
        
        
        function load_extracted_data(obj, expe, tracing_source)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            if nargin < 3 || isempty(tracing_source)
                tracing_source = 'swc';
            end
            
            cleanup = true;
            
            data_folders_current_exp = obj.general_info{expe}.original_path;            
            extracted_traces        = {}; % N recordings cell array of M analysis bins cell array of [Timepoints x ROI] matrices
            all_sr                  = {}; % N FLOAT of sampling rates
            arboreal_scans          = {}; % N recordings cell array of analysis_params, without data

            %% Load extracted data
            unextracted   = [];
            invalid       = [];
            
            for data_f_idx = 1:numel(data_folders_current_exp)
                paths = parse_paths([data_folders_current_exp(data_f_idx).folder, '/', data_folders_current_exp(data_f_idx).name]);
                source              = dir([paths,'/**/*.mat']);
                if ~isempty(source) % if folder exist and has an extracted mat file
                    extracted_data                  = load([source.folder, '/', source.name]);
                    extracted_data                  = extracted_data.obj;
                    
                    source_still_exist              = isfolder(extracted_data.data_folder);
                    if ~cleanup || (cleanup && source_still_exist)
            %             %% Reselect ROI subset
            %             if strcmp(tracing_source, 'swc')
            %                params = analysis_params('source'          ,extracted_data.data_folder,...
            %                                         'tracing_source'  , 'swc'         );    
            %                [~, ROIs] = get_branch_id_from_tracing_source(params);
            %             end
                        extracted_traces{data_f_idx}                = extracted_data.data;
                        extracted_data.data                         = []; % clear data content
                        extracted_data.analysis_params.data         = [];
                        extracted_data.analysis_params.simple_data  = [];
                        arboreal_scans{data_f_idx}           = extracted_data;
                        all_sr{data_f_idx}              = diff(extracted_data.analysis_params.timescale(1:2));  
                    elseif cleanup && ~source_still_exist
                        invalid                       = [invalid, data_f_idx];
                    end
                else
                    unextracted                       = [unextracted, data_f_idx];
                end
            end
            
            if ~isempty(invalid) % if cleanup
                obj.general_info{expe}.original_path(invalid) = [];
                obj.rec_list{expe}(invalid) = [];
            end

            %% Get rid of any bad file
            if ~isempty(unextracted)
                error_box('Curation and/or data extraction not clean. Consider reviewing this folder')
            end

            arboreal_scans                                              = arboreal_scans(~cellfun(@isempty, arboreal_scans));
            obj.general_info{expe}.sources                              = arrayfun(@(x) x.data_folder, [arboreal_scans{:}], 'UniformOutput' ,false);    
            obj.general_info{expe}.arboreal_scan                        = arboreal_scans{1};    
            obj.extracted_traces{expe}                                  = extracted_traces(~cellfun(@isempty , extracted_traces));
            
            %% Prepare timescale for each recording and concatenated timescale
            obj.timescale{expe}.sr  = all_sr(~cellfun(@isempty , all_sr));
            obj.timescale{expe}.tp  = cellfun(@(x) size(x, 1), obj.extracted_traces{expe},'UniformOutput', false);
            trial_timescale         = cellfun(@(x, y) linspace(0, y*x, x), obj.timescale{expe}.tp, obj.timescale{expe}.sr, 'UniformOutput', false);
            global_timescale        = cellfun(@(x) diff(x), trial_timescale, 'UniformOutput', false);
            global_timescale        = cellfun(@(x) [x(1), x], global_timescale, 'UniformOutput', false);
            obj.timescale{expe}.global_timescale = cumsum(horzcat(global_timescale{:}));
        end
        
        function failed = update_external_metrics(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = 1:numel(obj.expe_list);
            end
            
            failed = {};
%             for el = 1:numel(obj.general_info)
%                 current_sources = obj.general_info{el}.sources;
%                 current_external_variables = {};
%                 for rec = 1:numel(current_sources)
%                     fprintf([current_sources{rec}, '\n']);
%                     current_external_variables{rec} = [];                 
%                 end  
%                 obj.external_variables{el} = current_external_variables;
%             end        
%             
            failed = {};
            for el = expe
                current_sources = obj.general_info{el}.sources;
                current_external_variables = {};
                for rec = 1:numel(current_sources)
                    fprintf([current_sources{rec}, '\n']);
                    try
                    parameters = analysis_params('source',current_sources{rec},'data_type','raw');
                    [parameters, ~, ~, ~] = load_generic_scan(parameters);
                    parameters.external_var = structfun(@(x) struct('time',x.time(:),'data',x.speed(:)), parameters.external_var,'UniformOutput', false);
                    current_external_variables{rec} = parameters.external_var;
                    catch err
                        current_external_variables{rec} = [];
                        failed{end + 1} = [{err}, current_sources{rec}];
                    end                    
                end  
                obj.external_variables{el} = current_external_variables;
            end            
        end
        
        
        function [bins, metrics, bin_legend] = prepare_binning(obj, condition, expe)
            if nargin < 3 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            %% Define current binning rule. See arboreal_scan.get_ROI_groups for more info    
            % CC and peak extractions are based on this binning            
            [bins, metrics, bin_legend]         = obj.general_info{expe}.arboreal_scan.get_ROI_groups(condition);
            obj.binned_data{expe}.groups        = bins;
            obj.binned_data{expe}.bin_legend    = bin_legend;
            obj.binned_data{expe}.metrics       = metrics;    
            obj.binned_data{expe}.readmap       = sort(unique([bins{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
            
            %% Clear dependant analyses
            obj.crosscorr{expe}                 = [];
            obj.event_fitting{expe}             = [];
            obj.general_info{expe}.offset       = zeros(1,size(obj.extracted_traces{expe}{1}, 2)); 
            obj.general_info{expe}.scaling      = ones(1,size(obj.extracted_traces{expe}{1}, 2)); 
        end
        
        function rescale_traces(obj, expe) 
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            
            [best_scal_f, best_offset]          = scale_every_recordings(obj.extracted_traces{expe}, obj.demo);
            obj.general_info{expe}.scaling      = best_scal_f;
            obj.general_info{expe}.offset       = best_offset;
            if obj.rendering
                figure();plot(nanmedian(obj.get_rescaled_traces(),2));
                arrangefigures([1,2]);
            end
        end
        
        function rescaled_traces = get_rescaled_traces(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            all_traces_per_rec = obj.extracted_traces{expe};
            
            %% Now rescale each trace with a unique value across recordings (which would be spectific of that region of the tree).
            for trace = 1:size(obj.extracted_traces{expe}{1}, 2) % [1    37    83   102]%
                temp = cellfun(@(x) x(:, trace) ,all_traces_per_rec, 'UniformOutput', false); % 'end' or any group you want to test
                for rec = 1:numel(temp)
                    temp{rec}(1:2) = NaN; % first 2 points are sometimes a bit lower than the rest. We NaN them. This is also useful later on to identify trials
                end
                concat_trace    = vertcat(temp{:}); off = prctile(concat_trace, obj.general_info{expe}.offset(trace));
                temp = cellfun(@(x) x - off, temp, 'UniformOutput', false);
                temp = cellfun(@(x) x / obj.general_info{expe}.scaling(trace), temp, 'UniformOutput', false);
                for rec = 1:numel(temp)
                    all_traces_per_rec{rec}(:, trace) = temp{rec};
                end
            end 
            rescaled_traces = vertcat(all_traces_per_rec{:});           
        end
        
        function set_median_traces(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            %obj.extracted_traces{expe}           = vertcat(obj.extracted_traces{expe}{:}); %% QQ set as a get method
            rescaled_traces = obj.get_rescaled_traces(expe);
            
            %rescaled_traces = rescaled_traces - prctile(nanmedian(wavelet_denoise(rescaled_traces),2),5);
            
            %% Create median trace per bins
            all_traces_per_bin = cell(1, numel(obj.binned_data{expe}.groups));
            for gp = 1:numel(obj.binned_data{expe}.groups)
                columns                     = ismember(obj.binned_data{expe}.readmap, obj.binned_data{expe}.groups{gp});
                all_traces_per_bin{gp}      = nanmedian(rescaled_traces(:,columns), 2);
            end    
            
            obj.binned_data{expe}.median_traces = cell2mat(all_traces_per_bin);%cellfun(@(x) cell2mat(x) ,all_traces_per_bin, 'UniformOutput', false); 
            
            if obj.rendering
                %% Plot the mean trace for each bin  
                figure(1001);cla();
                plot(obj.timescale{expe}.global_timescale, obj.binned_data{expe}.median_traces); hold on;
                legend(obj.binned_data{expe}.bin_legend); hold on;
                title('mean scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
                arrangefigures([1,2]); 
            end
        end
    
        function similarity_plot(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            win = ceil(1./median(diff(obj.timescale{expe}.global_timescale)));
            comb = nchoosek(1:size(obj.binned_data{expe}.median_traces,2),2);
            corr_results = {};
            for pair = 1:size(comb,1)
                corr_results{pair} = movcorr(obj.binned_data{expe}.median_traces(:,comb(pair, 1)),obj.binned_data{expe}.median_traces(:,comb(pair, 2)),[win, 0]);
            end

            corr_results = cell2mat(corr_results);
            precision = 1./nanvar(corr_results,[], 2);
            out = nanmax(precision(:)) / 10;
            precision(precision > out) = NaN;

            if obj.rendering
                f = figure(1006);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                f.Tag = 'Similarity plot'; %for figure saving
                ax1 = subplot(3,1,1); hold on;plot(obj.binned_data{expe}.median_traces);title('Cell Signal');ylim([-50,range(obj.binned_data{expe}.median_traces(:))]);
                ax2 = subplot(3,1,2); hold on;plot(corr_results,'Color',[0.9,0.9,0.9]);
                hold on;plot(nanmean(corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(corr_results, 1)],[-1,-1],'k--');plot([0,size(corr_results, 1)],[1,1],'k--')
                ax3 = subplot(3,1,3); hold on;plot(precision);title('Precision (1/Var)');set(ax3, 'YScale', 'log');plot([0,size(precision, 1)],[1,1],'k--')
                linkaxes([ax1,ax2, ax3],'x')
                arrangefigures([1,2]); 
            end
        end
        
        function expe = identify(obj, filter)
            expe = find(cellfun(@(x) contains(x, filter), obj.expe_list));
            fprintf(['Required experiment is # ', num2str(expe),'\n'])
            
        end
        
        function norm_cumsum = plot_events_distribution(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            peaks = obj.event_fitting{expe}.post_correction_peaks;
            
            %% Detect and display peak histogram distribution (mean and individual groups)
            max_peaks = max(peaks(:));
            mean_pks  = mean(peaks,2);

            bin_size = max_peaks/30;
            figure(1003);cla();title('peak mean distribution'); xlabel('Amplitude'); ylabel('counts');set(gcf,'Color','w');
            if ~isempty(mean_pks)
                histogram(mean_pks,0:bin_size:max_peaks, 'FaceColor', 'k');
            end
            f = figure(1004);clf();hold on; title('peak distribution per group');hold on;set(gcf,'Color','w');hold on;
            f.Tag = 'peak distribution per group'; %for figure saving
            cmap = lines(size(obj.binned_data{expe}.median_traces, 2));
            [m,n] = numSubplots(size(obj.binned_data{expe}.median_traces, 2));
            for gp = 1:size(obj.binned_data{expe}.median_traces, 2)
                subplot(m(1), m(2), gp)
                %figure();plot(global_timescale,obj.binned_data{expe}.median_traces(:,gp));hold on;scatter(obj.event_fitting{expe}.peak_times,obj.event_fitting{expe}.post_correction_peaks(:,gp), 'kv')
                if ~isempty(peaks)
                    histogram(peaks(:,gp),0:bin_size:max_peaks, 'FaceAlpha', 0.8,'EdgeColor','none', 'FaceColor', cmap(gp, :));hold on;
                end
                title(obj.binned_data{expe}.bin_legend(gp));
            end

            %% Plot all individual peaks to see if they covary
            figure(1005);cla();plot(peaks);legend(obj.binned_data{expe}.bin_legend);title('peaks values per subgroups'); xlabel('event #'); ylabel('Amplitude');set(gcf,'Color','w');
            
            %% Plot the mean amplitude per subgroup of peaks, per distance          
            %bin_step = 10;
            %norm_cumsum = cumsum(obj.binned_data{expe}.median_traces) ./ nanmax(cumsum(obj.binned_data{expe}.median_traces));
            norm_cumsum = cumsum(peaks) ./ nanmax(cumsum(peaks));
            figure(1007);cla();plot(norm_cumsum); hold on;set(gcf,'Color','w');xlabel('Event #');ylabel('normalized cumulative amplitude')
            title('cumulative sum of peaks'); ylim([0, 1]);legend(obj.binned_data{expe}.bin_legend,'Location','southeast');
        end
        
        function assess_variability(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            %sr = nanmedian(diff(obj.timescale{expe}.global_timescale));
            vmr = nanvar(obj.event_fitting{expe}.post_correction_peaks,[],2)./nanmean(obj.event_fitting{expe}.post_correction_peaks, 2);
            %cv  = nanstd(obj.event_fitting{expe}.post_correction_peaks,[],2)./nanmean(obj.event_fitting{expe}.post_correction_peaks, 2); % (maybe chack snr at one point?  mu / sigma)
            %fano = []; % windowed VMR. usually for spike trains
            [~, idx] = sort(vmr,'descend');
            
            %figure(123); cla();ylim([0, nanmax(obj.event_fitting{expe}.post_correction_peaks(:))]); hold on;
            %     for event = idx'
            %         show_event(obj.binned_data{expe}.median_traces, round(obj.event_fitting{expe}.peak_times/sr), event);
            %         drawnow;%pause(0.1)
            %     end

            figure(1009);cla();plot(obj.event_fitting{expe}.post_correction_peaks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')

            R = max(range(obj.binned_data{expe}.median_traces));
            %norm_vmr = vmr/range(vmr);
            norm_vmr = vmr/mean(vmr);
            figure(1010);clf();
            ax1 = subplot(2,1,1);plot(obj.timescale{expe}.global_timescale, obj.binned_data{expe}.median_traces); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
            ax2 = subplot(2,1,2);plot(obj.event_fitting{expe}.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
            plot([obj.timescale{expe}.global_timescale(1), obj.timescale{expe}.global_timescale(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
            plot([obj.timescale{expe}.global_timescale(1), obj.timescale{expe}.global_timescale(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
            hold on;plot([obj.timescale{expe}.global_timescale(1), obj.timescale{expe}.global_timescale(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');
            linkaxes([ax1, ax2], 'x');

            %figure();histogram(norm_vmr, round(10*range(norm_vmr)/std(norm_vmr)))

            figure(1011);cla();scatter(nanmedian(obj.event_fitting{expe}.post_correction_peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
        end
        
        function plot_cc(obj, expe)  
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            cc = obj.crosscorr{expe};
            bin_legend = obj.binned_data{expe}.bin_legend;
            
            if obj.rendering
                figure(1008);cla();imagesc(cc); hold on;set(gcf,'Color','w');
                caxis([0,1]); hold on;
                xticks([1:size(cc, 1)])
                yticks([1:size(cc, 1)])
                colorbar; hold on;
                plot([1.5,1.5],[0.5,size(cc, 1)+0.5],'k-');xticklabels(['whole cell',bin_legend]);xtickangle(45);
                plot([0.5,size(cc, 1)+0.5],[1.5,1.5],'k-');yticklabels(['whole cell',bin_legend]);
                title('peak amplitude correlation per subgroup');
                arrangefigures([1,2]); 
            end
        end
        
        function get_dimensionality(obj, cross_validate, expe)
            if nargin < 2 || isempty(cross_validate)
                cross_validate = false;
            end
            if nargin < 3 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            rescaled_traces     = obj.get_rescaled_traces(expe);
            all_ROIs            = 1:size(rescaled_traces, 2);
            normal_n_NaN        = median(sum(isnan(rescaled_traces))) * 4;
            valid_trace_idx     = sum(isnan(rescaled_traces)) <= normal_n_NaN;
            rescaled_traces     = fillmissing(rescaled_traces(:, valid_trace_idx),'spline'); % removed funny traces
            valid_ROIs          = all_ROIs(valid_trace_idx);

            %% Get single or multiple factor estimate
            if ~cross_validate
                [LoadingsPM, specVarPM, T, stats, F] = factoran(double(rescaled_traces), 10);
            else
                %% Cross validation (thanks Harsha)
                [~,N]              = size(rescaled_traces);
                nFactors                    = 40;%round(0.66*N);
                nCV                         = 5;
                train                       = NaN(nFactors, nCV);
                test                        = NaN(nFactors, nCV);
                fac_steps                   = 1;

                for jj = 1:nCV
                    %% Partition data into Xtrain and Xtest
                    test_idx = randperm(size(rescaled_traces, 1));
                    Xtrain = rescaled_traces(test_idx(1:2:end), :);   
                    Xtest  = rescaled_traces(test_idx(2:2:end), :);

                    [jj, toc]
                    parfor factor_nb = 1:nFactors  
                        if ~rem(factor_nb-1, fac_steps) % every fac_steps steps, starting at 1
                            [LoadingsPM, specVarPM, ~, stats, F]    = factoran(Xtrain, factor_nb);    
                            train(factor_nb, jj)                 = stats.loglike;
                            test(factor_nb, jj)                  = testloglike_factorAnalysis(Xtest, nanmean(Xtrain,1), LoadingsPM, specVarPM);
                        end
                    end
                end
                toc

                tested_modes = 1:fac_steps:nFactors;
                figure();plot(tested_modes,train(tested_modes,1:nCV),'o');hold on;plot(tested_modes,nanmean(train(tested_modes,1:nCV), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('train');
                figure();plot(tested_modes,test(tested_modes,1:nCV),'o');hold on;plot(tested_modes,nanmean(test(tested_modes,1:nCV), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('test');

                error('Add optimal selection here')
            end

            %% Store results
            obj.dimensionality{expe}.LoadingsPM         = LoadingsPM;
            obj.dimensionality{expe}.specVarPM          = specVarPM;
            obj.dimensionality{expe}.T                  = T;                % Rotation matrix
            obj.dimensionality{expe}.stats              = stats;            % Factoran stats
            obj.dimensionality{expe}.F                  = F;                % components
            obj.dimensionality{expe}.all_ROIs           = all_ROIs;         % first occurences
            obj.dimensionality{expe}.valid_trace_idx    = valid_trace_idx;  % additional filter for recordings with many NaNs
        end
        
        function f = plot_value_tree(obj, expe, values, locations, ~)
            %% e.g.
            % values = max(meta_analysis_result.extracted_traces{1}{1});
            % meta_analysis_result.plot_value_tree(expe,values 1:numel(values));

            tree_info       = obj.general_info{expe}.arboreal_scan.tree_info;
            global all_traces data_folders mask_name
            all_traces      = obj.get_rescaled_traces(expe);
            data_folders    = obj.general_info{expe}.sources;
            mask_name       = obj.general_info{expe}.arboreal_scan.analysis_params.mask_value;
            %f_handle        = @(x) load_several_experiments(x, data_folders);
            f_handle        = @(x) plot_ROI(x);
            path = tree_info.path;
            path = strrep(path,'D:/V1_Curated_Data/' ,'Y:\general\Analysis\Antoine_and_Tommy\V1_Curated_Data\');
            f               = plot_summary_tree(path, values, locations, '', tree_info.soma_location, '', tree_info.excluded_branches, tree_info.pia_soma,1018, f_handle); set(gcf,'Color','w');arrangefigures([1,2]);
        end

        function [tree, soma_location, tree_values, mean_bin_cc] = plot_corr_tree(obj, expe)  
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            %% Identify valid set of traces
            valid_gp            = find(~all(isnan(obj.binned_data{expe}.median_traces))); % QQ not sure when we have full NaN series, but we can not denoise them with wdenoise
            if valid_gp(1) == 1
                valid_gp = valid_gp(2:end);
            end
            
            %% Reload CC
            cc              = obj.crosscorr{expe}(2:end,2:end);
            
            %% Build tree values per bin
            mean_bin_cc     = [];
            ROIs_list       = [];
            if ~isempty(cc)
                for gp = 1:numel(obj.binned_data{expe}.groups)
                    roi_of_gp           = obj.binned_data{expe}.groups{gp};
                    v_of_gp             = cc(valid_gp(1),gp);
                    ROIs_list           = [ROIs_list, roi_of_gp];
                    mean_bin_cc         = [mean_bin_cc, repmat(v_of_gp, 1, numel(roi_of_gp))];
                end
            end

            %% Map CC values on the tree
            f                       = obj.plot_value_tree(expe, mean_bin_cc, ROIs_list);caxis([0,1]);title('Correlation with most proximal segment');
            
            %% Prepare output
            tree                    = obj.general_info{expe}.arboreal_scan.original_tree;
            soma_location           = obj.general_info{expe}.arboreal_scan.soma_location;
            tree_values             = {f.XData; f.YData; f.ZData; f.CData};
        end
        
        function [tree, soma_location, tree_values, values] = plot_dim_tree(obj, comp, expe)
            if nargin < 3 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            if nargin < 2 || isempty(comp) || ~comp
                [tree, soma_location, tree_values, values] = obj.plot_strongest_comp_tree(expe);
                return
            end
            
            % check 58, % noise issue 64 'D:/Curated Data/2019-09-24/experiment_1/18-13-20/'

            %% Reload loadings
            LoadingsPM = obj.dimensionality{expe}.LoadingsPM;
            Valid_ROIs = obj.dimensionality{expe}.all_ROIs(obj.dimensionality{expe}.valid_trace_idx);

            values = NaN(1, numel(Valid_ROIs));
            if comp
                for roi = 1:numel(Valid_ROIs)
                    %ROI = valid_ROIs(roi);
                    values(roi) = LoadingsPM(roi, comp); 
                end
                titl = ['Weighted average for component ',num2str(comp),' per ROI'];
            else
                [~, loc]    = max(LoadingsPM(:,1:5)');
                for ROI = 1:numel(Valid_ROIs)
                    roi = Valid_ROIs(ROI);
                    values(roi) = loc(ROI);
                end 
                titl = 'Location of strongest component';
            end
            
            %% Map dimension weights on the tree
            f                       = obj.plot_value_tree(expe, values, Valid_ROIs);title(titl);            
            tree                    = obj.general_info{expe}.arboreal_scan.original_tree;
            soma_location           = obj.general_info{expe}.arboreal_scan.soma_location;
            tree_values             = {f.XData; f.YData; f.ZData; f.CData};
        end
        
        function [tree, soma_location, tree_values, values] = plot_strongest_comp_tree(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            [~, loc]    = max(obj.dimensionality{expe}.LoadingsPM(:,1:5)');
            values      = NaN(size(loc));
            Valid_ROIs  = obj.dimensionality{expe}.all_ROIs(obj.dimensionality{expe}.valid_trace_idx);
            for ROI = 1:numel(Valid_ROIs)
                roi     = obj.dimensionality{expe}.all_ROIs(Valid_ROIs(ROI));
                values(roi)  = loc(ROI);
            end
            f                       = obj.plot_value_tree(expe, values, Valid_ROIs);   title('Location of strongest component')             
            
            %% Map dimension weights on the tree
            tree                    = obj.general_info{expe}.arboreal_scan.original_tree;
            soma_location           = obj.general_info{expe}.arboreal_scan.soma_location;
            tree_values             = {f.XData; f.YData; f.ZData; f.CData};
            arrangefigures([1,2]);
        end
        
        function [tree, soma_location, tree_values, values] = plot_dist_tree(obj, bin, expe)   
            if nargin < 2 || isempty(bin)
                bin = [];
            end
            if nargin < 3 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end


            %% Build tree values per bin
            %tree = obj.general_info{expe}.arboreal_scan.original_tree{1};
            tree                    = skeleton_to_tree(obj.general_info{expe}.arboreal_scan.analysis_params.trees{1}); % entire skeleton as tree, no xcludion
            [~, info]               = get_branch_id_from_ROI(tree, 0, 0);            
            morpho                  = obj.general_info{expe}.arboreal_scan.tree_info;
            tree                    = fix_tree(tree{1}, morpho.primary_dendrites, morpho.manual_reconnections, morpho.excluded_branches, morpho.pia_soma, true);
            info                    = info(~ismember(info(:,1),morpho.excluded_branches'),:); % QQ INVERSIONS NOT HANDELD
            ROIs_list               = info(:, 2);
            v                       = Pvec_tree(tree{1});
            br_offset               = 0;
            for roi = 1:numel(ROIs_list)
                v_idx = roi + br_offset;
                if roi == numel(ROIs_list) || ~diff(info(roi:(roi+1), 1))
                    values(roi) = mean(v(v_idx:(v_idx+1)));
                else
                    values(roi) = mean(v(v_idx:(v_idx+1)));
                    br_offset = br_offset + 1;
                end
            end

            %% Map distance values on the tree
            f                       = obj.plot_value_tree(expe, values, ROIs_list);title('Distance from soma');
            
            %% Prepare output
            tree                    = obj.general_info{expe}.arboreal_scan.original_tree;
            soma_location           = obj.general_info{expe}.arboreal_scan.soma_location;
            tree_values             = {f.XData; f.YData; f.ZData; f.CData}; 
        end    
        
        function plot_weight_map(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end
            
            %% Recover traces
            rescaled_traces = obj.get_rescaled_traces(expe);
            
            %% Reload loadings
            LoadingsPM = obj.dimensionality{expe}.LoadingsPM;
            Valid_ROIs = obj.dimensionality{expe}.all_ROIs(obj.dimensionality{expe}.valid_trace_idx);
            rescaled_traces = rescaled_traces(:, Valid_ROIs);

            
            figure(1017);cla();imagesc(LoadingsPM(:,1:5)); colorbar;set(gcf,'Color','w');title('Components 1 to 5 per ROI');xlabel('component');ylabel('ROI');
            [~, loc] = max(LoadingsPM(:,1:5)');

            w1 = LoadingsPM(:,1)/sum(LoadingsPM(:,1));
            w2 = LoadingsPM(:,2)/sum(LoadingsPM(:,2));
            w3 = LoadingsPM(:,3)/sum(LoadingsPM(:,3));
            w4 = LoadingsPM(:,4)/sum(LoadingsPM(:,4));
            w5 = LoadingsPM(:,5)/sum(LoadingsPM(:,5));
            figure(1021);cla();
            hold on; plot(nanmean(rescaled_traces'.* w1, 1));
            hold on; plot(nanmean(rescaled_traces'.* w2, 1));
            hold on; plot(nanmean(rescaled_traces'.* w3, 1));
            hold on; plot(nanmean(rescaled_traces'.* w4, 1));
            hold on; plot(nanmean(rescaled_traces'.* w5, 1));
            legend();set(gcf,'Color','w');
            title('Weighted signal average per component')
        end
        
        function plot_gallery(obj, mode, range, proj_axis)
            if nargin < 2 || isempty(mode)
                mode = 'corr';
            end
            if nargin < 3 || isempty(range)
                range = find(~cellfun(@isempty, obj.dimensionality));
            end
            if nargin < 4 || isempty(proj_axis)
                proj_axis = 'xz';
            end
            
            %% Now generate the figures
            [all_trees,  all_soma, all_values]  = deal({});
            for expe = range
                if strcmp(mode, 'corr')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.plot_corr_tree(expe);
                elseif strcmp(mode, 'dim')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.plot_dim_tree(0, expe); 
                elseif strcmp(mode, 'dist')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.plot_dist_tree(0, expe); 
                elseif strcmp(mode, 'order')
                    [all_trees{expe}, all_soma{expe}, all_values{expe}] = obj.plot_ord_tree(0, expe); 
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
                    shifted_x_coor = all_values{expe}{req_x}(1,:);
                    shifted_z_coor = all_values{expe}{req_y}(1,:);
                    
                    zero_shifted_coor = -min(shifted_x_coor) + x_offset;
                    shifted_x_coor = shifted_x_coor + zero_shifted_coor;
                    X = shifted_x_coor;
                    Y = shifted_z_coor+y_offset;
                    Z = all_values{expe}{req_z}(1,:);
                    S = all_values{expe}{4}(1,:);
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
        
        function save_figures(obj, expe)
            if nargin < 2 || isempty(expe)
                expe = obj.current_expe;
            else
                obj.current_expe = expe;
            end

            p = get(groot,'DefaultFigurePosition');
            folder = [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/'];
            if isfolder(folder)
                rmdir(folder,'s');
            end
            mkdir(folder);
            for fig_idx = [1001:1021, 1081, 1082, (10020+1):(10020+size(obj.event_fitting{expe}.events, 1)), (20020+1):(20020+5)]
                f = figure(fig_idx);
                set(f, 'Position', p)
                try
                    saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                    saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                    savefig(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
                catch % for figures with subplots
                    try
                        saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                        saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                        savefig(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
                    end
                end
            end
            name = parse_paths([folder,'/',f.Tag,data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.mat']);
            batch_params = obj.general_info{expe}.arboreal_scan.batch_params;

            %% Save every time in case of failure.
            save(parse_paths([obj.export_folder, 'summary.mat']), 'obj','-v7.3')
            %save(name, 'obj.binned_data{expe}.median_traces','all_fits', 'all_pks', 'peak_times','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','-v7.3');   
            %save(name, 'all_traces_per_bin','all_traces_per_bin','all_traces_concat','all_fits', 'all_pks', 'peak_times','cc','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','ROIs_list','dimensionality','-v7.3');   
        end
        
        function extracted_traces = get.extracted_traces(obj)
            extracted_traces = obj.extracted_traces;
            if ~isempty(obj.current_expe) && ~isempty(extracted_traces{obj.current_expe}) && any(obj.filter_win)
                temp = extracted_traces{obj.current_expe};
                extracted_traces{obj.current_expe} = cellfun(@(x) smoothdata(x, 1, obj.filter_type, obj.filter_win), temp, 'UniformOutput', false);
            end
        end
    end
end

