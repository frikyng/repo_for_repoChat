%% TO DO :
% - improve update system when changing folder.
% - enable loading if we just have a folder with arboreal_scans objects
% - add warning if source folde ris the actual raw experiment

classdef arboreal_scan_experiment < handle & arboreal_scan_plotting & event_fitting
    properties
        %% Extraction/Re-extraction Settings
        extraction_method = 'median'; % the 1D -> 0D compression method
        source_folder           % The folder where individual arboreal_scans were located when you built the object
        update_folder     = ''; % If you need to update the arboral_scans, but the orginal path changed, set the folder containing the new files here
        extracted_data_paths    % The original location of the individual arboreal_scans when you built the object
        need_update

        %% Children
        arboreal_scans          % A copy of the individual arboreal_scans, but where the uncompressed signal per ROIs had ben deleted

        %% Across expe Timescale info
        timescale               % time and sampling rate info for the individual recordings

        %% Saving options
        demo            = 0;
        auto_save_analysis = false
        auto_save_figures = false

        %% Analysis/extraction settings
        filter_win      = [0, 0];
        filter_type     = 'gaussian';
        dim_red_type    = 'nnmf'
        peak_thr        = 2;
        bad_ROI_thr     = 0.1;
        cc_mode         = 'groups_peaks'; % or raw
        detrend         = false;
        is_rescaled     = false; % is set to True once you ran the rescaling step.
        default_handle
        rescaling_method= 'by_trials';
        breakpoints     = []; % if you had a disruptive event during th experiment, a first scaling is done with large blocks

        %% All the fields computed in
        binned_data             % Defines how ROIs are grouped
        rescaling_info          % Defines how each ROi get rescaled to match the cell median
        event                   % Event detection output (event times, amplitude, correlation etc....)
        variability             % Signal variability over time
        dimensionality          % Results of diemensionality reduction
        behaviours              % List of available behaviours
        spiketrains             % If available, spike inference results
        bad_ROI_list     = 'unset';  % list of uncorrelated ROIs (following event detection)
        updatable               % If arboreal_scan are still available, you could update the arboral_scan_experiment compression
        crosscorr               % Correlation of peaks/signal across ROIs/groups during/between bAps/activity_bouts
    end

    properties (Dependent = true, Transient = true)
        updated_path            % If update_folder is used, the updated filpath
        extracted_traces        % Concatenated version of each obj.arboral_scan.simple_data
        extracted_pop           % Concatenated version of each obj.arboral_scan.simple_pop_data
        extracted_traces_conc   % Concatenated version of extracted_traces
        rescaled_traces         % Rescaled Traces according to rescaling_info
        extracted_pop_conc      % Concatenated version of extracted_pop
        global_median_raw       % The median of extracted_traces_conc
        t                       % Pointer to obj.timescale.global_timescale
        n_ROIs                  % Total number of ROIs in the swc, including bad ones
        n_pop_ROIs              % Total number of population ROIs
        ref                     % A pointer to the first extracted arboreal scan, for conveniency
        batch_params            % Pointer to obj.ref.batch_params ; the info to rebuild and locate the tree
        logs                    % Display the logs for all experiments
    end

    properties (Dependent = true, Transient = true, Hidden = true)
        external_variables      % Pointer to behavioural variables of each arboreal scan --> set in obj.behaviours
    end

    methods
        function obj = arboreal_scan_experiment(source_folder, varargin)
            if nargin < 1
                return
            end
            %% Fix paths
            obj.source_folder       = parse_paths(source_folder);

            %% Get a list of extracted arboreal scans in the source folder
            obj.extracted_data_paths= list_sources(obj);
            obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update

            %% Load arboreal scans
            obj.update(true, varargin);
        end

        function extracted_data_paths = list_sources(obj)
            %% List available arboreal scans
            all_recordings          = dir([obj.source_folder,'/**/*-*-*_exp_*_*-*-*']);
            if isempty(all_recordings)
                obj.updatable       = false;
                warning('no extracted arboreal_scans found in this folder. extraction/re-extraction not available');
                return
            else
                obj.updatable       = true;
            end
            all_recordings          = all_recordings(~[all_recordings(:).isdir]);
            all_recordings          = all_recordings(~(arrayfun(@(x) strcmp(x.name, '.'), all_recordings) | arrayfun(@(x) strcmp(x.name, '..'), all_recordings)));

            names                   = [vertcat(all_recordings.folder), repmat('/',numel(all_recordings), 1), vertcat(all_recordings.name)];
            extracted_data_paths    = cellfun(@(x) parse_paths(x), cellstr(names)', 'UniformOutput', false);
        end

        function update(obj, bypass, varargin)
            if nargin < 2 || isempty(bypass) || ~bypass
                quest = questdlg('WARNING : UPDATING SOURCES WILL DELETE ALL PROCESS DATA. Continue?','Update?','Yes','No','No');
            else
                quest = 'Yes';
            end

            if strcmp(quest, 'Yes')
                obj.extracted_data_paths= list_sources(obj);
                if isempty(obj.extracted_data_paths)
                    quest = questdlg('Do you want to try to extract arboreal_scans? If yes, you will be able to select the export folder in the next step','Extract?','Yes','No','No');
                    if strcmp(quest, 'Yes')
                        fold = parse_paths(uigetdir(pwd, 'Export folder'));
                        optional_analysis_params = [];optional_settings_path = [];
                        if ~isempty(varargin{1}) %  when extracting directly
                            optional_analysis_params = analysis_params(varargin{1}{1});
                            if numel(varargin{1}) > 1
                                optional_settings_path = varargin{1}{2};
                            end
                        end

                        err = meta_batch_process_ribbon_scan(obj.source_folder,optional_settings_path, optional_analysis_params,fold);
                        if isempty(err{1})
                            obj.source_folder       = fold;
                            obj.update(true);
                        else
                            error('Error detected during extraction. Check that the settings.txt file is present in the top_folder or manually indicated, and that it contains the correct paths')
                        end
                    else
                        return
                    end
                end
                obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update
                for field = {'arboreal_scans','binned_data', 'rescaling_info','event', 'variability','dimensionality'}
                    obj.(field{1}) = {};
                end

%                 %% Update arboreal_scan objects. Identify recordings that changed
%                 obj.extracted_data_paths = cellfun(@(x) parse_paths(x), obj.extracted_data_paths, 'UniformOutput', false); %temporary
%                 new_list = obj.list_sources();
%                 old_list_clean = erase(obj.extracted_data_paths, find_common_path(obj.extracted_data_paths));
%                 new_list_clean = erase(new_list, parse_paths(obj.source_folder));
%                 added          = ~ismember(new_list_clean, old_list_clean);
%                 deleted        = ~ismember(old_list_clean, new_list_clean);
%
%                 %% Clear deleted trees
%                 obj.arboreal_scans(deleted)         = [];
%                 obj.need_update(deleted)            = [];
%                 obj.extracted_data_paths(deleted)   = [];

                %% Update required trees first (bc idx are from before the detection of new trees)
                for el = fliplr(find(obj.need_update))
                    add_tree(el);
                end
%
%                 %% Add new trees
%                 for el = new_list(added)
%                     obj.extracted_data_paths{end+1} = el{1};
%                     add_tree(numel(obj.extracted_data_paths));
%                 end
%
%                 %% Sort trees by name
%                 [obj.extracted_data_paths, new_order] = sort(obj.extracted_data_paths);
%                 obj.arboreal_scans = obj.arboreal_scans(new_order);
            end

            function add_tree(pos)
                obj.arboreal_scans{pos}                             = load(obj.extracted_data_paths{pos});
                if isa(obj.arboreal_scans{pos}.obj, 'arboreal_scan')
                    obj.arboreal_scans{pos}                         = obj.arboreal_scans{pos}.obj;
                    obj.arboreal_scans{pos}.full_data               = []; % clear full data.
                    obj.arboreal_scans{pos}.population_data         = []; % clear population_data.
                    obj.need_update(pos)                            = false;
                else
                    obj.arboreal_scans(pos)                         = [];
                    obj.need_update(pos)                            = [];
                end
            end
        end

        %% ########### GET METHODS ###################

        function updatable = get.updatable(obj)
            updatable = cellfun(@(x) ~isempty(obj.updated_path), obj.extracted_data_paths);
        end

        function updated_path = get.updated_path(obj)
            if isempty(obj.update_folder)
                updated_path = obj.extracted_data_paths;
            else
                updated_path = cellfun(@(x) get_fname(x), obj.extracted_data_paths, 'uni', false);
            end

            function fname = get_fname(in)
                [~, in_name, in_ext] = fileparts(in);
                fname = [in_name,in_ext];
                replacement = dir([obj.update_folder,'/**/',fname]);
                if ~isempty(replacement)
                    fname = [replacement(1).folder,'/',replacement(1).name];
                end
            end
        end

        function breakpoints = get.breakpoints(obj)
            if isfield(obj.batch_params, 'breakpoints') && ~isempty(obj.batch_params.breakpoints)
                breakpoints = obj.batch_params.breakpoints;%find(cellfun(@(x) contains(x, obj.batch_params.breakpoints),obj.extracted_data_paths));
            else
                breakpoints = [];
            end
        end
        
        function set.breakpoints(obj, breakpoints)
            for rec = 1:numel(obj.arboreal_scans)
                obj.arboreal_scans{rec}.batch_params.breakpoints = breakpoints;
            end
        end
        
        function extracted_traces = fix_changes_in_gain(obj,extracted_traces)
        	temp        = vertcat(extracted_traces{:}); % get traces
            real_low    = 0;%nanmin(temp(:)); % real dark noise level (won't change if gain change)
            temp        = temp - real_low; % remove dark noise
            tp_start    = cumsum([1, obj.timescale.tp]); % pt start of trials
            breakpoints_idx = find(cellfun(@(x) ~isempty(obj.breakpoints) && contains(x, obj.breakpoints),obj.extracted_data_paths));
            blocks      = [1, breakpoints_idx, numel(extracted_traces)+1]; % pt range of the blocks
            pts         = [1,tp_start(breakpoints_idx),size(temp,1)+1]; % pt start of blocks
            bsl_value   = [];
            slope_gain  = {};
            
            for block_idx = 1:(numel(pts)-1)
                tmp_traces = temp(pts(block_idx):pts(block_idx+1)-1,:);
                
                %% Fix slop first since we'll fix the gain later
                local_tp        = size(tmp_traces,1);

                if obj.detrend
                    if any(isnan(tmp_traces(:)))
                        tmp_traces  = fillmissing(tmp_traces,'movmean',5);
                        to_use      = ~all(isnan(tmp_traces), 1);
                    else
                        to_use      = true(1,size(tmp_traces, 2));
                    end
                    
%                     %% Get slope for current trace
%                     p = polyfit(repmat(1:local_tp,sum(to_use),1),movmin(tmp_traces(:,to_use),50)',1);
%                     slope = (local_tp*p(1)) ./ nanmean(tmp_traces(:)); % slope in percent
%                     slope_gain{block_idx} = linspace(-slope/2,slope/2,local_tp); 
                    
%                     x = 1:local_tp;
%                     a = nanmean(movmin(tmp_traces(:,to_use),50),2);
%                     out = fit(x',a,'a + b*log(c - x)','Lower',[0,1,numel(a)],'Upper',[max(a),100,numel(a)*2]);toc
%                     slope_gain{block_idx} = out.a + out.b*log(out.c-x);
%                     %slope_gain{block_idx}(slope_gain{block_idx} < nanmin(a)) = nanmin(a);
%                     slope_gain{block_idx} = slope_gain{block_idx} ./ mean(slope_gain{block_idx}) - 1;;

                    if islogical(obj.detrend)
                        poly_order = 1;
                    else
                        poly_order = obj.detrend;
                    end
                    x = 1:local_tp;
                    p = polyfit(repmat(x,sum(to_use),1),movmin(tmp_traces(:,to_use),[300, 0])',poly_order);
                   % figure(555);cla();plot(nanmean(movmin(tmp_traces(:,to_use),300)')); hold on;plot(x, polyval(p,x))
                    slope_gain{block_idx} = polyval(p,x);
                    slope_gain{block_idx} = slope_gain{block_idx} ./ nanmean(slope_gain{block_idx}) - 1;
                    
                    
%                     out = fit(x',a,'a + b*log(c - x)','Lower',[0,1,numel(a)],'Upper',[max(a),100,numel(a)*2]);toc
%                     slope_gain{block_idx} = out.a + out.b*log(out.c-x);
%                     %slope_gain{block_idx}(slope_gain{block_idx} < nanmin(a)) = nanmin(a);
%                     slope_gain{block_idx} = slope_gain{block_idx} ./ mean(slope_gain{block_idx}) - 1;;


                    %% Temporarily corect the block
                   
                    tmp_traces = tmp_traces - tmp_traces.* slope_gain{block_idx}';

                    %% Reslice per trial
                    slope_gain{block_idx} = mat2cell(slope_gain{block_idx},1, obj.timescale.tp(blocks(block_idx):(blocks(block_idx+1)-1)));
                end
                if ~isempty(breakpoints_idx)
                    bsl_value = [bsl_value; prctile(tmp_traces,1)];
                end
            end
            
            %% Get scaling if gain fix is required
            if ~isempty(breakpoints_idx)
                scaling     = bsl_value .\ bsl_value(1,:); % normalize to first group;
            end
            
            %% Now apply correction per trial
            for block_idx = 1:(numel(pts)-1)
                trial_idx                     = (blocks(block_idx)):(blocks(block_idx+1)-1);
                if obj.detrend
                    extracted_traces(trial_idx)   = cellfun(@(x, y) x - x.*y', extracted_traces(trial_idx),slope_gain{block_idx}, 'uni',false);
                end
                if ~isempty(breakpoints_idx)
                    extracted_traces(trial_idx)   = cellfun(@(x) (x.*scaling(block_idx,:)), extracted_traces(trial_idx), 'uni',false);
                end
            end
        end
        
        
        function extracted_traces = get.extracted_traces(obj) % checked
            tic
            extracted_traces = cellfun(@(x) x.simple_data, obj.arboreal_scans, 'UniformOutput', false);
            
            %% If expe was interrupted signal gain changed, we fix it here
            if ~isempty(obj.breakpoints) || obj.detrend
            	extracted_traces = obj.fix_changes_in_gain(extracted_traces);
            end

            %extracted_traces = cellfun(@(x, y) x - y, extracted_traces,obj.rescaling_info.individual_offset, 'UniformOutput', false);
            
            if any(obj.filter_win)
                extracted_traces = cellfun(@(x) smoothdata(x, 'gaussian', obj.filter_win), extracted_traces, 'UniformOutput', false);
            end
            toc
        end

        function extracted_pop = get.extracted_pop(obj) % checked
            extracted_pop = cellfun(@(x) x.simple_pop_data, obj.arboreal_scans, 'UniformOutput', false);
            if obj.detrend
                extracted_pop = cellfun(@(x) x - prctile(x, 1), extracted_pop, 'UniformOutput', false);
            end
        end

        function extracted_traces_conc = get.extracted_traces_conc(obj) % checked
            extracted_traces_conc = vertcat(obj.extracted_traces{:});
        end

        function extracted_pop_conc = get.extracted_pop_conc(obj) % checked
            extracted_pop_conc = vertcat(obj.extracted_pop{:});
            % figure(666);cla();plot(normalize_sig(smoothdata(extracted_pop_conc,'gaussian',[100,0])', '', 'norm_method','dF/F0','percentile',10)')
        end

        function timescale = get.timescale(obj)
            %% Prepare timescale for each recording and concatenated timescale
            timescale                   = {};
            timescale.sr                = 1./cellfun(@(x) x.analysis_params.points_per_s, obj.arboreal_scans);
            timescale.time_source       = cellfun(@(x) x.analysis_params.time_source, obj.arboreal_scans, 'UniformOutput', false);
            if any(~strcmp(timescale.time_source, 'encoder'))
                warning('some timescale are estimated and not measured');
            end
            timescale.tp                = cellfun(@(x) x.analysis_params.timepoints, obj.arboreal_scans);
            timescale.durations         = timescale.sr.*timescale.tp; % same as cellfun(@(x) x.analysis_params.duration, obj.arboreal_scans)
            timescale.rec_timescale     = arrayfun(@(x, y) linspace(0, y*x, x), timescale.tp, timescale.sr, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) diff(x), timescale.rec_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) [x(1), x], timescale.global_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cumsum(horzcat(timescale.global_timescale{:}));
            timescale.t_start_nogap     = cellfun(@(x) x(end), timescale.rec_timescale);
            timescale.t_start_nogap     = cumsum([0, timescale.t_start_nogap(1:end-1) + timescale.sr(1:end-1)]);
            timescale.t_start_real      = [];  % QQ to do
        end

        function t = get.t(obj) % checked
            t = obj.timescale.global_timescale;
        end

        function ref = get.ref(obj)  % checked
            ref = obj.arboreal_scans{1};
        end

        function batch_params = get.batch_params(obj) % checked
            batch_params = obj.ref.batch_params;
            if ~isfield(batch_params, 'breakpoints')
                batch_params.breakpoints = [];
                warning('breakpoint field added post hoc. plase re-extract')
            end
        end

        function n_ROIs = get.n_ROIs(obj) % checked
            n_ROIs = size(obj.ref.simple_data,2);
        end

        function n_ROIs = get.n_pop_ROIs(obj) % checked
            n_ROIs = size(obj.ref.population_data,2);
        end
        
        function logs = get.logs(obj) % checked
            logs = cellfun(@(x) x.log, obj.arboreal_scans, 'Uniformoutput', false); % vertcat(obj.logs{:})
        end

        function f_handle = get.default_handle(obj)
            use_mask = false;
            f_handle = @(x) load_several_experiments(x, cellfun(@(x) x.data_folder, obj.arboreal_scans, 'UniformOutput', false), use_mask);
        end

        function global_median_raw = get.global_median_raw(obj) % checked
            global_median_raw = obj.extracted_traces_conc;
            global_median_raw(:, obj.bad_ROI_list) = [];
            global_median_raw = nanmedian(global_median_raw, 2);
        end

        function binned_data = get.binned_data(obj)
            binned_data = obj.binned_data;
%             if isfield(obj.binned_data, 'median_traces')
%             	binned_data.median_traces = smoothdata(binned_data.median_traces,'gaussian',obj.filter_win);
%             end
%             if isfield(obj.binned_data, 'global_median')
%             	binned_data.global_median = smoothdata(binned_data.global_median,'gaussian',obj.filter_win);
%             end
        end

        function external_variables = get.external_variables(obj)
            %failed_encoder = cellfun(@(x) ~numel(x.analysis_params.external_var.encoder.time), obj.arboreal_scans);
            external_variables = cellfun(@(x) x.analysis_params.external_var, obj.arboreal_scans, 'UniformOutput', false);
        end

        function behaviours = get.behaviours(obj)
            behaviours                  = obj.behaviours;
            behaviours.external_var     = obj.external_variables;
            behaviours.valid_encoder    = cellfun(@(x) ~isempty(x.encoder.time), behaviours.external_var);
            behaviours.valid_mc_log     = cellfun(@(x) isfield(x,'RT3D_MC'), behaviours.external_var);
            [Max_var, Max_var_loc]      = max(cellfun(@(x) numel(fieldnames(x)), behaviours.external_var));
            behaviours.types            = fieldnames(behaviours.external_var{Max_var_loc})';
            behaviours.valid_behaviours = false(numel(behaviours.valid_encoder), Max_var);
            for beh = 1:Max_var
                behaviours.valid_behaviours(:,beh) = arrayfun(@(y) isfield(y, behaviours.types{beh}), behaviours.external_var)';
            end
        end
        
        function bad_ROI_list = get.bad_ROI_list(obj)
            bad_ROI_list = obj.bad_ROI_list;
            if ischar(bad_ROI_list) && strcmpi(bad_ROI_list, 'unset')
                warning('Uncorrelated ROIs were never set. Detecting now. To ignore this set value to []')
                %rendering = obj.rendering;obj.rendering = false;
                obj.detect_events(); %find uncorrelated ROIs based on correlation
                %obj.rendering = rendering;
                bad_ROI_list = obj.bad_ROI_list;
            end
        end
        

        %% ########### ... ###################

        function [raw_beh, downsampd_beh, concat_downsamp_beh] = get_behaviours(obj, type, rendering, detrend_sig)
            %% Return the selected behaviour and corresponding timescale
            if nargin < 2 || isempty(type)
                type = 'encoder';
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end
            %if nargin < 4 || isempty(detrend_sig)
                detrend_sig = false;
            %end            

            downsampd_beh               = {};
            concat_downsamp_beh.time    = [];
            concat_downsamp_beh.value   = [];
            temp = {};
            type                        = obj.behaviours.types(contains(obj.behaviours.types, type));
            if isempty(type)
                raw_beh = {};
                warning(['Type not detected. Valid behaviours are :\n', strjoin(obj.behaviours.types,'\n')]);
                return
            end
            for beh = 1:numel(type)
                raw_beh{beh}                     = arrayfun(@(x) x{1}.(type{beh}), obj.behaviours.external_var, 'UniformOutput', false, 'ErrorHandler', @cellerror_empty);

                temp.time  = [];
                temp.value = [];
                for rec = 1:numel(raw_beh{beh})                    
                    %% If behaviour timescale is longer than recording, clip it
                    if ~isempty(raw_beh{beh}{rec}.time)
                        clipping_idx = find(raw_beh{beh}{rec}.time - obj.timescale.durations(rec) > 0, 1, 'first');
                        if ~isempty(clipping_idx)
                            raw_beh{beh}{rec}.time = raw_beh{beh}{rec}.time(1:clipping_idx);
                            raw_beh{beh}{rec}.value = raw_beh{beh}{rec}.value(1:clipping_idx);
                        end
                    end
                        
                    downsampd_beh{beh}{rec}.time     = interpolate_to(raw_beh{beh}{rec}.time, obj.timescale.tp(rec));
                    downsampd_beh{beh}{rec}.value    = interpolate_to(raw_beh{beh}{rec}.value, obj.timescale.tp(rec));
                    if isempty(downsampd_beh{beh}{rec}.time)
                        downsampd_beh{beh}{rec}.time = linspace(0, obj.timescale.durations(rec), obj.timescale.tp(rec));
                        downsampd_beh{beh}{rec}.value = NaN(1,obj.timescale.tp(rec));
                    end
                if detrend_sig
                    downsampd_beh{beh}{rec}.value = downsampd_beh{beh}{rec}.value - movmin(downsampd_beh{beh}{rec}.value, [20, 0]);
                end

                    %downsampd_beh{beh}{rec}.value = smoothdata(downsampd_beh{beh}{rec}.value, 'movmean', [5, 0]);

                    temp.time = [temp.time, downsampd_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec)];
                    value = downsampd_beh{beh}{rec}.value;
                    if size(value, 1) > 1 % if your behavioural metrics is made of multiple arrays
                        value = nanmax(value, [], 1);
                    end
                    temp.value = [temp.value, value];
                end

                concat_downsamp_beh.time = [concat_downsamp_beh.time ;temp.time];
                concat_downsamp_beh.value = [concat_downsamp_beh.value ;temp.value];

            end

            if rendering
                figure(1026);clf();
                for beh = 1:numel(type)
                    subplot(numel(type),1,beh)
                    for rec = 1:numel(raw_beh{beh})
                        if ~isempty(downsampd_beh{beh}{rec}.time)                           
                            plot(downsampd_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec), downsampd_beh{beh}{rec}.value);hold on;
                        end
                    end
                end
            end

            function out = cellerror_empty(~,varargin)
                out = {};
                out.time = [];
                out.value = [];
            end
        end

        function [bouts, beh_sm, active_tp] = get_activity_bout(obj, beh_types, rendering, smoothing)
            if nargin < 2 || isempty(beh_types)
                beh_types = {'encoder'};
            elseif ischar(beh_types)
                beh_types = {beh_types};
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end
            if nargin < 4 || isempty(smoothing)
                smoothing = 1;
            elseif ~any(smoothing)
                smoothing(1) = 1;
            end
            if numel(smoothing) == 1
                smoothing = [smoothing, 0];
            end
            thr = 10 % threshold in % of max
            margin_duration = 3 %in sec
            detrend_win = 100
            
            if numel(margin_duration) == 1
                margin_duration = [margin_duration, margin_duration];
            end

            plts = {};
            for idx = 1:numel(beh_types)
                current_type = beh_types{idx};
                [~, ~, beh] = obj.get_behaviours(current_type, false);
                beh.value = beh.value - movmin(beh.value, [detrend_win, 0]);
                beh_sm = smoothdata(beh.value, 'gaussian', smoothing);
                %beh_sm = detrend(fillmissing(beh_sm,'nearest'),'linear',cumsum(obj.timescale.tp));
                %thr = prctile(beh_sm(beh_sm > 0), 20);
                current_thr = prctile(beh_sm,(100/thr)) + range(beh_sm)/(100/thr); % 5% of max

                %% Define bouts
                active_tp   = abs(beh_sm) > abs(current_thr);
                [starts, stops] = get_limits(active_tp, beh_sm);

                %% Add some pts before and after each epoch
                dt = nanmedian(diff(obj.t));
                for epoch = starts
                    active_tp(max(1, epoch-round(margin_duration(1)*(1/dt))):epoch) = 1;
                end
                for epoch = stops
                    active_tp(epoch:min(numel(active_tp), epoch+round(margin_duration(2)*(1/dt)))) = 1;
                end
                [starts, stops] = get_limits(active_tp, beh_sm); % update bouts edges now that we extended the range

                if rendering
                    figure(1027);hold on;
                    if idx == 1
                        clf();
                    end
                    plts{idx} = subplot(numel(beh_types),1,idx);hold on;
                    title(strrep(current_type,'_','\_'))
                    plot(obj.t, beh_sm);hold on;
                    for el = 1:numel(starts)
                        x = [starts(el),starts(el),stops(el),stops(el)];
                        y = [0,nanmax(beh_sm),nanmax(beh_sm),0];
                        patch('XData',obj.t(x),'YData',y,'FaceColor','red','EdgeColor','none','FaceAlpha',.1);hold on
                    end
                end

                bouts = sort([starts, stops]);
            end

            if rendering
                linkaxes([plts{:}],'x');
            end

            function [starts, stops] = get_limits(active_tp, beh_sm)
                if size(active_tp, 1) > 1
                    active_tp = nanmean(active_tp, 1) > 0;
                end
                starts      = find(diff(active_tp) == 1);
                stops       = find(diff(active_tp) == -1);
                if any(starts)
                    if starts(1) > stops(1)
                        starts = [1, starts];
                    end
                    if stops(end) < starts(end)
                        stops = [stops, size(beh_sm,2)];
                    end
                elseif all(active_tp)
                    starts = 1; stops = size(beh_sm,2);
                end
            end
        end

        %% #############################

        function [tree, soma_location, tree_values, values] = plot_dist_tree(obj, bin)
            if nargin < 2 || isempty(bin)
                bin = [];
            end
            [tree, soma_location, tree_values, values] = obj.ref.plot_dist_tree(bin);
        end

        %         function [tree, soma_location, tree_values, values] = plot_seg_length_tree(obj)
        %             %% Map dimension weights on the tree
        %             [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(Pvec_tree(obj.ref.simplified_tree{1}), '', '', 'Segment length');
        %         end
        

        function animate_tree(obj, values, timepoints, ROIs, color_range)
            if nargin < 2 || isempty(values)
                values = obj.rescaled_traces;
            end
            if nargin < 3 || isempty(timepoints)
                timepoints = 1:size(obj.rescaled_traces,1);
            end
            if nargin < 4 || isempty(ROIs)
            	ROIs =  obj.ref.indices.complete_ROIs_list(:,1);
            end
            if nargin < 5 || isempty(color_range)
            	color_range = [nanmin(values(:)),nanmax(values(:))];
            end 

            %% Get the ROIs to display
            ROI_idx           = obj.ref.indices.complete_ROIs_list(:,1);

            obj.ref.animate_tree(values, timepoints, ROI_idx, color_range)
        end

        %% #############################

        function prepare_binning(obj, condition)

            %% Clear fields depending on a different scaling
            obj.rescaling_info = {};
            obj.need_update(:)              = true;
            if  nargin < 2
                obj.binned_data.condition   = 'single group';
                obj.binned_data.groups      = {1:obj.n_ROIs};
                obj.binned_data.metrics     = 1;
                obj.binned_data.bin_legend  = {'all ROIs'};
            elseif iscell(condition) && ischar(condition{1})
                obj.binned_data.condition   = condition;

                %% Define current binning rule. See arboreal_scan.get_ROI_groups for more info % CC and peak extractions are based on this binning
                [obj.binned_data.groups, obj.binned_data.metrics, obj.binned_data.bin_legend] = obj.ref.get_ROI_groups(obj.binned_data.condition, obj.demo);
            elseif iscell(condition) && ismatrix(condition{1})
                obj.binned_data.condition   = 'custom';
                obj.binned_data.groups      = condition;
                obj.binned_data.metrics     = 1:numel(condition);
                legends = strcat('group ', num2str(1:numel(condition))');
                obj.binned_data.bin_legend  = cellstr(legends(1:3:end,:))';
            end
            obj.binned_data.readmap         = sort(unique([obj.binned_data.groups{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
            obj.set_median_traces(false);
        end

        function rescale_traces(obj, method, smoothing)
            if nargin >= 2 && ~isempty(method)
                obj.rescaling_method = method;
            end
            if nargin < 3 || isempty(smoothing)
                pk_width                = nanmedian(vertcat(obj.event.peak_width{:}));
                smoothing               = [pk_width*2,0];
            end

            %% Set some flag
            obj.is_rescaled             = false;
            invalid                     = ~ismember(1:obj.n_ROIs, obj.ref.indices.valid_swc_rois') | ismember(1:obj.n_ROIs, obj.bad_ROI_list);

            %% Get traces to rescale and Filter out excluded ROIs so they don't mess up the scaling process
            if strcmp(obj.rescaling_method, 'global')
                traces                   = obj.extracted_traces_conc;
                traces(:,invalid)          = NaN;
            elseif contains(obj.rescaling_method, 'by_trials')
                traces                      = obj.extracted_traces;
                for idx = 1:numel(traces)
                    traces{idx}(:,invalid) = NaN;
                end
            else
                error('rescaling method not valid, use "global" or "by_trials"')
            end

            %% Now rescale
            if strcmp(obj.rescaling_method, 'global')
                [~, obj.rescaling_info.offset, obj.rescaling_info.scaling] = tweak_scaling(traces, unique(vertcat(obj.event.peak_time{:})), smoothing);
                obj.rescaling_info.individual_scaling = repmat({obj.rescaling_info.scaling}, 1, numel(obj.extracted_traces));
                obj.rescaling_info.individual_offset = repmat({obj.rescaling_info.offset}, 1, numel(obj.extracted_traces));
            elseif contains(obj.rescaling_method, 'by_trials')
                t_peak_all              = unique(vertcat(obj.event.peak_time{obj.event.is_global}));
                t_for_baseline          = find(~ismember(1:numel(obj.t),unique([obj.event.t_win_no_overlap{:}])));
                [obj.rescaling_info.scaling, obj.rescaling_info.offset, obj.rescaling_info.individual_scaling, obj.rescaling_info.individual_offset, obj.rescaling_info.scaling_weights, obj.rescaling_info.offset_weights] = scale_every_recordings(traces, obj.demo, t_peak_all, t_for_baseline, smoothing); % qq consider checking and deleting "scale_across_recordings"
            end

            obj.is_rescaled = true;
            obj.set_median_traces(true);
            if obj.rendering
                obj.plot_rescaled_traces();
                obj.plot_rescaling_info();arrangefigures([1,2]);
            end
        end

        function rescaled_traces = get.rescaled_traces(obj)
            if isempty(obj.rescaling_info) || ~obj.is_rescaled
                warning('TRACES HAVE NOT BEEN RESCALED YET - CALL obj.rescale_traces() first.\n');
                obj.is_rescaled = false;
                rescaled_traces = [];
                return
            else
                %% Now rescale each trace with a unique value across recordings (which would be spectific of that region of the tree).
                rescaled_traces                         = obj.extracted_traces_conc;
                transition                              = cumsum([1,obj.timescale.tp(1:end-1)]);
                rescaled_traces([transition,transition+1],  obj.bad_ROI_list) = NaN; % remove bad ROis and transition timepoints
                rescaled_traces = rescaled_traces - diag(prctile(rescaled_traces, obj.rescaling_info.offset,1))'; %remove correct offst percentile for ach trace. much faster than  a loop
                rescaled_traces = rescaled_traces ./ obj.rescaling_info.scaling;
            end
        end

        function [global_median, all_traces_per_bin] = set_median_traces(obj, use_rescaled)
            if (nargin < 2 || isempty(use_rescaled) || use_rescaled) && ~isempty(obj.rescaling_info)
                traces = obj.rescaled_traces;
                obj.is_rescaled = true;
            else
                traces = obj.extracted_traces_conc;
                obj.is_rescaled = false;
            end

            %% Create median trace per bins
            all_traces_per_bin = cell(1, numel(obj.binned_data.groups));
            for gp = 1:numel(obj.binned_data.groups)
                columns                     = ismember(obj.binned_data.readmap, obj.binned_data.groups{gp}) & ~ismember(obj.binned_data.readmap, obj.bad_ROI_list);
                all_traces_per_bin{gp}      = nanmedian(traces(:,columns), 2);
            end

            global_median                   = nanmedian(traces, 2);
            obj.binned_data.global_median   = global_median;
            all_traces_per_bin              = cell2mat(all_traces_per_bin);
            obj.binned_data.median_traces   = all_traces_per_bin;

            if obj.rendering
                obj.plot_median_traces();arrangefigures([1,2]);
            end
        end

        function precision = compute_similarity(obj)
            win                         = ceil(1./median(diff(obj.t)));
            if size(obj.binned_data.median_traces,2) > 1
                comb                        = nchoosek(1:size(obj.binned_data.median_traces,2),2);
                corr_results                = {};
                for pair = 1:size(comb,1)
                    corr_results{pair}      = movcorr(obj.binned_data.median_traces(:,comb(pair, 1)),obj.binned_data.median_traces(:,comb(pair, 2)),[win, 0]);
                end

                corr_results                = cell2mat(corr_results);
                precision                   = 1./nanvar(corr_results,[], 2);
                out                         = nanmax(precision(:)) / 10;
                precision(precision > out)  = NaN;
                obj.variability.corr_results= corr_results;
                obj.variability.precision   = precision;
            else
                obj.variability.corr_results= NaN(size(obj.binned_data.median_traces));
                obj.variability.precision   = NaN(size(obj.binned_data.median_traces));
            end
            if obj.rendering
                obj.plot_similarity();arrangefigures([1,2]);
            end
        end

        function norm_cumsum = get_events_statistics(obj)

            peaks = obj.event.fitting.post_correction_peaks;

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
            cmap = lines(size(obj.binned_data.median_traces, 2));
            [m,n] = numSubplots(size(obj.binned_data.median_traces, 2));
            for gp = 1:size(obj.binned_data.median_traces, 2)
                subplot(m(1), m(2), gp)
                %figure();plot(global_timescale,obj.binned_data.median_traces(:,gp));hold on;scatter(obj.event.fitting.peak_times,obj.event.fitting.post_correction_peaks(:,gp), 'kv')
                if ~isempty(peaks)
                    histogram(peaks(:,gp),0:bin_size:max_peaks, 'FaceAlpha', 0.8,'EdgeColor','none', 'FaceColor', cmap(gp, :));hold on;
                end
                title(obj.binned_data.bin_legend(gp));
            end

            %% Plot all individual peaks to see if they covary
            figure(1006);cla();plot(peaks);legend(obj.binned_data.bin_legend);title('peaks values per subgroups'); xlabel('event #'); ylabel('Amplitude');set(gcf,'Color','w');

            %% Plot the mean amplitude per subgroup of peaks, per distance
            %bin_step = 10;
            %norm_cumsum = cumsum(obj.binned_data.median_traces) ./ nanmax(cumsum(obj.binned_data.median_traces));
            norm_cumsum = cumsum(peaks) ./ nanmax(cumsum(peaks));
            figure(1007);cla();plot(norm_cumsum); hold on;set(gcf,'Color','w');xlabel('Event #');ylabel('normalized cumulative amplitude')
            title('cumulative sum of peaks'); ylim([0, 1]);legend(obj.binned_data.bin_legend,'Location','southeast');
        end

        function norm_vmr = assess_variability(obj)
            %sr = nanmedian(diff(obj.timescale{expe}.global_timescale));
            vmr = nanvar(obj.event.fitting.post_correction_peaks,[],2)./nanmean(obj.event.fitting.post_correction_peaks, 2);
            %cv  = nanstd(obj.event.fitting.post_correction_peaks,[],2)./nanmean(obj.event.fitting.post_correction_peaks, 2); % (maybe chack snr at one point?  mu / sigma)
            %fano = []; % windowed VMR. usually for spike trains
            [~, idx] = sort(vmr,'descend');

            %figure(123); cla();ylim([0, nanmax(obj.event.fitting.post_correction_peaks(:))]); hold on;
            %     for event = idx'
            %         show_event(obj.binned_data.median_traces, round(obj.event.fitting.peak_times/sr), event);
            %         drawnow;%pause(0.1)
            %     end

            figure(1009);cla();plot(obj.event.fitting.post_correction_peaks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')

            R = max(range(obj.binned_data.median_traces));
            %norm_vmr = vmr/range(vmr);
            norm_vmr = vmr/mean(vmr);
            obj.variability.index_of_disp = norm_vmr;
            figure(1010);clf();
            ax1 = subplot(2,1,1);plot(obj.t, obj.binned_data.median_traces); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
            ax2 = subplot(2,1,2);plot(obj.event.fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
            plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
            plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
            hold on;plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');
            linkaxes([ax1, ax2], 'x');

            %figure();histogram(norm_vmr, round(10*range(norm_vmr)/std(norm_vmr)))


            %             figure()
            %             [~, ~, beh] = obj.get_behaviours('encoder');
            %             behaviour = smoothdata(beh.value, 'gaussian', [50, 0]);
            %             hold on;plot(obj.t, behaviour,'b');
            %
            %             plot(obj.event.fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
            %             plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
            %             plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
            %             hold on;plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');
            %
            %
            %             temp = behaviour(round(obj.event.fitting.peak_pos));
            %             hold on;plot(obj.t(obj.event.fitting.peak_pos), temp,'-vr');

            figure(1011);cla();scatter(nanmedian(obj.event.fitting.post_correction_peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
        end

        %% #############################################
        
        function set.cc_mode(obj, cc_mode)
        	if ~strcmp(obj.cc_mode,cc_mode)
                obj.crosscorr = []; %if you change the cc mod, clear the previous correlation results
            end
            obj.cc_mode = cc_mode;
        end
        
        
        function cross_corr = get_correlations(obj, cc_mode)
            if nargin > 1
                obj.cc_mode = cc_mode; % change cc mode
            end

            %% Get CC (and build the correlation matrix if the settings changed)
            cross_corr = obj.crosscorr;
            
            if obj.rendering
                obj.plot_correlation_results(cross_corr);
            end
        end
        
        function set_crosscorr(obj)
            mode = obj.cc_mode;

            %% Get signal time range
            tp_of_events        = sort(unique([obj.event.t_win{:}]));
            if contains(obj.cc_mode, 'peaks') %% event time
                tp          = obj.event.fitting.peak_pos;% could be using obj.event.fitting.pre_correction_peaks
            elseif contains(obj.cc_mode, 'quiet') %% low corr window
                tp              = true(1,size(obj.binned_data.median_traces,1));
                tp(tp_of_events)= false;
            elseif contains(obj.cc_mode, 'active') %% high corr window
                tp              = false(1,size(obj.binned_data.median_traces,1));
                tp(tp_of_events)= true;
            else %% all tp
                tp    = deal(1:size(obj.rescaled_traces,1));
            end
            mode = erase(mode, {'peaks','quiet','active'});

            %% Get ref ROIs
            somatic_ROIs = obj.ref.indices.somatic_ROIs;

            %% Get signal to use
            if contains(obj.cc_mode, 'groups')
                signal = obj.binned_data.median_traces(tp,:);
                ref    = nanmean(obj.rescaled_traces(tp, somatic_ROIs),2);
            else% if contains(obj.cc_mode, 'ROIs')
                signal = obj.rescaled_traces(tp,obj.ref.indices.valid_swc_rois);
                ref    = nanmean(obj.rescaled_traces(tp, somatic_ROIs),2);
            end
            mode = erase(mode, {'ROIs','groups'});

            %% Add population signal if needed
            if contains(obj.cc_mode, 'pop')
                if isempty(obj.extracted_pop_conc)
                    warning('No population data for this rcording')
                    pop = [];
                else
                    pop = obj.extracted_pop_conc(tp,:);
                end
            else
                pop = [];
            end
            mode = erase(mode, {'pop','_'});

            beh = obj.behaviours.types(find(contains(obj.behaviours.types, mode)));
            if ~isempty(beh)
                [~, ~, beh] = obj.get_behaviours(beh);
                ref         = beh.value(tp)';
            end

            variable    = [ref, signal, pop];
            obj.crosscorr   = corrcoef(variable,'Rows','Pairwise')';
        end

        function crosscorr = get.crosscorr(obj)
            %% Add a few more check
            if isempty(obj.binned_data)
                crosscorr = [];
                return
            elseif isempty(obj.crosscorr) || size(obj.crosscorr, 1) ~= numel(obj.binned_data)
                obj.set_crosscorr(); 
            end
            crosscorr = obj.crosscorr;
        end

        function [tree, soma_location, tree_values, mean_bin_cc] = plot_corr_tree(obj, cc)
            if nargin < 2 || isempty(cc)
                cc = obj.crosscorr(2:end,2:end);
                if contains(obj.cc_mode,'pop')
                    pop_sz = size(obj.extracted_pop_conc,2);
                    cc = cc(1:(end-pop_sz),1:(end-pop_sz));
                end
            end

            if size(cc, 1) == numel(obj.binned_data.bin_legend)
                %% Identify valid set of traces
                valid_gp            = find(~all(isnan(obj.binned_data.median_traces))); % You get NaN'ed bins if the soma location is not scanned (eg a big pyramidal cell)

                %% Build tree values per bin
                mean_bin_cc     = [];
                ROIs_list       = [];
                if ~isempty(cc)
                    for gp = 1:numel(obj.binned_data.groups)
                        roi_of_gp           = obj.binned_data.groups{gp};
                        roi_of_gp = roi_of_gp(~ismember(roi_of_gp, obj.bad_ROI_list)); %% COMMENT OUT TO INCLUDE BAD ROIS
                        v_of_gp             = cc(valid_gp(1),gp);
                        ROIs_list           = [ROIs_list, roi_of_gp];
                        mean_bin_cc         = [mean_bin_cc, repmat(v_of_gp, 1, numel(roi_of_gp))];
                    end
                end
            else
                ROIs_list   = obj.ref.indices.valid_swc_rois;
                mean_bin_cc = cc(:,1);
            end

            %% Map CC values on the tree
            if contains(obj.cc_mode, 'group')
                ROIs_list = obj.binned_data.groups;
            end
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(mean_bin_cc, ROIs_list, obj.default_handle, 'Correlation with most proximal segment','',1018);
            caxis([0,1]);
            col = colorbar; col.Label.String = 'Spatial correlation between ROIs/groups with soma';
        end

        %% ###################

        function get_spike_trains(obj)
            %% Use ML spike to get a spike train estimate
            if ~exist('spk_autocalibration.m','file') || ~exist('fn_getfile.m','file')
                error('You need to download the "bricks" and "ml_sikes" toolboxes, and add thn to the path')
            end

            %% Formatting calcium. It seems that signal need to be normalized
            calcium = num2cell(fillmissing(normalize_sig(obj.binned_data.median_traces',[], 'norm_method','signal/max','percentile',10)','constant',0), 1);
            dt      = median(diff(obj.t));
            valid   = find(cellfun(@(x) any(x), calcium));
            %first_valid = find(~cellfun(@(x) all(isnan(x)), calcium),1,'first');

            AMAX        = 10; % wtf
            p = obj.event.fitting.global_fit;f = p(2) / (p(2)+p(4));estimated_tau = p(3)*f+p(5)*(1-f);

            %% Auto-calibration
            pax         = spk_autocalibration('par');
            pax.dt      = dt/2; % trick to avoid some error message.
            pax.maxamp  = AMAX;
            pax.amin    = AMAX/100;
            pax.amax    = AMAX; % pax.taumin = 0.01 ; % pax.taumax = 5
            pax.eventa  = AMAX; % increase if baseline gets too high ; try 10-15?
            pax.eventtau= estimated_tau*2; % ok-ish
            pax.saturation = 0.00195;   % for gcamp6f
            pax.driftparam = 0.002;     % that's enough as baseline is clean
            pax.hill    = 2.05;         % for gcamp6f  - may be 2.27

            %% Autoestimate doesn't work on single runs. We run batches of 1000pts, ignore failed estimates and get a consensus estimate
            batch_size  = 1000;
            tauest      = [];
            aest        = [];
            sigmaest    = [];
            % qq we could use the mean to get a better tau estimate
            for r = 1:(floor(numel(calcium{valid(1)})/batch_size) - 1)
                try
                    [tauest(r), aest(r), sigmaest(r), evt, par] = spk_autocalibration(calcium{valid(1)}((r*batch_size):((r+1)*batch_size),1),pax);
                catch
                    tauest(r) = NaN; aest(r) = NaN;  sigmaest(r) = NaN;
                end
            end

            %% Now get MAP
            par                 = spk_est('par');
            par.dt              = dt/2;
            par.a               = nanmedian(aest);
            par.tau             = nanmedian(tauest);
            par.finetune.sigma  = nanmedian(sigmaest)/10; %sigmeest is always too high
            par.ton             = 0;%0.045; --> 0 improves the detection quality, but I'm not sure why
            par.saturation      = 0.00195;  % from their example code with gcamp6f
            par.hill            = 2.05; % from their example code with gcamp6f
            par.F0              = [-0.01,0.02]; % tight range around 0 (are baseline is flat). works better than fixed value
            par.drift.parameter = 0.002; % gives a bit of slack, but not too much since baseline is already good
            obj.spiketrains     = {};
            obj.spiketrains.settings = par;
            for bin = valid
                [obj.spiketrains.spike_estimate{bin},obj.spiketrains.fit{bin},obj.spiketrains.drift{bin}] = spk_est(calcium{bin},par);
                obj.spiketrains.spike_estimate{bin} = obj.spiketrains.spike_estimate{bin}*2; %bc of the signal upsampling
            end
            obj.plot_inferred_spikes();
        end

        function plot_inferred_spikes(obj)
            calcium = num2cell(normalize_sig(obj.binned_data.median_traces',[], 'norm_method','signal/max','percentile',10)', 1);
            valid   = cellfun(@(x) ~all(isnan(x)), calcium);

            figure(1025);cla();plot(obj.t,calcium{1},'k');hold on;title('Spike inference per bin');xlabel('time(s)')
            colors = NaN(numel(valid),3);colors(valid,:) = lines(sum(valid));
            for bin = find(valid)
                plot(obj.t, obj.spiketrains.fit{bin},'Color',colors(bin, :));hold on
                scatter(obj.spiketrains.spike_estimate{bin}, repmat(bin/10, 1, numel(obj.spiketrains.spike_estimate{bin})),'o','MarkerEdgeColor',colors(bin, :));hold on;
            end
        end

        %% ###################
        function weighted_averages = get_dimensionality(obj, cross_validate, n_factors, timepoints, dim_red_type, weigthed_average)
            if nargin < 2 || isempty(cross_validate)
                cross_validate                  = false;
            end
%             if ~isfield(obj.dimensionality, 'n_factors')
%                 obj.dimensionality.n_factors = 5; % temp fix until we regenerate all recordings
%             end
            if (nargin >= 3 && ~isempty(n_factors)) && (isempty(obj.dimensionality) || n_factors ~= obj.dimensionality.n_factors) % if you change the value
                obj.dimensionality           = {};
                obj.dimensionality.n_factors = n_factors;
            elseif nargin < 3 || isempty(n_factors) && ~isfield(obj.dimensionality, 'n_factors') || (isfield(obj.dimensionality, 'n_factors') && isempty(obj.dimensionality.n_factors))
                obj.dimensionality.n_factors = 5;                
            end
            if nargin < 4 || isempty(timepoints)
                timepoints                  = true(size(obj.timescale.global_timescale));
            end
            if nargin < 5 || isempty(dim_red_type)
                dim_red_type                  = 'pca';
            end
            if nargin < 6 || isempty(weigthed_average)
                weigthed_average              = true;
            end
            
            t_mode = 'peaks'
            if strcmp(t_mode, 'peaks')
            	timepoints     = vertcat(obj.event.peak_time{:});
            elseif strcmp(t_mode, 'beh')
                [timepoints, beh_sm, active_tp] = obj.get_activity_bout('BodyCam_Wheel', true, 5);
            end            
            data                = obj.rescaled_traces(timepoints,:);
            
            %rescaled_traces     = rescaled_traces - obj.binned_data.global_median;
            
            all_ROIs            = 1:size(data, 2);
            normal_n_NaN        = median(sum(isnan(data(:,~all(isnan(data)))))) * 4; % get an indicative number of NaN in a normal traces, and set acceptable thr at 4 times that
            valid_trace_idx     = sum(isnan(data)) <= normal_n_NaN; % eclude traces with too many NaNs (eg. traces that got masked completely)
            data                = fillmissing(data(:, valid_trace_idx),'spline'); % removed funny traces
            
            %% Need to set it now
            obj.dimensionality.valid_trace_idx = valid_trace_idx;

            %% Get single or multiple factor estimate
            if ~cross_validate
                T = [];
                stats = {};
                specVarPM = [];
                
                switch dim_red_type
                    case 'pca'
                        %[LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(rescaled_traces),'NumComponents',obj.dimensionality.n_factors);
                        [LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(data) - obj.binned_data.global_median(vertcat(obj.event.peak_time{:})),'NumComponents',obj.dimensionality.n_factors);

                        %[LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(rescaled_traces) - F(:,1)*LoadingsPM(:,1)','NumComponents',obj.dimensionality.n_factors)
                        T = []
                    case 'nnmf'
                        [F,LoadingsPM, D] = nnmf(double(data),obj.dimensionality.n_factors,'replicates',20,'algorithm','als');
                        LoadingsPM = LoadingsPM';
                        T = [];
                        stats = {};
                        specVarPM = [];
                    case 'factoran'
                        [LoadingsPM, specVarPM, T, stats, F] = factoran(double(data), obj.dimensionality.n_factors,'rotate','varimax'); % varimax, quartimax
                end
                

                %% If no weighted average, Find clusters from latent variables
                if ~weigthed_average
                    cluster_idx = kmeans_opt(LoadingsPM,ceil(sqrt(size(LoadingsPM, 1))),0.8,10);
                    %idx = kmeans(LoadingsPM, obj.dimensionality.n_factors);

                    %% Sort clusters by number of elements
                    [~, gp] = sort(hist(cluster_idx,unique(cluster_idx)), 'descend');
                    idx_sorted = NaN(size(cluster_idx));
                    count = 1;
                    for gp_idx = gp  
                        idx_sorted(cluster_idx == gp_idx) = count;
                        count = count + 1;
                    end
                    obj.dimensionality.cluster_idx = idx_sorted;
                else
                    obj.dimensionality.cluster_idx = [];
                end

                %% Store results
                obj.dimensionality.LoadingsPM         = LoadingsPM;
                obj.dimensionality.specVarPM          = specVarPM;
                obj.dimensionality.T                  = T;                % Rotation matrix
                obj.dimensionality.stats              = stats;            % Factoran stats
                obj.dimensionality.F                  = F;                % components
                obj.dimensionality.all_ROIs           = all_ROIs;         % first occurences
                obj.dimensionality.mask               = timepoints;

                %% Plot weight-tree for each component
                for comp = 1:obj.dimensionality.n_factors
                    obj.plot_dim_tree(comp);
                end

                %% Plot a map of tree weights by ROI number for each component
                weighted_averages = obj.get_weight_map();
                if obj.rendering
                    obj.plot_factor_tree();
                end

                %% Plot clusters
                obj.plot_cluster_tree();


                    
%                     %% To assign to strongest component
%                     for row = 1:size(LoadingsPM,1)
%                         [~, maxloc] = max(LoadingsPM(row, :));
%                         LoadingsPM(row, :) = 0;
%                         LoadingsPM(row, maxloc) = 1;                        
%                     end

                    %% To assign to strongest component
%                     for row = 1:size(LoadingsPM,1)
%                         LoadingsPM(row, :) = 0;
%                         LoadingsPM(row, idx(row)) = 1;    
%                         if ~ismember(idx(row), [1,4,5])
%                             LoadingsPM(row, :) = 0;
%                         end                            
%                     end
                
            else
                %% Cross validation (thanks Harsha)
                %[~,N]                       = size(rescaled_traces);
                nFactors                    = 20;%round(0.66*N);
                n_iter                      = 5;
                train                       = NaN(nFactors, n_iter);
                test                        = NaN(nFactors, n_iter);
                fac_steps                   = 1;
                weighted_averages = [];
                tic
                for jj = 1:n_iter
                    %% Partition data into Xtrain and Xtest
                    %% random indexing
                    blocks = false
                    if ~blocks
                        test_idx = randperm(size(data, 1));
                        Xtrain = double(data(test_idx(1:2:end), :));
                        Xtest  = double(data(test_idx(2:2:end), :));
                    else
                        test_idx = randperm(size(data, 1)/10)*10;
                        test_idx(test_idx == max(test_idx)) = [];
                        idx_train = cell2mat(arrayfun(@(x) x:x+9, test_idx(1:2:end-1), 'UniformOutput', false));
                        idx_test = cell2mat(arrayfun(@(x) x:x+9, test_idx(2:2:end-1), 'UniformOutput', false));
                        Xtrain = double(data(idx_train, :));
                        Xtest  = double(data(idx_test, :));
                    end

                    jj
                    parfor factor_nb = 1:nFactors
                        if ~rem(factor_nb-1, fac_steps) % every fac_steps steps, starting at 1
                            try
                            [LoadingsPM, specVarPM, ~, stats, F] = factoran(Xtrain, factor_nb,'maxit',1500);
                            train(factor_nb, jj)                 = stats.loglike;
                            test(factor_nb, jj)                  = testloglike_factorAnalysis(Xtest, nanmean(Xtrain,1), LoadingsPM, specVarPM);
                            catch
                            test(factor_nb, jj)                  = NaN;
                            train(factor_nb, jj)                 = NaN;
                            
                            end
                        end
                    end
                end

                tested_modes = 1:fac_steps:nFactors;
                figure(1028);cla();plot(tested_modes,train(tested_modes,1:n_iter),'o');hold on;plot(tested_modes,nanmean(train(tested_modes,1:n_iter), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('train');
                figure(1029);cla();plot(tested_modes,test(tested_modes,1:n_iter),'o');hold on;plot(tested_modes,nanmean(test(tested_modes,1:n_iter), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('test');

                [~, n_factor] = max(nanmean(test, 2));
                obj.get_dimensionality(false, n_factor, timepoints, dim_red_type, weigthed_average)
            end
        end

        function [tree, soma_location, tree_values, values] = plot_dim_tree(obj, comp, fig_handle)
            if nargin < 2 || isempty(comp) || ~comp
                [tree, soma_location, tree_values, values] = obj.plot_strongest_comp_tree();
                return
            end
            if nargin < 3 || isempty(fig_handle)
                fig_handle = 10200 + comp; % fig number or fig hande
            end


            % check 58, % noise issue 64 'D:/Curated Data/2019-09-24/experiment_1/18-13-20/'

            %% Reload loadings
            LoadingsPM = obj.dimensionality.LoadingsPM;
            Valid_ROIs = obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx);

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
            if obj.rendering || ishandle(fig_handle)
                [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, titl, '',  fig_handle, 'regular');
            end
        end

        function weighted_averages = get_weight_map(obj, weigths_to_show)
            if nargin < 2 || isempty(weigths_to_show) % number or list of factor to display
                weigths_to_show = 1:obj.dimensionality.n_factors;
            end

            %% Recover traces
            rescaled_traces = obj.rescaled_traces;

            %% Reload loadings
            LoadingsPM = obj.dimensionality.LoadingsPM;
            Valid_ROIs = obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx);
            rescaled_traces = rescaled_traces(:, Valid_ROIs);
            [~, loc] = max(LoadingsPM(:,1:weigths_to_show)');

            all_weights         = {};
            weighted_averages   = [];
            for w = weigths_to_show
                L                       = LoadingsPM(:,w);
                %L(L<0.2) = 0;
                all_weights{w}          = L/sum(LoadingsPM(:,w));
                weighted_averages(w, :) = nanmean(rescaled_traces'.* all_weights{w}, 1);
            end
        end

        function [tree, soma_location, tree_values, values] = plot_strongest_comp_tree(obj, n_dim, fig_handle)
            if nargin < 2 || isempty(n_dim)
                n_dim = size(obj.dimensionality.LoadingsPM, 2);
            end
            if nargin < 3 || isempty(fig_handle)
                fig_handle = 10200;
            end
            

            [~, loc]    = nanmax(obj.dimensionality.LoadingsPM(:,1:n_dim),[],2);
            Valid_ROIs  = find(obj.dimensionality.valid_trace_idx);
            values      = NaN(size(Valid_ROIs));
            for ROI = 1:numel(Valid_ROIs)
                values(ROI)     = loc(ROI);
            end
            values = values(~isnan(values));

            %% Map dimension weights on the tree
            if obj.rendering
                [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, 'Location of strongest component','',fig_handle, 'regular', 'jet');
                colorbar('Ticks',1:nanmax(values));colormap(jet(nanmax(values)))
            end
        end

        %% #############################################

        function update_all_signals(obj, new_method)
            if nargin < 2 || isempty(new_method)
                new_method = obj.extraction_method;
            elseif ~any(strcmp(new_method, {'max', 'mean', 'median','min'}))
                error('Only max, min, mean and median method are supported')
            else
                obj.extraction_method = new_method;
            end

            %% Update raw signals changing the compression procedure
            for rec = 1:numel(obj.arboreal_scans)
                obj.arboreal_scans{rec}.extraction_method = new_method;
                if isempty(obj.arboreal_scans{rec}.simple_data)
                    temp = load(obj.extracted_data_paths{rec});
                    obj.arboreal_scans{rec}.simple_data = temp.obj.simple_data;
                    clear temp;
                end

                obj.arboreal_scans{rec}.update_segment_signals(new_method);
                obj.arboreal_scans{rec}.simple_data = [];
                obj.need_update(rec)                    = true;
            end
            error_box('COMPRESSION MODES WERE UPDATED BUT YOU NEED TO SAVE THE ARBOREAL SCANS TO KEEP THIS CHANGE FOR NEXT RELOADING')
        end

        function detect_events(obj, idx_filter, corr_window, thr_for_detection)
            if nargin < 2 || isempty(idx_filter)
                idx_filter = 1:obj.n_ROIs;
            elseif ischar(idx_filter) && strcmp(idx_filter, 'soma')
                idx_filter = obj.ref.indices.somatic_ROIs;
            end
            if nargin < 3 || isempty(corr_window)
                corr_window = [];
            end
            if nargin < 4 || isempty(thr_for_detection)
                thr_for_detection = 0.2;
            end
            
            fprintf(['Now detecting events with global pairwise correlation at least > ',num2str(thr_for_detection),' %% \n'])     

            %% Get original traces (unscaled)
            raw_traces              = obj.extracted_traces_conc(:, idx_filter);

            %% Get original traces
            [obj.event, correlation_res] = detect_events(raw_traces, obj.t, 'corr', 0.2, corr_window, obj.rendering);
            %[obj.event]             = detect_events(raw_traces, obj.t, 'global_amp', 60, corr_window);

            %% Identify and log poorly correlated ROIs
            obj.find_bad_ROIs(correlation_res.corr_results, correlation_res.comb, corr_window, idx_filter);
        end

        function [bad_ROIs, mean_corr_with_others] = find_bad_ROIs(obj, corr_results, comb, corr_window, ROIs)
            if nargin < 4
                corr_window = obj.get_ideal_corr_window();
            end
            if nargin < 5
                ROIs = 1:obj.n_ROIs;
            end
            if nargin < 3
                [corr_results, comb] = get_pairwise_correlations(ROIs, corr_window); % same as in detect_events
            end

            THR_FOR_CONNECTION      = 0.2; % Defines what level of minimal pairwise correlation means "these two ROIs are connected"
            fprintf(['* Now detecting ROIs that are either,\n'...
                     '      - so poorly correlated to the rest of the tree that they probably belong to another cell (or have no signal).\n',...
                     '      - are member of batch_params.excluded_branches \n',...
                     '* Threshold for exclusion is  ',num2str(THR_FOR_CONNECTION),' %% \n'])  
            THR_FOR_CONNECTION
            
            %% Show mean correlation with each ROI
            mean_corr_with_others   = [];
            max_corr_for_cell       = max(obj.event.globality_index(2:end)); % QQ 1st point sometimes show some artifacts
            for key = 1:numel(ROIs)
                corr_results_sub    = corr_results(obj.event.t_corr(obj.event.is_global), comb(:,2) == key | comb(:,1) == key);
                corr_results_sub    = [corr_results_sub(:,1:(key-1)), NaN(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];
                mean_corr           = nanmean(corr_results_sub,1);
                mean_corr_with_others(key) = sum(mean_corr > THR_FOR_CONNECTION);
            end

            %% Normalize to 100%
            mean_corr_with_others = mean_corr_with_others / numel(mean_corr_with_others); % renormalize to max possible corr for this cell
            
            figure(88888);cla();hist(100*mean_corr_with_others,0:2:100); hold on; 
            xlabel('%% of correlation with all other ROIs (from the tree)'); ylabel('counts')
            
            %% Renormalize to max possible corr for this cell
            mean_corr_with_others_norm = mean_corr_with_others / max(mean_corr_with_others); 
            figure(88889);cla();hist(100*mean_corr_with_others_norm,0:2:100); hold on; 
            xlabel('%% of correlation with all other ROIs (from the tree)'); ylabel('counts')
            
            %% Update class variables
            bad_ROIs = mean_corr_with_others_norm < obj.bad_ROI_thr;
            obj.bad_ROI_list = bad_ROIs; % below threshold% of max correlation
            obj.bad_ROI_thr = THR_FOR_CONNECTION;
            bad_ROIs = find(bad_ROIs);
            
            %% ROIs that were manually excluded
            excl = ismember(obj.ref.indices.swc_list(:,4), obj.batch_params.excluded_branches);
            
            RECOVERY_THR = 1-THR_FOR_CONNECTION
            excl_but_not_bad    = excl & ~obj.bad_ROI_list';
            excl_but_good       = excl & (mean_corr_with_others_norm > RECOVERY_THR)';
            if any(excl_but_good)
                 fprintf(['!!! ROIs ',num2str(find(excl_but_good')),' was/were excluded but seem highly correlated\n'])
            end
                            
            if obj.rendering     
                bad = obj.extracted_traces_conc(:, find(obj.bad_ROI_list));
               
                              
                ref = nanmedian(obj.extracted_traces_conc,2);
                ref = ref - prctile(ref, 1);
                
                figure(1031);clf();title(['Bad ROIs (NEVER above ',num2str(obj.bad_ROI_thr*100),' % correlation with the rest of the tree)']);hold on;

                %% Plot normalized excluded traces
               %plot_many_traces(smoothdata(normalize(bad),'gaussian',[20,0]),'','k')
                
                %% Plot normalized excluded traces
              %  plot_many_traces(bad, 1031, 'Color', [0.8,0.8,0.8]);
                
                hold on; plot(ref/nanmax(ref),'r')
                
                regroup_traces(bad', 80, 'pca') 
                
                bad = bad - prctile(bad, 1);
                %plot(smoothdata(bad./nanmax(bad),'gaussian',obj.filter_win),'r');hold on;
%                 plot(bad./nanmax(bad),'r');hold on;
%                 plot(ref/nanmax(ref),'k');hold on;
                missed = obj.extracted_traces_conc(:, find(excl_but_not_bad));
                missed = missed - prctile(missed, 1);
                figure(1031);hold on;
                plot(missed/nanmax(missed),'b');hold on;
                invalid = zeros(1,size(obj.ref.indices.swc_list,1));
                invalid(obj.bad_ROI_list) = 1;
                f = obj.ref.plot_value_tree(invalid, 1:numel(invalid), obj.default_handle, 'Uncorrelated ROIs', '',  1032,'','RedBlue'); hold on;
                if any(excl_but_not_bad)
                    obj.ref.plot_value_tree(repmat(0.7,1,sum(excl_but_not_bad)), find(excl_but_not_bad), obj.default_handle, 'Uncorrelated ROIs', '',  f.Parent,'','RedBlue'); hold on;
                end
%                 recovered = ((excl | obj.bad_ROI_list') & ~excl_but_good);
%                 if any(recovered)
%                     obj.ref.plot_value_tree(repmat(0.2,1,sum(recovered)), find(recovered), obj.default_handle, 'Uncorrelated ROIs', '',  f.Parent,'','RedBlue'); hold on;
%                 end
                caxis([0,1]); % otherwise if all values are the same you get a white tree on a white bkg
            end
            %obj.bad_ROI_list = find((excl | obj.bad_ROI_list') & ~excl_but_good);
            obj.bad_ROI_list = find((excl | obj.bad_ROI_list'));  
            %obj.bad_ROI_list = find(excl); 
            %obj.bad_ROI_list = find(obj.bad_ROI_list); 
        end

        function [corr_results, comb] = get_pairwise_correlations(idx_filter, corr_window)
            if nargin < 2 || isempty(idx_filter)
                idx_filter = 1:obj.n_ROIs;
            elseif ischar(idx_filter) && strcmp(idx_filter, 'soma')
                idx_filter = obj.ref.indices.somatic_ROIs;
            end
            if nargin < 3 || isempty(corr_window)
                med = nanmedian(obj.binned_data.median_traces,2);
                med = med(~isnan(med));
                bsl_guess = rms(med)*2;
                [~,~,w] = findpeaks(nanmedian(obj.binned_data.median_traces,2),'SortStr','descend','MinPeakProminence',bsl_guess);
                corr_window = ceil(nanmean(w)); % value set as an asymetrical filter in generate_pairwise_correlations ([corr_window, 0])
            end
            [corr_results, comb]    = generate_pairwise_correlations(obj.extracted_traces_conc(:, idx_filter), obj.event.corr_window); % same as in detect_events
        end


        function process(obj, condition, filter_win, rendering)
            %% Call all processing steps
            if nargin < 2 || isempty(condition)
                condition = {'distance',Inf};
            end
            if nargin < 3 || isempty(filter_win)
                obj.filter_win = [0,0];
            else
                obj.filter_win  = filter_win;
            end
            if nargin >= 4 && ~isempty(rendering)
                obj.rendering = rendering;
            end

            %% Remove previous plots
            obj.clear_plots();

            %% Load and concatenate traces for the selected experiment
            % obj.load_extracted_data();   % this also sets the current expe #
            % quality_control(1, obj.extracted_traces)
            % quality_control(2, cell2mat(all_sr))

            %% Prepare binning of ROIs based on specific grouping condition
            obj.prepare_binning(condition);

            %% Find peaks based on amplitude AND correlation
            obj.detect_events();

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces(); % note that signal rescaling is computed on peaks, not noise

            %% Create median trace per bins
            obj.set_median_traces()

            %% Get summary covariance plots of the raw traces
            obj.compute_similarity();

            %% Correct for decay to avoid overestimating peak amplitude
            obj.fit_events();

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.get_events_statistics();

            %% Study bAP heterogeneity
            obj.assess_variability()

            %% Check how correlated are the peaks between different parts of the tree
            obj.get_correlations();

            %% Dimensionality reduction
            obj.get_dimensionality(false,'','','pca'); % set to true for cross validation

            %% Optionally, if external variables need an update, do it here
            %obj.update_external_metrics(60)
            
            %% Extract behaviours
            obj.get_behaviours();

            %% Spike inference
            try
            %    obj.get_spike_trains();
            end

            %% Save figure and/or analysis
            obj.auto_save_figures = true
            obj.auto_save_analysis = true
            if obj.auto_save_figures
                obj.save_figures();
            end
            if obj.auto_save_analysis
            	obj.save(true);
            end
        end

        function save(obj, auto)
            %% Save all data in a mat file
            if nargin < 2 || isempty(auto) || ~auto
                save_folder = uigetdir(obj.source_folder);
            else
                save_folder = obj.source_folder;
            end
            if any(save_folder) || obj.auto_save_analysis
                name = strsplit(obj.source_folder, '/');
                name = name{end-1};
                save([obj.source_folder, name],'obj','-v7.3')
            end
        end

        function save_figures(obj, file_type)
            if nargin < 2 || isempty(file_type)
                file_type = {'pdf'};
            elseif ischar(file_type)
                file_type = {file_type};
            end
            %% See arboreal_scan_plotting for an index of the figures
            p = get(groot,'DefaultFigurePosition');
            folder = parse_paths([obj.source_folder, '/figures/']);
            if isfolder(folder)
                rmdir(folder,'s');
            end
            mkdir(folder);
            [~, tag] = fileparts(fileparts(obj.source_folder));
            list = obj.get_fig_list();
            for fig_idx = list
                f = figure(fig_idx);
                set(f, 'Position', p)
                try
                    if any(contains(file_type, 'pdf'))
                        saveas(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.pdf']);
                    end
                    if any(contains(file_type, 'png'))
                        saveas(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.png']);
                    end
                    if any(contains(file_type, 'fig'))
                        dcm_obj = datacursormode(f);
                        bkp_callback = get(dcm_obj,'UpdateFcn');
                        set(dcm_obj,'UpdateFcn',[], 'enable', 'off');
                        savefig(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.fig']);
                        if ~isempty(bkp_callback)
                            set(dcm_obj,'UpdateFcn',bkp_callback, 'enable', 'on');
                        end
                    end
                catch % for figures with subplots
                    try
                        if any(contains(file_type, 'pdf'))
                            saveas(f, [folder,'/',f.Tag,' ', tag,'.pdf']);
                        end
                        if any(contains(file_type, 'png'))
                            saveas(f, [folder, '/',tag,'/',f.Tag,' ', tag,'.png']);
                        end
                        if any(contains(file_type, 'fig'))
                            savefig(f, [folder, '/',tag,'/',f.Tag,' ', tag,'.fig']);
                        end
                    end
                end
            end
        end
    end
end
