classdef signal_manipulation < handle
    %% subclass of arboreal_scan_experiment 
    properties
        %% Analysis/extraction settings
        time_smoothing          = [0, 0];       % Smoothing kernet to apply to the extracted Ca2+ signals. 2x1 Int indicate asymatrical kernel. Use negative value for smoothing windo in seconds
        filter_type             = 'gaussian';   % Type of time filter kernel (one of the valid method is smoothdata())
        bad_ROI_thr             = 0.2;          % Threshold to delimit bad vs good ROIs. Average cross correlation across ROIs < bad_ROI_thr are excluded from the analysis
        detrend                 = false;
        is_detrended            = false;        % is set to True once you ran the detrending once.
        is_rescaled             = false;        % is set to True once you ran the rescaling step.
        rescaling_method        = 'by_trials_on_peaks'; 
        breakpoints             = []; % if you had a disruptive event during th experiment, a first scaling is done with large blocks
        use_hd_data             = false;
        bad_ROI_list            = 'unset';      % list of uncorrelated ROIs (following event detection)
    end
    
    properties (Dependent = true, Transient = true)
        rescaled_traces         % Rescaled Traces according to rescaling_info
    end

    methods
        function breakpoints = get.breakpoints(obj)
            %% Get breakpoints from the batch_params field
            % -------------------------------------------------------------
            % Syntax:
            %   breakpoints = obj.breakpoints
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   breakpoints (1xN CELL ARRAY OF CHAR)
            %   List of breakpoints if any
            % -------------------------------------------------------------
            % Extra Notes:
            %   * Breakpoints can be manually provided to indicate abrupt
            %   changes in the experiment signal, which are usually caused
            %   by change in data acquisition variables, such as PMT gain,
            %   water level, or if you interrupted and restarted the
            %   experiment.
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            if isfield(obj.batch_params, 'breakpoints') && ~isempty(obj.batch_params.breakpoints)
                breakpoints = obj.batch_params.breakpoints;
            else
                breakpoints = [];
            end
        end

        function set.breakpoints(obj, breakpoints)
            %% Update breakpoints if not initially provided
            % -------------------------------------------------------------
            % Syntax:
            %   obj.breakpoints = breakpoints;
            % -------------------------------------------------------------
            % Inputs:
            %   breakpoints (1xN INT OR 1xN CELL ARRAY OF CHAR)
            %       Numerical list of experiment number that interrupted the
            %   experiment, or expe tag. for example :
            %       obj.breakpoints = [2,6,9];
            %       obj.breakpoints = {'12-59-59,'13-05-01'};
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %   * Breakpoints can be manually provided to indicate abrupt
            %   changes in the experiment signal, which are usually caused
            %   by change in data acquisition variables, such as PMT gain,
            %   water level, or if you interrupted and restarted the
            %   experiment.
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            obj.is_detrended = false; %if breakpoints are changed, the detrended needs a refresh too
            if isnumeric(breakpoints) && ~isempty(breakpoints)
                breakpoints = sort(breakpoints);
                breakpoints = obj.extracted_data_paths(breakpoints);
                breakpoints = cellfun(@(x) strsplit(fileparts(x),'/'), breakpoints, 'UniformOutput', false);
                breakpoints = cellfun(@(x) x{end}(end-7:end), breakpoints, 'UniformOutput', false);
                %find(cellfun(@(x) contains(x, obj.batch_params.breakpoints),obj.updated_data_path));
            end
            for rec = 1:numel(obj.arboreal_scans)
                obj.arboreal_scans{rec}.batch_params.breakpoints = sort(breakpoints);
            end
        end

        function set.bad_ROI_thr(obj, value)
            %% Defines the exclusion threshold based on correlation
            % This defines the lowest acceptable correlation coefficient
            % between one ROI and the rest of the tree.
            % -------------------------------------------------------------
            % Syntax:
            %   obj.bad_ROI_thr = bad_ROI_thr;
            % -------------------------------------------------------------
            % Inputs:
            %   bad_ROI_thr (0 > FLOAT > 1)
            %   Minimal correlation coefficient between one ROI and the
            %   rest of the tree to be considered part of the tree
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            %
            % See also : find_bad_ROIs

            if value < 0 || value > 1
                obj.disp_info('Cutoff must be between 0 and 1',4)
            elseif isequal(obj.bad_ROI_thr,value)
                % pass
            else
                obj.bad_ROI_thr = value;
                obj.disp_info('Changing bad_ROI_thr may discard/included new ROIs. It is recommended to call obj.reset() and reprocess the data',3)
                a = dbstack();
%                 if ~contains([a.name], 'arboreal_scan_experiment.reset') % not useful when resetting or initializing
%                     obj.find_bad_ROIs();
%                 end
            end
        end
        
        function set.time_smoothing(obj, time_smoothing)
            %% Defines the gaussian filtering window to use across analyses
            % -------------------------------------------------------------
            % Syntax:
            %   obj.time_smoothing = time_smoothing;
            % -------------------------------------------------------------
            % Inputs:
            %   time_smoothing (FLOAT OR 2x1 FLOAT)
            %   symetrical or asymetrical gaussian filter applied to all
            %   traces. negative values indicates that the value is in
            %   second, and conversion into timepoints is done
            %   automatically
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   19/10/2022
            %
            % See also :             
            
            if all(isnumeric(time_smoothing)) && numel(time_smoothing) == 2 && all(time_smoothing == obj.time_smoothing)
                obj.disp_info('No change detected in time smoothing',1)
                %% no change, pass                
            elseif all(isnumeric(time_smoothing)) && numel(time_smoothing) == 1 || numel(time_smoothing) == 2 
                try
                    time_smoothing(time_smoothing < 0)  = time_smoothing(time_smoothing < 0) * nanmedian(1./obj.timescale.sr);
                    time_smoothing                  = abs(round(time_smoothing));
                    if numel(time_smoothing) == 2 && all(time_smoothing == obj.time_smoothing)
                        %% no change, but the value was initially in seconds so we only see it here. Then do nothing
                        obj.disp_info('No change detected in time smoothing',1)
                        return
                    end
                end
                time_smoothing                  = abs(round(time_smoothing));                
                if numel(time_smoothing)    == 1                 
                    time_smoothing = [time_smoothing, time_smoothing];
                end
                obj.time_smoothing              = time_smoothing;               
                obj.reset();
            else
                obj.disp_info('filter window must be a set of one (for symmetrical gaussian kernel) or 2 (for asymetrical gaussian kernel) values. Values are rounded. If values are < 1, window is converted in seconds',4)
            end
        end
        
        function bad_ROI_list = get.bad_ROI_list(obj)
            bad_ROI_list = obj.bad_ROI_list;
            if ischar(bad_ROI_list) && strcmpi(bad_ROI_list, 'unset')
                obj.disp_info('bad_ROI_list has never been set. Detecting now. If events were not detected before, this step may take a few seconds. To ignore this step "obj.bad_ROI_list = []"');
                %rendering = obj.rendering;obj.rendering = false;
                obj.find_events(); %find uncorrelated ROIs based on correlation
                %obj.rendering = rendering;
                bad_ROI_list = obj.bad_ROI_list;
            end

            %% If analyzing every pixel, convert bad ROIs to bad pixels
            if obj.use_hd_data
                bad_ROI_list = obj.get_voxel_for_ROI(bad_ROI_list');
            end
        end
        
        function set.detrend(obj, value)
            if value ~= obj.detrend
                obj.is_detrended = false;
                obj.detrend = value;
                obj.rescaled_traces; % force detrending update
            end
        end


        function set.use_hd_data(obj, use_hd_data)
            %% Set use_hd_data variable. This changes the data used for computations
            % -------------------------------------------------------------
            % Syntax:
            %   external_variables = obj.external_variables;
            % -------------------------------------------------------------
            % Inputs:
            %   use_hd_data (BOOL)
            %       if true, and if obj.ref.full_data is present, set
            %       obj.use_hd_data to true.
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %   * if obj.use_hd_data is true, all voxels are used for
            %   computation instead of one value per ROI. This is possible
            %   only if the arboreal scan_experiment contains the full_data
            %   (which is not the default extraction behaviour, as it makes
            %   the files much larger). if you intend to use it, you need
            %   to build your objects using the keep_2D flag :
            %   arboreal_scan_experiment('',true)
            %   * Using use_hd_data makes the computations significantly
            %   slower.
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   17/06/2022
            
            if use_hd_data && ~isempty(obj.ref.full_data)
                obj.use_hd_data = true;
            elseif use_hd_data && isempty(obj.ref.full_data)
                obj.use_hd_data = false;
                obj.disp_info({'unuable to use HD data as it is not embedded in the current arboreal_scan_experiment object.','Rebuild the object with HD data using expe = arboreal_scan_experiment([arboreal_scans folder PATH],true).','obj.use_hd_data was set to false.'},3);
            else
                obj.use_hd_data = false;
            end
        end
        
        function rescale_traces(obj, method, smoothing)
            if nargin >= 2 && ~isempty(method)
                obj.rescaling_method = method;
            end
            smoothing = 0;
            
            if isempty(obj.event)
                obj.disp_info('LD RESCALING REQUIRES DETECTED EVENTS. You must run obj.find_events() first',3);
                answ = questdlg({'RESCALING REQUIRES DETECTED EVENTS','To rescale traces, you must run obj.find_events() first; Run peak detection now?'},'','Yes','No','Yes');
                if strcmp(answ, 'Yes')
                    obj.find_events();
                else
                    return
                end
            end
            
%             if (nargin < 3 || isempty(smoothing)) && ~obj.use_hd_data
%                 if isempty(obj.event)
%                     warning('LD RESCALING REQUIRES DETECTED EVENTS. You must run obj.find_events()')   
%                     return
%                 end
%                 pk_width                = nanmedian(vertcat(obj.event.peak_width{:}));
%                 smoothing               = [pk_width*2,0];
%             elseif (nargin < 3 || isempty(smoothing)) && obj.use_hd_data
%                 %pass
%             elseif numel(smoothing) == 1
%                 smoothing = [smoothing, 0];
%             end

            %% Prepare the trace to use (whoe trace cocnatenated, or individidual trials)
            invalid                     = ~ismember(1:obj.n_ROIs, obj.ref.indices.valid_swc_rois') | ismember(1:obj.n_ROIs, obj.bad_ROI_list);

            %% Get traces to rescale and Filter out excluded ROIs so they don't mess up the scaling process
            if contains(obj.rescaling_method, 'global')
                traces                   = obj.extracted_traces_conc;
                traces(:,invalid)        = NaN;
            elseif contains(obj.rescaling_method, 'trials')
                traces                    = obj.extracted_traces;
                for idx = 1:numel(traces)
                    traces{idx}(:,invalid) = NaN;
                end
            else
                obj.disp_info('rescaling method not valid, It must contain "global" or "trials"',4);
            end
            
            %% Set some flag
            obj.is_rescaled             = false;
%             if ~obj.use_hd_data && contains(obj.rescaling_method, 'peaks')
                %% Now rescale
                if contains(obj.rescaling_method, 'global')
                    [~, obj.rescaling_info.offset, obj.rescaling_info.scaling] = tweak_scaling(traces, unique(vertcat(obj.event.peak_time{:})), smoothing);
                    obj.rescaling_info.individual_scaling = repmat({obj.rescaling_info.scaling}, 1, numel(obj.extracted_traces));
                    obj.rescaling_info.individual_offset = repmat({obj.rescaling_info.offset}, 1, numel(obj.extracted_traces));
                elseif contains(obj.rescaling_method, 'trials')
                    t_peak_all              = unique(vertcat(obj.event.peak_time{obj.event.is_global}));
                    t_for_baseline          = find(~ismember(1:obj.tp,unique([obj.event.t_win_no_overlap{:}])));
                    [obj.rescaling_info.scaling, obj.rescaling_info.offset, obj.rescaling_info.individual_scaling, obj.rescaling_info.individual_offset, obj.rescaling_info.scaling_weights, obj.rescaling_info.offset_weights] = scale_every_recordings(traces, obj.demo, t_peak_all, t_for_baseline, smoothing); % qq consider checking and deleting "scale_across_recordings"
                end                
%             else
%                 warning('to finish')
%                 obj.bad_ROI_list         = [];
%                 bsl                      = mode(traces);
%                 temp                     = sort(traces, 1);
%                 for idx = 1:size(traces, 2)
%                     obj.rescaling_info.offset(idx) = NaN;
%                     if ~all(isnan(temp(:,idx)))
%                         obj.rescaling_info.offset(idx) = 100* find(temp(:,idx) > bsl(idx), 1, 'first') / size(traces, 1);
%                         traces(:, idx)                   = traces(:, idx)  - prctile(traces(:, idx) , obj.rescaling_info.offset(idx), 1);
%                     end
%                 end
%                 ref                      = nanmedian(traces, 2);
%                 ref                      = ref - prctile(ref, 5);
%                 for idx = 1:size(traces, 2)
%                     valid = ~isnan(traces(:,idx));
%                     obj.rescaling_info.scaling(idx) = 1/(traces(valid,idx) \ ref(valid));
%                 end
%                 obj.rescaling_info.scaling = obj.rescaling_info.scaling';
%                 obj.rescaling_info.individual_scaling = repmat({obj.rescaling_info.scaling}, 1, numel(obj.extracted_traces));
%                 obj.rescaling_info.individual_offset = repmat({obj.rescaling_info.offset}, 1, numel(obj.extracted_traces));
%             end

            obj.is_rescaled = true;
            obj.set_median_traces(true);
            if obj.rendering
                obj.plot_rescaled_traces();
                obj.plot_rescaling_info();arrangefigures([1,2]);
            end
        end

        function rescaled_traces = get.rescaled_traces(obj)
            %% Rescale traces using precomputed rescaling info
            % NaN the ROIs in bad_ROI_list 
            % Nan the points between two recordings to prevent some stitching artefact or help you find these locations
            if isempty(obj.rescaling_info) || ~obj.is_rescaled
                obj.disp_info('TRACES HAVE NOT BEEN RESCALED YET - CALL obj.rescale_traces() first',3);
                obj.is_rescaled = false;
                rescaled_traces = [];
                return
            else
                %% Now rescale each trace with a unique value across recordings (which would be spectific of that region of the tree).
                rescaled_traces                         = obj.extracted_traces_conc;
                transition                              = cumsum([1,obj.timescale.tp(1:end-1)]);
                rescaled_traces([transition,transition+1],  obj.bad_ROI_list) = NaN; % remove bad ROis and transition timepoints
                rescaled_traces = rescaled_traces - diag(prctile(rescaled_traces, obj.rescaling_info.offset,1))'; %remove correct offset percentile for each trace. much faster than  a loop
                try
                    obj.rescaling_info.scaling; % debug hack. somtimes it needs to be called twice at initialisation --> to be fixed
                end
                scaling = obj.rescaling_info.scaling;
                scaling(isinf(scaling)) = NaN;
                try
                    rescaled_traces = rescaled_traces ./ scaling;
                catch
                    rescaled_traces = rescaled_traces ./ scaling'; %qq to fix --> hapens with hd data only
                end
                
                %rescaled_traces = rescaled_traces - nanmedian(rescaled_traces, 2);                
            end
        end
        
        function [bad_ROIs, mean_corr_with_others] = find_bad_ROIs(obj, correlation_res, ROIs)
            if nargin < 3 || isempty(ROIs)
                ROIs = 1:obj.n_ROIs;
            end      
            mean_corr_with_others   = [];
            bad_ROIs                = [];
            if nargin < 2 || isempty(correlation_res)
                obj.disp_info('BAD ROI IDENTIFICATION RELIES ON EVENT DETECTION. EVENT DETECTION FIELD SEEMS EMPTY.  run obj.find_events() first',3);
                return
            else
                events = obj.event;
            end
            
            
            THR_FOR_CONNECTION      = obj.bad_ROI_thr; % Defines what level of minimal pairwise correlation means "these two ROIs are connected"
            obj.disp_info({'Now detecting ROIs that are either,'...
                            '- so poorly correlated to the rest of the tree that they probably belong to another cell (or have no signal).',...
                            '- are member of batch_params.excluded_branches ',...
                            ['Threshold for exclusion is ',num2str(THR_FOR_CONNECTION),' %%']},1)

            %% Show mean correlation with each ROI
            max_corr_for_cell       = max(events.globality_index(2:end)); % QQ 1st point sometimes show some artifacts
            for key = 1:numel(ROIs)
                corr_results_sub    = correlation_res.corr_results(events.t_corr(events.is_global), correlation_res.comb(:,2) == key | correlation_res.comb(:,1) == key);
                corr_results_sub    = [corr_results_sub(:,1:(key-1)), NaN(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];
                mean_corr           = nanmean(corr_results_sub,1);
                mean_corr_with_others(key) = sum(mean_corr > THR_FOR_CONNECTION); %nanmean(mean_corr);%
            end

            %% Normalize to 100% (i.e. all ROIs)
            mean_corr_with_others = mean_corr_with_others / numel(mean_corr_with_others); % renormalize to max possible corr for this cell

            %% Renormalize to max possible corr for this cell
            mean_corr_with_others_norm = mean_corr_with_others / max(mean_corr_with_others); 
            
            if obj.rendering
                figure(88888);cla();hist(100*mean_corr_with_others,0:2:100); hold on; 
                title('% of correlation with all other ROIs');set(gcf, 'Color','w')
                xlabel('% of correlation'); ylabel('counts')
                
                figure(88889);cla();hist(100*mean_corr_with_others_norm,0:2:100); hold on;
                title('% of correlation with all other ROIs (Normalized to max)');set(gcf, 'Color','w')
                xlabel('% of correlation'); ylabel('counts')
            end
            
            %% Update the obj.bad_ROI_thr field
            if obj.bad_ROI_thr ~= THR_FOR_CONNECTION
                obj.bad_ROI_thr = THR_FOR_CONNECTION;
            end            
            bad_ROIs            = mean_corr_with_others_norm < obj.bad_ROI_thr;            
            obj.bad_ROI_list    = bad_ROIs; % below threshold % of max correlation
            bad_ROIs            = find(bad_ROIs);
            
            %% ROIs that were manually excluded
            was_excluded        = ismember(obj.ref.indices.swc_list(:,4), obj.batch_params.excluded_branches);            
            RECOVERY_THR        = 1 - THR_FOR_CONNECTION;
            excl_but_not_bad    = was_excluded & ~obj.bad_ROI_list';
            excl_but_good       = was_excluded & (mean_corr_with_others_norm > RECOVERY_THR)';
            if any(excl_but_good)
                 disp_info(['!!! ROIs ',num2str(find(excl_but_good')),' was/were excluded but seem highly correlated\n'],2)
            end
                            
            if obj.rendering     
                %% Get the bad traces
                bad_traces          = obj.extracted_traces_conc(:, find(obj.bad_ROI_list));
                bad_traces          = bad_traces - prctile(bad_traces, 1);
                
                %% Get the reference trace 
                reference_trace     = obj.global_median_raw;
                reference_trace     = reference_trace - prctile(reference_trace, 1);
                
                %% Get traces that we may want to recover
                recoverable         = obj.extracted_traces_conc(:, find(excl_but_not_bad));
                recoverable         = recoverable - prctile(recoverable, 1);
                
                figure(1031);clf();subplot(1,2,1);set(gcf, 'Color','w');
                title(['Bad ROIs (NEVER above ',num2str(obj.bad_ROI_thr*100),' % correlation with the rest of the tree)']);
                plot(smoothdata(reference_trace,'gaussian',[20,0]),'k'); hold on;
                plot(smoothdata(bad_traces,'gaussian',[20,0]),'r');hold on;
                plot(smoothdata(recoverable,'gaussian',[20,0]),'b');
                
                %% Plot normalized excluded traces
                ax = subplot(1,2,2);
                color_code = repmat([0.5,0.5,0.5], obj.n_ROIs, 1);
                color_code(obj.bad_ROI_list,:) = repmat([1,0,0], sum(obj.bad_ROI_list), 1);
                excl = ~ismember(obj.ref.indices.swc_list(:,1), obj.ref.indices.valid_swc_list(:,1));
                color_code(excl,:) = repmat([0.8,0.8,0.8], sum(excl), 1);
                plot_many_traces(smoothdata(obj.extracted_traces_conc,'gaussian',[20,0]), ax);
                colororder(ax, color_code);
                
                %% Plot location of excluded traces
                f = obj.ref.plot_value_tree(obj.bad_ROI_list, 1:numel(obj.bad_ROI_list), obj.default_handle, 'Uncorrelated ROIs', '',  1032,'','redblue'); hold on;
                if any(excl_but_not_bad)
                    obj.ref.plot_value_tree(repmat(0.7,1,sum(excl_but_not_bad)), find(excl_but_not_bad), obj.default_handle, 'Uncorrelated ROIs', '',  f(1).Parent); hold on;
                end
                %                 recovered = ((excl | obj.bad_ROI_list') & ~excl_but_good);
                %                 if any(recovered)
                %                     obj.ref.plot_value_tree(repmat(0.2,1,sum(recovered)), find(recovered), obj.default_handle, 'Uncorrelated ROIs', '',  f.Parent,'','RedBlue'); hold on;
                %                 end
                caxis([0,1]); % otherwise if all values are the same you get a white tree on a white bkg
                
                %% Plot (if possible) the excluded traces, and group them by activity pattern if possible
                regroup_traces(bad_traces', 80, 'pca')  
            end
            obj.bad_ROI_list = find((was_excluded | obj.bad_ROI_list'));
        end

        function [corr_results, comb] = get_pairwise_correlations(obj, idx_filter, corr_window)
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
        
        function [global_median, all_traces_per_bin] = set_median_traces(obj, use_rescaled)
            if (nargin < 2 || isempty(use_rescaled) || use_rescaled) && ~isempty(obj.rescaling_info)
                traces          = obj.rescaled_traces;
                obj.is_rescaled = true;
            else
                traces          = obj.extracted_traces_conc;
                obj.is_rescaled = false;
            end

            %% Create median trace per bins
            all_traces_per_bin = cell(1, numel(obj.binned_data.groups));
            for gp = 1:numel(obj.binned_data.groups)
                columns                     = ismember(obj.binned_data.readmap, obj.binned_data.groups{gp}) & ~ismember(obj.binned_data.readmap, obj.bad_ROI_list);
                if obj.use_hd_data
                     columns = obj.get_voxel_for_ROI(find(columns));
                end

                all_traces_per_bin{gp}      = nanmedian(traces(:,columns), 2);
            end

            global_median                   = nanmedian(traces, 2);
            obj.binned_data.global_median   = global_median;
            all_traces_per_bin              = cell2mat(all_traces_per_bin);
            obj.binned_data.median_traces   = all_traces_per_bin;

            if obj.rendering
                obj.plot_median_traces(obj.is_rescaled);arrangefigures([1,2]);
            end
            if obj.is_rescaled
                obj.disp_info('obj.binned_data.median traces were computed using rescaled signals',1)
            else
                obj.disp_info('obj.binned_data.median traces were computed using raw signals',1)
            end
        end
    end
end

