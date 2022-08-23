classdef behaviours_analysis < handle
    %% Subclass of arboreal_scan_experiment 
    properties
        behaviours                      % List of available behaviours
        detrend_behaviour   = false;    % If true, behaviours are detrended before applying threshold
        detrend_win         = 100       % Defines a moving min subtraction of the behaviour. Only if detrend_behaviour is true        
        beh_thr             = 10        % Threshold in Percent of Max behavioural value, after detrending
        bout_extra_win      = [3, 3]    % Enlarge bouts windows by [before, after] seconds. 
        
    end

    methods
        function behaviours = get.behaviours(obj)    
            %% Get method for the behaviours variables
            % -------------------------------------------------------------
            % Syntax:
            %   behaviours = obj.behaviours
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   behaviours (STRUCT)
            %   behavioural information with the following fields:
            %       * external_var : a copy of obj.external_var
            %       * valid_encoder : a list of booleans indicating encoder
            %           data availability 
            %       * valid_mc_log : a list of booleans indicating movement
            %           correction log availability 
            %       * types : a 1xN cell array of available behaviours
            %       * valid_behaviours : a PxN validity boolean, with N 
            %           the behaviours listed in behaviours.types and P the
            %           number of recordings 
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            
            behaviours                  = obj.behaviours;
            behaviours.external_var     = obj.external_variables;
            behaviours.valid_encoder    = cellfun(@(x) ~isempty(x.encoder.time), behaviours.external_var);
            behaviours.valid_mc_log     = cellfun(@(x) isfield(x,'RT3D_MC'), behaviours.external_var);
            [Max_var, Max_var_loc]      = max(cellfun(@(x) numel(fieldnames(x)), behaviours.external_var));
            behaviours.types            = fieldnames(behaviours.external_var{Max_var_loc})';
            behaviours.valid_behaviours = false(numel(behaviours.valid_encoder), Max_var);
            t_starts_no_gap_real        = cellfun(@(x) x.analysis_params.t_starts - cumsum(x.analysis_params.inter_trials), obj.arboreal_scans, 'UniformOutput', false);
            behaviours.triggers         = cellfun(@(x, y) y + x.header.TTL_delay, obj.arboreal_scans, t_starts_no_gap_real, 'UniformOutput', false);
            behaviours.triggers(cellfun(@(x) ~any(x), behaviours.triggers)) = {[]};
            
            
            for beh = 1:Max_var
                behaviours.valid_behaviours(:,beh) = cellfun(@(y) isfield(y, behaviours.types{beh}), behaviours.external_var)' & cellfun(@(y) ~isempty(y.(behaviours.types{beh}).value), behaviours.external_var)' & cellfun(@(y) ~all(isnan(y.(behaviours.types{beh}).value(:))), behaviours.external_var)';
            end
            if isa(obj.detrend_behaviour, 'function_handle')
                for rec = 1:numel(behaviours.external_var)
                    for beh = 1:Max_var
                        temp = behaviours.external_var{rec}.(behaviours.types{beh}).value;
                        behaviours.external_var{rec}.(behaviours.types{beh}).value = temp - obj.detrend_behaviour(temp);
                    end
                end
            end
        end
        
        function [raw_beh, downsamp_beh, concat_downsamp_beh] = get_behaviours(obj, type, rendering, detrend_sig)
            %% Return the selected behaviour and corresponding timescale
            % -------------------------------------------------------------
            % Syntax:
            %   [raw_beh, downsamp_beh, concat_downsamp_beh] =
            %           EXPE.get_behaviours(type, rendering, detrend_sig)
            % -------------------------------------------------------------
            % Inputs:
            %   type (CHAR OR CELL ARRAY OF CHAR)
            %       Set the behaviour(s) to load. For more than one
            %       behaviour, use a cell array. Any behaviour listed in
            %       obj.behaviours.types will be selected. The filtering 
            %       uses the contains() function (not strcmp) and is case
            %       insensitive. 
            %   rendering (BOOL)
            %       If true, display the selected behaviours
            %   detrend_sig (BOOL, INT or function_handle) - Optional -
            %       Default is false
            %       * If true, a moving min value of 1s is removed from the 
            %       entire signal. 
            %       * If Int, a moving min value of the set number of point 
            %       is removed from the entire signal 
            %       * If function_handle, the function is applied to every
            %       recording and every behaviour. see get.beaviours  
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_traces (1xN CELL ARRAY of 1xP CELLS ARRAY of STRUCT)
            %       For each N behaviours, a structure of P recordings with 
            %       a structure containing a value and a time field, 
            %       as extracted in
            %       arboreal_scan.analysis_params.external_var
            %   downsamp_beh (1xN CELL ARRAY of 1xP CELLS ARRAY of STRUCT)
            %       For each N behaviours, a structure of P recordings with 
            %       a structure containing a value and a time field,
            %       downsampled to match the imaging data
            %   concat_downsamp_beh (STRCUT of NxT timepoints)
            %       For each N behaviour, the concatenated behaviour and
            %       timescale
            % -------------------------------------------------------------
            % Extra Notes:
            %  * If you make a typo in the behaviour selection and at least 
            %   one variable is returned, you won't be informed
            %  * Behaviours are extracted from individual
            %    arboreal_scan.analysis_params.external_var. see
            %    get.behaviours doc
            %  * Default detrending is a min
            %   [~, ~, ori] = obj.get_behaviours('EyeCam_R_forelimb', true, true)
            %   [~, ~, detrend] = obj.get_behaviours('EyeCam_R_forelimb', true, false)
            %   figure();plot(ori.value'); hold on; plot(detrend.value')
            %  * If a behaviour has multiple variables, the average is 
            %    returned. You can change that section of the code
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            if nargin < 2 || isempty(type)
                type = 'encoder';
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end            
            if nargin < 4 || isempty(detrend_sig)
                detrend_sig = false;
            end   
            obj.detrend_behaviour = detrend_sig;

            %% Initialize variables
            downsamp_beh                = {};
            concat_downsamp_beh.time    = [];
            concat_downsamp_beh.value   = [];
            temp                        = {};
            
            %% Check if at least one variable name is valid
            all_beh                     = obj.behaviours;
            type                        = all_beh.types(contains(all_beh.types, type));
            if isempty(type)
                raw_beh = {};
                warning(['Type not detected. Valid behaviours are :\n', strjoin(all_beh.types,'\n')]);
                return
            end
            
            %% Extract variables
            for beh = 1:numel(type)
                raw_beh{beh}                     = arrayfun(@(x) x{1}.(type{beh}), all_beh.external_var, 'UniformOutput', false, 'ErrorHandler', @cellerror_empty);

                temp.time  = [];
                temp.value = [];
                for rec = 1:numel(raw_beh{beh})                    
                    %% If behaviour timescale is longer than recording, clip it
                    if ~isempty(raw_beh{beh}{rec}.time)
                        
                        clipping_idx = find((raw_beh{beh}{rec}.time - obj.timescale.durations_w_gaps(rec) ) > 0, 1, 'first');
                        if ~isempty(clipping_idx)
                            raw_beh{beh}{rec}.time = raw_beh{beh}{rec}.time(1:clipping_idx);
                            raw_beh{beh}{rec}.value = raw_beh{beh}{rec}.value(1:clipping_idx);
                        end
                    end                        
                    downsamp_beh{beh}{rec}.time     = interpolate_to(raw_beh{beh}{rec}.time, obj.timescale.tp(rec));
                    downsamp_beh{beh}{rec}.value    = interpolate_to(raw_beh{beh}{rec}.value, obj.timescale.tp(rec));
                    if isempty(downsamp_beh{beh}{rec}.time)
                        downsamp_beh{beh}{rec}.time = linspace(0, obj.timescale.durations(rec), obj.timescale.tp(rec));
                        downsamp_beh{beh}{rec}.value = NaN(1,obj.timescale.tp(rec));
                    end

                    %downsampd_beh{beh}{rec}.value = smoothdata(downsampd_beh{beh}{rec}.value, 'movmean', [5, 0]);
                    temp.time = [temp.time, downsamp_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec)];
                    value = downsamp_beh{beh}{rec}.value;
                    
                    if size(value, 1) > 1 % if your behavioural metrics is made of multiple arrays
                        value = nanmean(value, 1);
                    end
                    temp.value = [temp.value, value];
                end

                concat_downsamp_beh.time = [concat_downsamp_beh.time ;temp.time];
                concat_downsamp_beh.value = [concat_downsamp_beh.value ;temp.value];
            end

            %% Render extracted traces
            if rendering
                figure(1026);clf();
                for beh = 1:numel(type)
                    subplot(numel(type),1,beh)
                    for rec = 1:numel(raw_beh{beh})
                        if ~isempty(downsamp_beh{beh}{rec}.time)                           
                            plot(downsamp_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec), downsamp_beh{beh}{rec}.value);hold on;
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

        function [bouts, beh_sm, active_tp] = get_activity_bout(obj, beh_types, rendering, smoothing, invert, thr)
            %   invert (BOOL) - Optional - Default is false
            %       * If true, the detected behaviours is inverted  
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
            if nargin < 5 || isempty(invert)
                invert = false;
            end 
            if nargin < 6 || isempty(thr)
                % pass
            else
                obj.beh_thr = thr;
            end 
            
            plts = {};
            for idx = 1:numel(beh_types)
                current_type    = beh_types{idx};
                [~, ~, beh]     = obj.get_behaviours(current_type, false);
                if ~isempty(beh.value)
                    if any(obj.detrend_win)
                        beh.value       = beh.value - movmin(beh.value, [obj.detrend_win, 0]);
                    end
                    beh_sm{idx}          = smoothdata(beh.value, 'gaussian', smoothing);
                    beh_sm{idx}          = nanmean(beh_sm{idx},1);
                    %beh_sm{idx}         = detrend(fillmissing(beh_sm{idx},'nearest'),'linear',cumsum(obj.timescale.tp));
                    %thr            = prctile(beh_sm{idx}(beh_sm{idx} > 0), 20);
                    current_thr     = prctile(beh_sm{idx},(100/obj.beh_thr)) + range(beh_sm{idx})/(100/obj.beh_thr); % 5% of max

                    %% Define bouts
                    if ~invert
                        active_tp{idx}       = abs(beh_sm{idx}) > abs(current_thr);
                    else
                        active_tp{idx}       = abs(beh_sm{idx}) < abs(current_thr);
                    end
                    [starts, stops] = get_limits(active_tp{idx}, beh_sm{idx});

                    %% Add some pts before and after each epoch
                    dt = nanmedian(diff(obj.t));
                    for epoch = starts
                        active_tp{idx}(max(1, epoch-round(obj.bout_extra_win(1)*(1/dt))):epoch) = 1;
                    end
                    for epoch = stops
                        active_tp{idx}(epoch:min(numel(active_tp{idx}), epoch+round(obj.bout_extra_win(2)*(1/dt)))) = 1;
                    end
                    [starts, stops] = get_limits(active_tp{idx}, beh_sm{idx}); % update bouts edges now that we extended the range

                    if rendering
                        figure(1027);hold on;
                        if idx == 1
                            clf();
                        end
                        plts{idx} = subplot(numel(beh_types),1,idx);hold on;
                        title(strrep(current_type,'_','\_'))
                        plot(obj.t, beh_sm{idx});hold on;
                        for el = 1:numel(starts)
                            x = [starts(el),starts(el),stops(el),stops(el)];
                            y = [0,nanmax(beh_sm{idx}),nanmax(beh_sm{idx}),0];
                            patch('XData',obj.t(x),'YData',y,'FaceColor','red','EdgeColor','none','FaceAlpha',.1);hold on
                        end
                    end

                    bouts{idx} = sort([starts, stops]);
                end
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
        
        function set.detrend_behaviour(obj, detrend_behaviour)
            if islogical(detrend_behaviour)
                if ~detrend_behaviour
                    obj.detrend_behaviour = false;
                else
                    sr = obj.timescale.sr;
                    obj.detrend_behaviour = @(x) movmin(x, [ceil(1/nanmedian(sr)), 0]);
                end                
            elseif isnumeric(detrend_behaviour)
                sr = obj.timescale.sr;
                obj.detrend_behaviour = @(x) movmin(x, [ceil(detrend_behaviour/nanmedian(sr)), 0]);
            elseif isa(detrend_behaviour, 'function_handle')
                obj.detrend_behaviour = detrend_behaviour;
            else
                error('detrend_behaviour must be a boolean, a value in seconds or a function handle')
            end
        end
        
        function set.bout_extra_win(obj, bout_extra_win)
            if numel(bout_extra_win) == 1
                bout_extra_win = [bout_extra_win, bout_extra_win];
            end
            obj.bout_extra_win = bout_extra_win;
        end
        
        function detrend_win = get.detrend_win(obj)
            detrend_win = double(obj.detrend_behaviour) * obj.detrend_win;
        end
        
        function [stim_time, stim_pt] = get_stim_epochs(obj, rendering)
            if nargin < 2 || isempty(rendering)
                rendering = false;
            end

            stim_time   = cellfun(@(x, st) x + st, obj.behaviours.triggers, num2cell(obj.timescale.t_start_nogap),'UniformOutput',false);
            stim_pt     = cellfun(@(x, st, sr) floor(x/sr) + st, obj.behaviours.triggers, num2cell(cumsum([0, obj.timescale.tp(1:end-1)])), num2cell(obj.timescale.sr),'UniformOutput',false);

            if ~isempty(rendering) && any(rendering)
                if ismatrix(rendering) && numel(rendering) == numel(obj.t)
                    traces  = rendering;
                    m       = mode(traces(:));
                else
                    traces  = obj.extracted_traces_conc;
                    m       = nanmin(nanmedian(traces,2));   
                    traces  = nanmedian(traces,2);
                end
                figure();plot(obj.t, traces);hold on; scatter([stim_time{:}], repmat(m, size([stim_time{:}])), 'k^', 'filled');
            end
        end
        
    end
end

