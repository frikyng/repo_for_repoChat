classdef behaviours_analysis < handle
    %% Subclass of arboreal_scan_experiment
    properties
        behaviours                      % List of available behaviours
        detrend_behaviour   = false;    % If true, behaviours are detrended using movmin function before applying threshold
        detrend_win         = -100      % Defines a moving min subtraction of the behaviour. Only if detrend_behaviour is true
        beh_thr             = 10        % Threshold in Percent of Max behavioural value, after detrending
        bout_extra_win      = [3, 3]    % Enlarge bouts windows by [before, after] seconds.
        beh_smoothing       = [-1,0]    % The smoothing window for behaviour. Values < 1 are in seconds
        beh_sm_func         = 'gaussian'% any valid smoothdata option. If contains "robust", outlier are filtered before smoothing
        multi_beh_func      = @nanmean  % The function applied to behaviours that contains multiple arrays (along dim 1)
        shuffling_block_size=0         % number of timepoints for the behaviour shuffling window (applied if a behaviour contains "shuffle" in its name
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
            %behaviours.valid_behaviours(:,beh) = cellfun(@(y) isfield(y, behaviours.types{beh}), behaviours.external_var)' & cellfun(@(y) ~isempty(y.(behaviours.types{beh}).value), behaviours.external_var)' & cellfun(@(y) ~all(isnan(y.(behaviours.types{beh}).value(:))), behaviours.external_var)';
%             if isa(obj.detrend_behaviour, 'function_handle') && any(obj.detrend_win)
%                 for rec = 1:numel(behaviours.external_var)
%                     for beh = 1:Max_var
%                         if isfield(behaviours.external_var{rec}, behaviours.types{beh})
%                             temp = behaviours.external_var{rec}.(behaviours.types{beh}).value;
%                             if ~isempty(temp)
%                                 behaviours.external_var{rec}.(behaviours.types{beh}).value = temp - obj.detrend_behaviour(temp);
%                             end
%                         else
%                            % warning(['Missing behaviour "',behaviours.types{beh},'" for recording ',num2str(rec)], 'backtrace' )
%                         end
%                     end
%                 end
%             end
        end
        
        function set.beh_smoothing(obj, beh_smoothing)
            %% Defines the gaussian filtering window for behaviours
            % -------------------------------------------------------------
            % Syntax:
            %   obj.beh_smoothing = beh_smoothing;
            % -------------------------------------------------------------
            % Inputs:
            %   beh_smoothing (FLOAT OR 2x1 FLOAT)
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

            if all(isnumeric(beh_smoothing)) && numel(beh_smoothing) == 2 && all(beh_smoothing == obj.beh_smoothing) && ~any(obj.beh_smoothing < 0)
                %% no change, pass                
            elseif all(isnumeric(beh_smoothing)) && numel(beh_smoothing) == 1 || numel(beh_smoothing) == 2 
                try
                    beh_smoothing(beh_smoothing < 0)  = beh_smoothing(beh_smoothing < 0) * nanmedian(1./obj.timescale.sr);
                end
                beh_smoothing                  = abs(round(beh_smoothing));                
                if numel(beh_smoothing)    == 1                 
                    beh_smoothing = [beh_smoothing, beh_smoothing];
                end
                obj.beh_smoothing              = beh_smoothing;
            else
                error('filter window must be a set of one (for symmetrical gaussian kernel) or 2 (for asymetrical gaussian kernel) values. Values are rounded. Values < 0 indicate a window in points')
            end
        end
        
        function set.shuffling_block_size(obj, shuffling_block_size)
            if all(isnumeric(shuffling_block_size)) && all(shuffling_block_size == obj.shuffling_block_size)
                %% No change, pass                
            elseif all(isnumeric(shuffling_block_size))
                if shuffling_block_size < 0
                    shuffling_block_size = round(abs(shuffling_block_size) * nanmedian(1./obj.timescale.sr));
                else
                    shuffling_block_size = round(shuffling_block_size);                
                end
                obj.shuffling_block_size              = shuffling_block_size;
            else
                error('filter window must be one value < 0 (block size in s) or >0 (block size in points)')
            end
        end

        function [raw_beh, downsamp_beh, concat_downsamp_beh] = get_behaviours(obj, type, rendering, detrend_sig, ignorecase, average, smoothing)
            %% Return the selected behaviour and corresponding timescale
            % -------------------------------------------------------------
            % Syntax:
            %   [raw_beh, downsamp_beh, concat_downsamp_beh] =
            %           EXPE.get_behaviours(type, rendering, detrend_sig,
            %           ignorecase, average)
            % -------------------------------------------------------------
            % Inputs:
            %   type (CHAR OR CELL ARRAY OF CHAR)
            %       Set the behaviour(s) to load. For more than one
            %       behaviour, use a cell array. Any behaviour listed in
            %       obj.behaviours.types will be selected. The filtering
            %       uses the contains() function (not strcmp) and is case
            %       insensitive.
            %   rendering (BOOL) - Optional - Default is false
            %       If true, display the selected behaviours
            %   detrend_sig (BOOL, INT or function_handle) - Optional -
            %       Default is false
            %       * If true, a moving min value of 1s is removed from the
            %       entire signal.
            %       * If Int, a moving min value is subtracted from the
            %       entire signal. Use a value > 0 to set the moving windo
            %       in points, and < 0 to set a moving window in s
            %       * If function_handle, the function is applied to every
            %       recording and every behaviour. see get.beaviours
            %   ignorecase (BOOL) - Optional - Default is false
            %       If true, the case in the behaviour names is ignored
            %   average (BOOL) - Optional - Default is false
            %       If true, all requested beahviours are averaged
            %   smoothing (INT or 2x1 INT) - Optional - Default is 
            %       obj.beh_smoothing.
            %       * If > 0, defines a smoothing window in points
            %       * If < 0 defines a smoothing window in
            %       seconds. Default smoothing is a gaussian filter, 
            %       unless you change obj.beh_sm_func. If 2x1 INT,
            %       an asymmetrical filter is set.
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
            %  * To plot all behaviours
            %    detrend = true;smoothing = false
            %    obj.get_behaviours({''},true,detrend,'','','',smoothing)
            %  * If you make a typo in the behaviour selection and at least
            %    one variable is returned, you won't be informed
            %  * Behaviours are extracted from individual
            %    arboreal_scan.analysis_params.external_var. see
            %    get.behaviours doc
            %  * Default detrending is a min
            %   [~, ~, ori] = obj.get_behaviours('EyeCam_R_forelimb', true, true)
            %   [~, ~, detrend] = obj.get_behaviours('EyeCam_R_forelimb', true, false)
            %   figure();plot(ori.value'); hold on; plot(detrend.value')
            %  * If a behaviour has multiple variables, the *average* is
            %    returned. You can change that section of the code
            %  * obj.beh_smoothing is applied to every behaviour unless you
            %   specify a value
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            if nargin < 2 || isempty(type)
                type = 'encoder';
            else
                if ischar(type)
                    type = cleanup_type(type);
                elseif iscell(type)
                    type = cellfun(@(x) cleanup_type(x), type , 'UniformOutput', false);
                end
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end            
            if nargin >= 4 && ~isempty(detrend_sig)
                initial_detrend         = obj.detrend_behaviour;  % backup initial detrending
                obj.detrend_behaviour   = detrend_sig;    % handle all the special cases     
                cleanupObj = onCleanup(@() restore_detrend(obj, initial_detrend));
            end
            if nargin < 5 || isempty(ignorecase)
                ignorecase = false;
            end
            if nargin < 6 || isempty(average)
                average = false;
            end
            if nargin < 7 || isempty(smoothing)
                smoothing = obj.beh_smoothing;
            else
                bkp_smoothing = obj.beh_smoothing;
                obj.beh_smoothing = smoothing;
                smoothing = obj.beh_smoothing;
                obj.beh_smoothing = bkp_smoothing;
            end

            %% Initialize variables
            downsamp_beh                = {};
            concat_downsamp_beh.time    = [];
            concat_downsamp_beh.value   = [];
            temp                        = {};

%             %% Identify if shuffling will be needed
%             shuffle_size                = double(shuffle_size * any(contains(type, 'shuffle')));
%             type                        = erase(type,  {'shuffled',' shuffled','shuffled ','_shuffled','shuffled_',' shuffle','shuffle ','_shuffle','shuffle_','shuffle'});
%             
            %% Check if at least one variable name is valid
            all_beh                     = obj.behaviours;
            type                        = all_beh.types(contains(all_beh.types, type, 'IgnoreCase', ignorecase));
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
                        clipping_idx        = find((raw_beh{beh}{rec}.time - obj.timescale.durations_w_gaps(rec) ) > 0, 1, 'first');
                        if ~isempty(clipping_idx)
                            raw_beh{beh}{rec}.time  = raw_beh{beh}{rec}.time(1:clipping_idx);
                            raw_beh{beh}{rec}.value = raw_beh{beh}{rec}.value(1:clipping_idx);
                        end
                    end
                    downsamp_beh{beh}{rec}.time     = interpolate_to(raw_beh{beh}{rec}.time, obj.timescale.tp(rec));
                    downsamp_beh{beh}{rec}.value    = interpolate_to(raw_beh{beh}{rec}.value, obj.timescale.tp(rec));
                    if isempty(downsamp_beh{beh}{rec}.time)
                        downsamp_beh{beh}{rec}.time     = linspace(0, obj.timescale.durations(rec), obj.timescale.tp(rec));
                        downsamp_beh{beh}{rec}.value    = NaN(1,obj.timescale.tp(rec));
                    end
                    if any(smoothing)
                        beh_sm_func = obj.beh_sm_func;
                        if contains(beh_sm_func, 'robust')                            
                            beh_sm_func = erase(beh_sm_func, {'_','robust'});
                            downsamp_beh{beh}{rec}.value = filloutliers(downsamp_beh{beh}{rec}.value, 'pchip','movmedian',[100,0]);
                            %figure();plot(downsamp_beh{beh}{rec}.value);hold on;plot(filloutliers(downsamp_beh{beh}{rec}.value, 'linear','movmedian',[100,0]))
                        end
                        downsamp_beh{beh}{rec}.value    = smoothdata(downsamp_beh{beh}{rec}.value, beh_sm_func, smoothing);
                    end
                    if isa(obj.detrend_behaviour, 'function_handle') && any(obj.detrend_win)
                        %baseline_estimate(downsamp_beh{beh}{rec}.value', obj.detrend_win(1)/10, obj.detrend_win(1))
                        downsamp_beh{beh}{rec}.value = downsamp_beh{beh}{rec}.value - obj.detrend_behaviour(downsamp_beh{beh}{rec}.value);
                    end
                    
                    temp.time               = [temp.time, downsamp_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec)];
                    value                   = downsamp_beh{beh}{rec}.value;

                    if size(value, 1) > 1 % if your behavioural metrics is made of multiple arrays
                        value                   = obj.multi_beh_func(value); 
                    end
                    temp.value              = [temp.value, value];
                end

                concat_downsamp_beh.time    = [concat_downsamp_beh.time     ;temp.time];
                concat_downsamp_beh.value   = [concat_downsamp_beh.value    ;temp.value];
            end
            
            %% Average of all selected behaviours
            if average
                new = {};
                for rec = 1:numel(downsamp_beh{1})
                    new{rec}.value = nanmean(cell2mat(cellfun(@(x) x{rec}.value', downsamp_beh, 'UniformOutput', false)), 2)';
                    new{rec}.time = downsamp_beh{1}{rec}.time;
                end
                downsamp_beh        = {new};
                concat_downsamp_beh = structfun(@(x) nanmean(x,1), concat_downsamp_beh, 'UniformOutput', false);
                type                = {strjoin(type)};
            end

            %% Plot extracted behaviours
            if rendering
                figure(1026);clf();
                ax_list = {};
                for beh = 1:numel(type)                 
                    ax_list{end+1} = subplot(numel(type),1,beh);
                    for rec = 1:numel(raw_beh{beh})
                        if ~isempty(downsamp_beh{beh}{rec}.time)
                            plot(downsamp_beh{beh}{rec}.time + obj.timescale.t_start_nogap(rec), downsamp_beh{beh}{rec}.value);hold on;
                        end
                    end
                end
                linkaxes([ax_list{:}], 'x')
            end

            function out = cellerror_empty(~,varargin) % missing behaviour management
                out = {};
                out.time = [];
                out.value = [];
            end
            
            function type = cleanup_type(type)
                type                = erase(type, '~');
                range_idx           = strfind(type, '[');
                type(range_idx:end) = [];
            end

            function restore_detrend(obj, initial_detrend)
                obj.detrend_behaviour   = initial_detrend;  % restore initial default option
            end
        end

        function [bouts, beh_sm, active_tp] = get_activity_bout(obj, beh_types, rendering, invert, thr)
            %% Return the selected behaviour and corresponding timescale
            % -------------------------------------------------------------
            % Syntax:
            %   [bouts, beh_sm, active_tp] =
            %       EXPE.get_activity_bout(
            %       beh_types, rendering, smoothing, invert, thr, window)
            % ------------------------------------------------------------- 
            % Inputs:
            %   type (CHAR OR CELL ARRAY OF CHAR) - Optional - default is
            %   'encoder'
            %       Set the behaviour(s) to load. For more than one
            %       behaviour, use a cell array. Any behaviour listed in
            %       obj.behaviours.types will be selected. The filtering
            %       uses the contains() function (not strcmp) and is case
            %       insensitive.
            %   rendering (BOOL) - Optional - Default is false
            %       If true, display the selected behaviours
            %   invert (BOOL) - Optional - Default is false
            %       If true, the detected behaviours is inverted
            %   thr (FLOAT) - Optional - Default is obj.beh_thr. 
            %       Defines the threshold used to start a bout. The
            %       threshold is defined as a value exceeding obj.beh_thr
            %       of the signal soft range (i.e. on the 1st-99th
            %       percentile range).
            % -------------------------------------------------------------
            % Outputs:
            %   bouts (1xN CELL ARRAY of 1x2P ARRAY)
            %       For each N behaviour, the start and stop point of
            %       the P activity bouts
            %   beh_sm (1xN CELL ARRAY of 1xP FLOAT)
            %       For each N behaviours, the corresponding behaviour
            %       (smoothed if smoothing was used)
            %   active_tp (1xN CELL ARRAY of 1xP LOGICAL)
            %       For each N behaviours, an array of boolean indicating
            %       if a given timepoint is active or inactive.
            % -------------------------------------------------------------
            % Extra Notes:
            %   * A "window" flag can be added to the behaviour (eg
            %       "encoder" becomes "encoder[10,-3]". The flag format is 
            %       (INT or 2x1 or 3x3 INT). If absent the value used is
            %       obj.bout_extra_win.
            %       This is an expension/erosion factor for each bout. See 
            %       Extra notes for some examples. If you provide a single
            %       INT, the factor is applied to the beginning and to the 
            %       end of each bout. Two value control independently the
            %       expension/shrinkage of the onset and offset of each bout.
            %       A 3rd value force the bout duration to this value
            %   "window" can be use to manipulate the limits of the bouts.
            %   - If you specifies window = 0, the bouts are unchanged. 
            %   - If you specified 3, the beginning and end of each bout is
            %   extended by that many seconds. If 2 bouts touch each other, 
            %   they are merged
            %   - If you specifies -3, the bouts are eroded on each side.
            %   If the bout if fully eroded, it is removed.
            %   - You could specifies [-10, 20] To ignore the first 10
            %   points of a bout, but study the 20 points following the end
            %   - Specifying a symmetrical window, but with inverted sign
            %   will introduce a lag. For example [-5, 5] study the data
            %   follwing bout onset by 5 seconds, until 5 second after 
            %   offset . [5, -5] study the data preceding the bout by 5
            %   seconds, until 5 seconds before the end of the bout.
            %   - If you specify a 3rd value, it forced the bout duration.
            %   Positive value align it to the onset, negative value align
            %   it to the offset. For example :
            %       * [0,0,10] create bouts of 10 points aligned on the
            %       onset.
            %       * [0,0, -10] create bouts of 10 points aligned on the
            %       offset.
            %       * [-5,0, 10] would create bout of 5 points preceding
            %       the onset by 10 points (i.e. it would detecting a lag).
            %       Note that in this case the second value has no affect
            %       on the offset
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   17/08/2022

            if nargin < 2 || isempty(beh_types)
                beh_types = {'encoder'};
            elseif ischar(beh_types)
                beh_types = {beh_types};
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end
            if nargin < 4 || isempty(invert)
                invert = contains(beh_types, '~');   
            elseif numel(invert) == 1
                invert = repmat(invert, 1, numel(beh_types));
            end
            if nargin < 5 || isempty(thr)
                % pass
            else
                obj.beh_thr = thr;
            end 
            
            plts = {};
            for idx = 1:numel(beh_types)
                current_type    = beh_types{idx};
                [~, ~, beh]     = obj.get_behaviours(current_type, false);
                
                %% Check if there is window suffix
                analysis_mode_no_win = strsplit(current_type,{'[',']'});
                analysis_mode_no_win = analysis_mode_no_win(1:3:end);
                to_test = [];
                for el = obj.behaviours.types
                    if contains(analysis_mode_no_win, el{1},'IgnoreCase',true) || contains(el{1},analysis_mode_no_win,'IgnoreCase',true)
                        to_test(end+1) = 1;
                    else
                        to_test(end+1) = 0;
                    end
                end
                beh_name = {obj.behaviours.types{find(to_test)}};
                
                
                beh_end_loc = strfind(current_type, beh_name{1}) + numel(beh_name{1});
                if ~isempty(beh_end_loc) && beh_end_loc < numel(current_type) && strcmp(current_type(beh_end_loc), '[')
                    loc             = strfind(current_type,']');
                    loc             = loc(find(loc > beh_end_loc, 1, 'first'));
                    window          = str2num(current_type(beh_end_loc:loc));
                else
                    window          = obj.bout_extra_win;
                end
                
                if ~isempty(beh.value)
                    beh_sm{idx}         = beh.value;
                    beh_sm{idx}         = nanmean(beh_sm{idx},1);
                    %beh_sm{idx}         = detrend(fillmissing(beh_sm{idx},'nearest'),'linear',cumsum(obj.timescale.tp));
                    
                    %% Define threshold. We used a "soft range" (the signal 1-99th percentile range) on the absolute value to determine a threshold
                    temp_beh            = abs(beh_sm{idx});
                    noise               = range(temp_beh)/100;
                    temp_beh            = abs(temp_beh + randn(size(temp_beh))*noise);                    
                    current_thr         = prctile(temp_beh,10) + ((prctile(temp_beh,99) - prctile(temp_beh,1)) * (obj.beh_thr/100));   %prctile(beh_sm{idx},(100/obj.beh_thr)) + range(beh_sm{idx})/(100/obj.beh_thr); % obj.beh_thr 'th percentile of the behaviour + obj.beh_thr % of the range

                    %% Define bouts starts and stops
                    dt                  = nanmedian(diff(obj.t));
                    active_tp{idx}      = false(size(beh_sm{idx}));
                    [starts, stops]     = get_limits(abs(beh_sm{idx}) > abs(current_thr), beh_sm{idx});
                    
                    if numel(window) == 3 && window(3) == 0
                        obj.disp_info('You set 3 values for the window control but the third value is 0. This will select only one timepoint. Make sure this is the intended behaviour',3)
                    end
                    if numel(window) == 3 && window(3) < 0 && window(1) ~=0
                        obj.disp_info('You set a negative value for the window control setting #3, which means you want to use events occuring before/after the OFFSET. However, you still have a setting# 1 value that will expend/shrink the ONSET. Make sure this is the intended behaviour',3)
                    end
                    if numel(window) == 3 && window(3) > 0 && window(2) ~=0
                        obj.disp_info('You set a positive value for the window control setting #3, which means you want to use events occuring before/after the ONSET. However, you still have a setting# 2 value that will expend/shrink the OFFSET. Make sure this is the intended behaviour',3)
                    end
                    
                    
                    %% Extend the bouts limits
                    if window(1) > 0
                        starts              = starts - round(window(1)*(1/dt));
                        starts(starts < 1)  = 1; starts(starts > numel(active_tp{idx})) = numel(active_tp{idx});
                    end
                    if window(2) > 0
                        stops               = stops + round(window(2)*(1/dt));
                        stops(stops < 1)    = 1; stops(stops > numel(active_tp{idx})) = numel(active_tp{idx});
                    end
                    [starts, stops] = merge_epochs(starts, stops);
                    
                    %% Now shrink bouts limits
                    if window(1) < 0
                        starts              = starts - round(window(1)*(1/dt));
                        starts(starts < 1)  = 1; starts(starts > numel(active_tp{idx})) = numel(active_tp{idx});
                    end
                    if window(2) < 0
                        stops               = stops + round(window(2)*(1/dt));
                        stops(stops < 1)    = 1; stops(stops > numel(active_tp{idx})) = numel(active_tp{idx});
                    end
                    [starts, stops] = merge_epochs(starts, stops);
                    obj.disp_info(['N bouts after extension / shrinkage and merging = ',num2str(numel(starts))],1)

                    %% Force bouts duration if required
                    if numel(window) == 3
                        if window(3) >= 0
                            stops = starts + round(window(3)*(1/dt));
                            window(2) = 0;
                            obj.disp_info(['all detected bouts adjusted to a duration of = ',num2str(numel(abs(window(3))*(1/dt))),' seconds before/after onset'],1)
                        elseif window(3) < 0
                            starts = stops + round(abs(window(3))*(1/dt));
                            window(1) = 0;
                            obj.disp_info(['all detected bouts adjusted to a duration of = ',num2str(numel(abs(window(3))*(1/dt))),' seconds before/after offset'],1)
                        end                        
                    end

                    %% Build range                    
                    for epoch = 1:numel(starts)
                        if numel(window) == 3
                            pts = sort([starts(epoch), stops(epoch)]); % in some complicated lag + fixed windo cases, the acceptance windo can be inverted
                            starts(epoch) = pts(1);stops(epoch) = pts(2);
                        end
                        R = (starts(epoch)-round(window(1)*(1/dt))):(stops(epoch)+round(window(2)*(1/dt)));
                        R(R < 1) = []; R(R > numel(active_tp{idx})) = [];
                        active_tp{idx}(R) = 1;
                    end
%                     for epoch = starts
%                         active_tp{idx}(max(1, epoch-round(window(1)*(1/dt))):epoch) = 1;
%                     end
%                     for epoch = stops
%                         active_tp{idx}(epoch:min(numel(active_tp{idx}), epoch+round(window(2)*(1/dt)))) = 1;
%                     end
                    
                    %% If required, invert
                    if invert(idx)
                        active_tp{idx}       = ~active_tp{idx};
                    end    
                    
%                     %% Update bouts edges now that we extended the range
%                     [starts, stops] = get_limits(active_tp{idx}, beh_sm{idx}); 

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
                %% INTERNAL FUNCTION
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
            
            function [merged_starts, merged_stops] = merge_epochs(starts, stops)
                % Combine starts and stops into a single array with tags
                timepoints = [starts(:); stops(:)];
                tags = [ones(size(starts(:))); 2*ones(size(stops(:)))]; % 1 for start, 2 for stop

                % Sort timepoints and tags
                [sorted_timepoints, sort_idx] = sort(timepoints);
                sorted_tags = tags(sort_idx);

                % Initialize output arrays and other variables
                merged_starts = [];
                merged_stops = [];
                active_epochs = 0;
                current_start = NaN;

                % Loop through sorted timepoints
                for i = 1:length(sorted_timepoints)
                    if sorted_tags(i) == 1 % start
                        if active_epochs == 0
                            current_start = sorted_timepoints(i);
                        end
                        active_epochs = active_epochs + 1;
                    elseif sorted_tags(i) == 2 % stop
                        active_epochs = active_epochs - 1;
                        if active_epochs == 0
                            merged_starts = [merged_starts; current_start];
                            merged_stops = [merged_stops; sorted_timepoints(i)];
                        end
                    end
                end
            end

            
            
        end

        function set.detrend_behaviour(obj, detrend_behaviour)
            if isa(detrend_behaviour, 'function_handle')
                % keep as such
            elseif islogical(detrend_behaviour) && detrend_behaviour
                obj.detrend_win     = ceil(100/nanmedian(obj.timescale.sr));
            	detrend_behaviour   = @(x) movmin(x, [obj.detrend_win, 0]);
            elseif islogical(detrend_behaviour) || ~any(detrend_behaviour)
                % keep as such
            elseif isnumeric(detrend_behaviour) && any(detrend_behaviour > 0)
                obj.detrend_win     = detrend_behaviour(1);
                detrend_behaviour   = @(x) movmin(x, [obj.detrend_win, 0]);
            elseif isnumeric(detrend_behaviour) && any(detrend_behaviour < 0)
                obj.detrend_win     = ceil(abs(detrend_behaviour(1))/nanmedian(obj.timescale.sr));
                detrend_behaviour   = @(x) movmin(x, [obj.detrend_win, 0]);
            else
                error('detrend_behaviour must be a boolean, a value in seconds or a function handle')
            end
            obj.detrend_behaviour = detrend_behaviour;
        end

        function set.bout_extra_win(obj, bout_extra_win)
            if numel(bout_extra_win) == 1
                bout_extra_win = [bout_extra_win, bout_extra_win];
            end
            obj.bout_extra_win = bout_extra_win;
        end
    end
end
