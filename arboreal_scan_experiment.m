
%% TO DO :
% - improve update system when changing folder.
% - enable loading if we just have a folder with arboreal_scans objects
% - add warning if source folde ris the actual raw experiment
% - add doc to load_Several_Experiment

%% TO CHECK:
% - loading from raw data
% - Breakpoints not iintorucing artefacts
% - Check if bheviours are fine
% why are camera behaviours downsampled. Do we want that?

classdef arboreal_scan_experiment < handle & arboreal_scan_plotting & event_fitting & behaviours_analysis
    properties
        %% Extraction/Re-extraction Settings
        extraction_method = 'median'; % the 1D -> 0D compression method. If None, XxT data is kept, but this make the files heavier
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
        dim_red_type    = 'pca'
        peak_thr        = 2;
        bad_ROI_thr     = 0;
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
        spiketrains             % If available, spike inference results
        bad_ROI_list     = 'unset';  % list of uncorrelated ROIs (following event detection)
        updatable               % If arboreal_scan are still available, you could update the arboreal_scan_experiment compression
        crosscorr               % Correlation of peaks/signal across ROIs/groups during/between bAps/activity_bouts
    end

    properties (Dependent = true, Transient = true)
        updated_data_path            % If update_folder is used, the updated filpath
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
        function obj = arboreal_scan_experiment(source_folder, keep_2D, varargin)
            %% arboreal_scan_experiment Constructor
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE = 
            %   arboreal_scan_experiment(source_folder, keep_2D, varargin)
            % -------------------------------------------------------------
            % Inputs:
            %   source_folder (STR)
            %       Path to the raw or extracted recordings. It can be :
            %           * a folder containing extracted arboreal_scans
            %             objects (recommended)
            %           * an experiment folder, in which case data will be
            %             extracted using additional varargin
            %   keep_2D (BOOL) - Optional, default is False
            %       If true, the full_data field is kept for each ROI. Not
            %       that this will considerably increase the object size.
            %   varargin
            %       When passing an experiment folder, varargin{1} can 
            %       contain the analysis_params to use for extraction, and 
            %       varargin{2} the path to a settings.txt file
            % -------------------------------------------------------------
            % Outputs:
            %   obj (arboreal_scan_experiment object)
            %       arboreal_scan_experiment object handle, that contains
            %       individual SIMPLIFIED arboreal_scans (see extra Notes),
            %       and built-in analysis features.
            % -------------------------------------------------------------
            % Extra Notes:
            %   * You can build the object from RAW data every time, but it
            %   is recommended to do it in two steps, in particular if you
            %   have several experiment to process.
            %   * arboreal_scan_experiment objects also contains the
            %   arboreal_scan objects used for extraction BUT to reduce
            %   the file size, each ROI data is compressed into a single 
            %   1xT time serie, while it is a XxT Timeseries in the
            %   arboreal_scan object. The compression is done using 
            %   obj.extraction_method. To keep XxT data, set Keep_2D = true
            %   * See Documentation for examples
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            if nargin < 1
                return
            end
            if nargin < 2 || isempty(keep_2D)                
                keep_2D = false;
            end
            
            %% Fix paths
            obj.source_folder       = parse_paths(source_folder);

            %% Get a list of extracted arboreal scans in the source folder
            obj.extracted_data_paths= list_sources(obj);
            obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update

            %% Load arboreal scans
            obj.update(true, keep_2D, varargin);
        end

        function extracted_data_paths = list_sources(obj)
            %% Return the location of the original data used for extraction
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_data_paths = EXPE.list_sources
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_data_paths (1 x N CELL ARRAY of CHAR)
            %       Each cell contains the path to the original data sued
            %       for extraction
            % -------------------------------------------------------------
            % Extra Notes:
            %   * If data is not available at the original location and not
            %   available at the updated location (if obj.update_folder is
            %   non empty), a warning message is printed. In the absence of
            %   the original data, you cannot :
            %       - re-extract data
            %       - Update the compression method,
            %       - reload the original signal from the analyzed trees
            %   Other analyses should work properly
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            %% List available arboreal scans
            if isempty(obj.update_folder) % normal case
                all_recordings          = dir([obj.source_folder,'/**/*-*-*_exp_*_*-*-*']);
            else
                all_recordings          = dir([obj.update_folder,'/**/*-*-*_exp_*_*-*-*']);
            end
            if isempty(all_recordings)
                warning(['Original arboreal_scans not found in : ',[obj.source_folder,'/**/*-*-*_exp_*_*-*-*'],' . extraction/re-extraction not available. If you moved the files to a new top folder, change obj.update_folder accordingly']);
                extracted_data_paths= obj.extracted_data_paths;
                return
            end
            all_recordings          = all_recordings(~[all_recordings(:).isdir]);
            all_recordings          = all_recordings(~(arrayfun(@(x) strcmp(x.name, '.'), all_recordings) | arrayfun(@(x) strcmp(x.name, '..'), all_recordings)));

            names                   = [vertcat(all_recordings.folder), repmat('/',numel(all_recordings), 1), vertcat(all_recordings.name)];
            extracted_data_paths    = cellfun(@(x) parse_paths(x), cellstr(names)', 'UniformOutput', false);
        end

        function update(obj, bypass, keep_2D, varargin)
            %% Update the arboreal scans
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.update(bypass, keep_2D, varargin)
            % -------------------------------------------------------------
            % Inputs:
            %   bypass (BOOL) - Optional - default is False
            %       If true, update is done without a prompt (this is used
            %       by the constructor for example when you build the obj)
            %   keep_2D (BOOL) - Optional, default is False
            %       If true, the full_data field is kept for each ROI. Not
            %       that this will considerably increase the object size.
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %   * Note that all analyses will be cleared. Do this only if
            %   you added / removed some recordings, or fundamentally
            %   changed something in the extracted arboreal_scans. 
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            if nargin < 2 || isempty(keep_2D)                
                keep_2D = false;
            end
            if nargin < 3 || isempty(bypass) || ~bypass
                quest = questdlg('WARNING : UPDATING SOURCES WILL DELETE ALL PROCESSED DATA. Continue?','Update?','Yes','No','No');
            else
                quest = 'Yes';
            end

            if strcmp(quest, 'Yes')
                obj.extracted_data_paths= list_sources(obj);
                if isempty(obj.extracted_data_paths) % extract from RAW
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

                        err = meta_batch_process_ribbon_scan(obj.source_folder,optional_settings_path, optional_analysis_params, fold);% QQ PLEASE CHECK IF THIS STILL WORKS
                        if isempty(err{1})
                            obj.source_folder       = fold;
                            obj.update(true,keep_2D);
                        else
                            error('Error detected during extraction. Check that the settings.txt file is present in the top_folder or manually indicated, and that it contains the correct paths')
                        end
                    else
                        return
                    end
                end
                obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update
                
                %% Clear all fields
                obj.reset();
                for field = {'arboreal_scans','binned_data','rescaling_info','event','variability','dimensionality'}
                    obj.(field{1}) = {};
                end
                
                %% Rebuild from extracted arboreal_scans
                for el = fliplr(find(obj.need_update))
                    add_tree(el, keep_2D);
                end
            end
            
            function add_tree(pos, keep_2D)
                %% INTERNAL FUNCTION THAT LOADS THE ARBOREAL_SCAN OBJECT
                obj.arboreal_scans{pos}                             = load(obj.extracted_data_paths{pos});
                if isa(obj.arboreal_scans{pos}.obj, 'arboreal_scan')
                    obj.arboreal_scans{pos}                         = obj.arboreal_scans{pos}.obj;
                    if ~keep_2D
                        obj.arboreal_scans{pos}.full_data           = []; % clear full data to save space.
                    end
                    obj.arboreal_scans{pos}.population_data         = []; % clear population_data.
                    obj.need_update(pos)                            = false;
                else
                    obj.arboreal_scans(pos)                         = [];
                    obj.need_update(pos)                            = [];
                end
            end
        end
        
        function reset(obj)
            %% Reset all analyzed fields
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.reset()
            % -------------------------------------------------------------
            % Inputs:
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
            
            %% Saving options
            obj.demo                = 0;
            obj.auto_save_analysis  = false;
            obj.auto_save_figures   = false;

            %% Analysis/extraction settings
            obj.filter_win      = [0, 0];
            obj.filter_type     = 'gaussian';
            obj.dim_red_type    = 'pca';
            obj.peak_thr        = 2;
            obj.bad_ROI_thr     = 0;
            obj.cc_mode         = 'groups_peaks'; % or raw
            obj.detrend         = false;
            obj.is_rescaled     = false; % is set to True once you ran the rescaling step.
            obj.rescaling_method= 'by_trials';
            obj.breakpoints     = []; % if you had a disruptive event during the experiment, a first scaling is done with large blocks

            %% All the fields computed in
            obj.binned_data     = [];	% Defines how ROIs are grouped
            obj.rescaling_info  = [];   % Defines how each ROi get rescaled to match the cell median
            obj.event           = [];   % Event detection output (event times, amplitude, correlation etc....)
            obj.variability     = [];   % Signal variability over time
            obj.dimensionality  = [];   % Results of diemensionality reduction
            obj.behaviours      = [];   % List of available behaviours
            obj.spiketrains     = [];   % If available, spike inference results
            obj.bad_ROI_list    = [];   % list of uncorrelated ROIs (following event detection)
            obj.updatable       = [];   % If arboreal_scan are still available, you could update the arboreal_scan_experiment compression
            obj.crosscorr       = [];   % Correlation of peaks/signal across ROIs/groups during/between bAps/activity_bouts
        end
        
        %% ########### GET METHODS ###################

        function updatable = get.updatable(obj)
            %% Get method for the updatable field
            % -------------------------------------------------------------
            % Syntax:
            %   updatable = obj.updatable
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   updatable (1xN BOOL)
            %   List of source folders that exist, and can be updated if
            %   required
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            updatable = cellfun(@(x) ~isfolder(x), obj.updated_data_path); % which is either the extracted data path or the updated one
        end

        function updated_data_path = get.updated_data_path(obj)
            %% Get method for the updatable field
            % -------------------------------------------------------------
            % Syntax:
            %   updated_data_path = obj.updated_data_path
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   updatable (1xN CELL ARRAY OF CHAR)
            %   List of update folders if obj.update_folder was provided,
            %   or the original extracted_data_paths if not. 
            % -------------------------------------------------------------
            % Extra Notes:
            %   * This doesn't tell you if the folder is valid. Use
            %   obj.updatable for that.
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            if isempty(obj.update_folder)
                updated_data_path = obj.extracted_data_paths;
            else
                updated_data_path = cellfun(@(x) get_fname(x), obj.extracted_data_paths, 'uni', false);
            end

            function fname = get_fname(in)
                %% INTERNAL FUNCTION
                [~, in_name, in_ext] = fileparts(in);
                fname = [in_name,in_ext];
                replacement = dir([obj.update_folder,'/**/',fname]);
                if ~isempty(replacement)
                    fname = [replacement(1).folder,'/',replacement(1).name];
                end
            end
        end

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
                breakpoints = obj.batch_params.breakpoints; %find(cellfun(@(x) contains(x, obj.batch_params.breakpoints),obj.updated_data_path));
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
            %   breakpoints (1xN CELL ARRAY OF CHAR)
            %   Numerical list of experiment number that interrupted the
            %   experiment. for example obj.breakpoints = [2,6,9];
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
                error('Cutoff must be between 0 and 1')
            else
                obj.bad_ROI_thr = value;
                a = dbstack(); 
                if ~contains([a.name], 'arboreal_scan_experiment.reset') % not useful when resetting
                    obj.find_bad_ROIs();
                end
            end
        end
        
        function extracted_traces = get.extracted_traces(obj)
            %% Get extracted traces from individual arboreal_scans
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_traces = obj.extracted_traces;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_traces (1xN CELL ARRAY of TxROI SINGLE)
            %   Original traces (detrended and smoothed if required) from
            %   the tree
            % -------------------------------------------------------------
            % Extra Notes:
            %  * If obj.detrend > 0, detrending is performed here
            %  * If filter_win  > 0, time smoothing is done here
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            extracted_traces = cellfun(@(x) x.simple_data, obj.arboreal_scans, 'UniformOutput', false);
            
            %% If expe was interrupted signal gain changed, we fix it here
            if ~isempty(obj.breakpoints) || obj.detrend
            	extracted_traces = obj.fix_changes_in_gain(extracted_traces);
            end

            %% Time smoothing if required
            if any(obj.filter_win)
                extracted_traces = cellfun(@(x) smoothdata(x, 'gaussian', obj.filter_win), extracted_traces, 'UniformOutput', false);
            end
        end

        function extracted_pop = get.extracted_pop(obj) % checked
            %% Get population signal from individual arboreal_scans
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_pop = obj.extracted_pop;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_pop (1xN CELL ARRAY of TxROI SINGLE)
            %   Original traces (detrended and smoothed if required) from
            %   the population recording
            % -------------------------------------------------------------
            % Extra Notes:
            %  * If filter_win  > 0, time smoothing is done here
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            extracted_pop = cellfun(@(x) x.simple_pop_data, obj.arboreal_scans, 'UniformOutput', false);
            if obj.detrend
                extracted_pop = cellfun(@(x) x - prctile(x, 1), extracted_pop, 'UniformOutput', false);
            end
        end

        function extracted_traces_conc = get.extracted_traces_conc(obj)
            %% All traces from the tree concatenated
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_traces_conc = obj.extracted_traces_conc;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_traces_conc (TxROI SINGLE)
            %   Original traces (detrended and smoothed if required) from
            %   the tree recording, concatenated
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            extracted_traces_conc = vertcat(obj.extracted_traces{:});
            %figure(666);cla();plot(normalize_sig(smoothdata(extracted_traces_conc,'gaussian',obj.filter_win)', '', 'norm_method','dF/F0','percentile',10)')
        end

        function extracted_pop_conc = get.extracted_pop_conc(obj)
            %% All traces from the population concatenated
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_pop_conc = obj.extracted_pop_conc;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_pop_conc (TxROI SINGLE)
            %   Original traces (detrended and smoothed if required) from
            %   the population recording, concatenated
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            extracted_pop_conc = vertcat(obj.extracted_pop{:});
            % figure(666);cla();plot(normalize_sig(smoothdata(extracted_pop_conc,'gaussian',obj.filter_win)', '', 'norm_method','dF/F0','percentile',10)')
        end

        function timescale = get.timescale(obj)
            %% Return timescale strcuture, with loads of timing info
            % -------------------------------------------------------------
            % Syntax:
            %   timescale = obj.timescale;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   timescale (STRUCT)
            %   Timescale information with the following fields:
            %       * sr : sampling rate per recording
            %       * time_source : source of timescale. Most accurate is
            %           "encoder" any in ({'command','computed','encoder'})
            %       * tp : timepoints per recording  
            %       * durations : duration per recording in seconds, based
            %           on time_source estimate
            %       * rec_timescale : 1xT timescale for each recording
            %       * global_timescale : 1xT timescale, for all recordings
            %           concatenated, ignoring any gap between recordings.
            %       * t_start_nogap : t start for each recording
            %       * datetime_start : real start time of each recording
            %       * real_timescale : 1xT timescale, for all recordings
            %           concatenated, including gaps between recordings.
            %       * t_start_real : t start for each recording, with start 
            %           of experiment at 0 
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            
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
            
            t_start_real_posix          = cellfun(@(x) posixtime(x.header.recording_t_start), obj.arboreal_scans);
            timescale.t_start_real      = t_start_real_posix - t_start_real_posix(1);
            timescale.datetime_start    = cellfun(@(x) x.header.recording_t_start, obj.arboreal_scans, 'UniformOutput', false);
            
            %timescale.real_timescale   = cellfun(@(x) diff(x), timescale.rec_timescale, 'UniformOutput', false);
            timescale.real_timescale   = cellfun(@(x, y) [x+y], timescale.rec_timescale, num2cell(timescale.t_start_real), 'UniformOutput', false);
            timescale.real_timescale   = horzcat(timescale.real_timescale{:});
        end

        function t = get.t(obj)
            %% Quick handle for global_timescale
            % -------------------------------------------------------------
            % Syntax:
            %   timescale = obj.timescale;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   t (1xT DOUBLE)
            %       equivalent to obj.timescale.global_timescale. Use it to 
            %       simplify your code
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            t = obj.timescale.global_timescale;
        end

        function ref = get.ref(obj)
            %% Quick handle for the first arboreal_scan
            % -------------------------------------------------------------
            % Syntax:
            %   ref = obj.ref;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   ref (arboreal_scan object)
            %       equivalent to obj.arboreal_Scan{1}. Use it to simplify
            %       your code
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            ref = obj.arboreal_scans{1};
        end

        function batch_params = get.batch_params(obj) 
            %% Quick handle for the batch_paramns of the first arboreal_scan
            % -------------------------------------------------------------
            % Syntax:
            %   batch_params = obj.batch_params;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   batch_params (STRUCT)
            %       equivalent to obj.arboreal_Scan{1}.batch_params. Use it
            %       to simplify your code
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            batch_params = obj.ref.batch_params;
            if ~isfield(batch_params, 'breakpoints')
                batch_params.breakpoints = [];
                warning('breakpoint field added post hoc. plase re-extract')
            end
        end

        function n_ROIs = get.n_ROIs(obj)
            %% Total Number of ROIs in the tree
            % -------------------------------------------------------------
            % Syntax:
            %   n_ROIs = obj.n_ROIs;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   n_ROIs (INT)
            %       Number of ROIs in the tree alone
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            n_ROIs = size(obj.ref.simple_data,2);
        end

        function n_pop_ROIs = get.n_pop_ROIs(obj)
            %% Total Number of ROIs in the population
            % -------------------------------------------------------------
            % Syntax:
            %   n_pop_ROIs = obj.n_pop_ROIs;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   n_pop_ROIs (INT)
            %       Number of ROIs in the population recording
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            n_pop_ROIs = size(obj.ref.population_data,2);
        end
        
        function logs = get.logs(obj)
            %% Individual arboreal scan logs (the comments you wrote
            % during the recording !)
            % -------------------------------------------------------------
            % Syntax:
            %   logs = obj.logs;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   logs (1xN CELL ARRAY OF CHAR)
            %       Number of ROIs in the population recording
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            logs = cellfun(@(x) x.log, obj.arboreal_scans, 'Uniformoutput', false); % vertcat(obj.logs{:})
        end

        function f_handle = get.default_handle(obj)
            %% Default handle used when you click on a tree
            % -------------------------------------------------------------
            % Syntax:
            %   f_handle = obj.f_handle;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   f_handle (FUNCTION HANDLE)
            %       This determine what happens when you click on one of
            %       the trees after doing all the meta analysis. Default
            %       handle enable you to reload the original Ribbon data
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            %
            % See also : load_several_experiments
            
            use_mask = false;
            f_handle = @(x) load_several_experiments(x, cellfun(@(x) x.data_folder, obj.arboreal_scans, 'UniformOutput', false), use_mask);
        end

        function global_median_raw = get.global_median_raw(obj)
            %% Median trace of all non-excluded RAW data
            % -------------------------------------------------------------
            % Syntax:
            %   global_median_raw = obj.global_median_raw;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   global_median_raw (Tx1 SINGLE)
            %       median RAW signal
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
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
            %% Return all external variables from original arboreal scans
            % -------------------------------------------------------------
            % Syntax:
            %   external_variables = obj.external_variables;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   external_variables (1xN STRUCT of external variables)
            %       Each external variable is a structure with a "value" 
            %       and a "time" field. Value can have more than one column
            %       (eg : X,Y,Z values for 3dMC, multiple ROIs for Motion
            %       index...)
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            %failed_encoder = cellfun(@(x) ~numel(x.analysis_params.external_var.encoder.time), obj.arboreal_scans);
            external_variables = cellfun(@(x) x.analysis_params.external_var, obj.arboreal_scans, 'UniformOutput', false);
        end

        function bad_ROI_list = get.bad_ROI_list(obj)
            bad_ROI_list = obj.bad_ROI_list;
            if ischar(bad_ROI_list) && strcmpi(bad_ROI_list, 'unset')
                warning('Uncorrelated ROIs were never set. Detecting now. To ignore this set value to []')
                %rendering = obj.rendering;obj.rendering = false;
                obj.find_events(); %find uncorrelated ROIs based on correlation
                %obj.rendering = rendering;
                bad_ROI_list = obj.bad_ROI_list;
            end
        end
        

        %% ########### ... ###################
        
        function extracted_traces = fix_changes_in_gain(obj, extracted_traces)
            %% Rescale gain changes between pairs of breakpoints
            % by correcting changes in F0 if obj.detrend > 0
            % -------------------------------------------------------------
            % Syntax:
            %   extracted_traces = EXPE.fix_changes_in_gain(extracted_traces);
            % -------------------------------------------------------------
            % Inputs:
            %   extracted_traces (1xN CELL ARRAY of TxROI SINGLE)
            %   Original traces, for each N recording extracted from each
            %   arboreal_scan.simple_data field
            % -------------------------------------------------------------
            % Outputs:
            %   extracted_traces (1xN CELL ARRAY of TxROI SINGLE)
            %   Original traces rescaled between each breakpoints to
            %   minimize change of baseline signal.
            % -------------------------------------------------------------
            % Extra Notes:
            %  * see fix_gain_changes() for methodology. Briefly, signal is
            %  rescaled for each ROI, between each pair of breakpoints to
            %  correct for linear change of F0 (eg from water drying out)
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            extracted_traces = fix_gain_changes(obj, extracted_traces);
        end

        %% #############################

        function [tree, soma_location, tree_values, values] = plot_distance_tree(obj, bin)
            if nargin < 2 || isempty(bin)
                bin = [];
            end
            [tree, soma_location, tree_values, values] = obj.ref.plot_dist_tree(bin);
        end

        %         function [tree, soma_location, tree_values, values] = plot_seg_length_tree(obj)
        %             %% Map dimension weights on the tree
        %             [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(Pvec_tree(obj.ref.simplified_tree{1}), '', '', 'Segment length');
        %         end
        

        function animate_experiment(obj, values, timepoints, ROIs, color_range)
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
            elseif numel(smoothing) == 1
                smoothing = [smoothing, 0];
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
            %% Set cross correlation mode
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.get_correlations(cc_mode)
            % -------------------------------------------------------------
            % Inputs:
            %   cc_mode (STR) - See Description for details - Optional - 
            %           Default is your current obj.cc_mode value
            %       update EXPE.cc_mode with existing value
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %   * Correlation is computed using:
            %       - A reference trace (located at index 1)
            %       - The data you want to correlate (individual ROI, 
            %           OR binned data, see below)
            %       - Optionally, you can append population data
            %   * cc_mode is build by concatenating strings defining the
            %   timepoints to use (all of them or a subset), the traces to  
            %   use (all of them or the binned data), whether you want to 
            %   use population data too and whether you want to use some 
            %   behavioural metrics. for example : 
            %      - 'groups_pop' : correlation between binned median
            %        traces and population data using all timepoints 
            %      - 'active' : correlation between all ROIs, using all
            %      timepoints when cell is active (i.e. during bAPs)
            %      - 'encoder_peaks' : correlation between all ROIS, only
            %        at peak time, and using encoder as ref  
            %      - 'encoder_peaks_refilter' : correlation between all 
            %        ROIS, only at peak time, AND only when running speed
            %        is high    
            %      - '~encoder_peaks_refilter' : correlation between all 
            %        ROIS, only at peak time, AND only when running speed
            %        is low  
            %   * The reference (located at index 1 of the matrix) is
            %    either 
            %       - The signal data during active behaviour IF you pass a
            %       valid behaviour in CC mode (list available when typing
            %       obj.behaviours.types)
            %       - The averaged somatic data when available (as defined
            %         by ref.indices.somatic_ROIs), or the nexus signal 
            %         when somatic ROIS are not available
            %   * Correlation matrix can be generated using the entire 
            %     trace (default) or a subset of timepoints : 
            %       - If cc_mode contains 'peaks', only the values at peak
            %       times (as defined by event.fitting.peak_pos) are used
            %       - If cc_mode contains 'active', all the timepoints
            %       containing bAPs (as defined by event.fitting.t_win)
            %       - If cc_mode contains 'quiet', all the timepoints
            %       between bAPs (as defined by ~event.fitting.t_win)
            %   * Correlation matrix can be generated using
            %       - every ROIs (except for excluded ROIs, indicated in
            %     ref.indices.valid_swc_rois)
            %       - If cc_mode contains 'groups', using the median traces 
            %       after categorical binning  (as defined by
            %        binned_data.median_traces). Use the 
            %   * To include population data, include 'pop' in cc_mode
            %   * If a behaviour type is added, the reference won't be the
            %   somatic average but the indicated behaviour for the same
            %   timepoints
            %   * If cc_mode contains refilter AND a valid behaviour, the
            %   timepoints used will be refiltered to include only events 
            %   during active behaviour (as defined by 
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022
            
        	if ~strcmp(obj.cc_mode, cc_mode)
                obj.crosscorr = []; %if you change the cc mode, clear the previous correlation results
            end
            
            msg = '\n';
            if contains(cc_mode, 'peaks')
                msg = [msg, 'Correlation done using all peaks, '];
            elseif contains(cc_mode, 'active')
                msg = [msg, 'Correlation done using timepoints during active behaviour, '];
            elseif contains(cc_mode, 'quiet')
                msg = [msg, 'Correlation done using timepoints during quiet behaviour, '];             
            else
                msg = [msg, 'Correlation done using all timepoints, '];
            end            
            if contains(cc_mode, 'groups')
                msg = [msg, 'Correlation done by groups using your bins, '];
            else
                msg = [msg, 'Correlation done between all ROIs, '];
            end            
            if contains(obj.cc_mode, 'pop')
                msg = [msg, 'including population data'];
            end
            
            fprintf([msg ,  '\n'])
            obj.cc_mode = cc_mode;
        end        
        
        function cross_corr = get_correlations(obj, cc_mode)
            %% Returns cross correlation (and triggers computation if required)
            % -------------------------------------------------------------
            % Syntax:
            %   cross_corr = EXPE.get_correlations(cc_mode)
            % -------------------------------------------------------------
            % Inputs:
            %   cc_mode (STR) - See Description for details - Optional - 
            %       Default is [];
            %       If provided, update EXPE.cc_mode. see set.cc_mode for
            %       more details
            % -------------------------------------------------------------
            % Outputs:
            %   crosscorr (NxN DOUBLE)
            %       Cross correlation between ROIs OR groups based on
            %       traces OR Peaks, depending on the settings. see
            %       set.cc_mode doc for the deails
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022            
            
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
            %% Compute cross correlation 
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.set_crosscorr()
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %       Cross correlation is computed between ROIs OR groups 
            %       based on traces OR Peaks, depending on the settings. 
            %       see set.cc_mode doc for the details
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022
            
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

            %% Get signal to use
            if contains(obj.cc_mode, 'groups')
                signal = obj.binned_data.median_traces;
            else% if contains(obj.cc_mode, 'ROIs')
                signal = obj.rescaled_traces(:,obj.ref.indices.valid_swc_rois);                
            end
            mode = erase(mode, {'ROIs','groups'});

            %% Add population signal if needed
            if contains(obj.cc_mode, 'pop')
                if isempty(obj.extracted_pop_conc)
                    warning('No population data for this recording')
                    pop = [];
                else
                    pop = obj.extracted_pop_conc;
                end
            else
                pop = [];
            end
            mode = erase(mode, {'pop','_', 'refilter', '~'});

            %% Get ref ROIs and trace
            beh = obj.behaviours.types(find(contains(obj.behaviours.types, mode)));
            if ~isempty(beh)
                [~, ~, beh] = obj.get_behaviours(beh);                
                if contains(obj.cc_mode, 'refilter')
                    [~,~, active_tp] = obj.get_activity_bout('encoder', true);
                    if contains(obj.cc_mode, '~')
                        active_tp = ~active_tp;
                    end
                    tp               = find(ismember(tp, find(active_tp)));
                end
                ref         = beh.value(tp)';
            else
                somatic_ROIs= obj.ref.indices.somatic_ROIs;
                ref         = nanmean(obj.rescaled_traces(:, somatic_ROIs),2);
            end

            %% Build the arrays used for the correlation matrix
            variable        = [ref, signal(tp, :)];
            if ~isempty(pop)
                variable = [variable, pop(tp,:)];
            end
            obj.crosscorr   = corrcoef(variable,'Rows','Pairwise')';
        end

        function crosscorr = get.crosscorr(obj)
            %% Returns cross correlation (and triggers computation if required)
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.crosscorr()
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   crosscorr (NxN DOUBLE)
            %       Cross correlation between ROIs OR groups based on
            %       traces OR Peaks, depending on the settings. see
            %       set.cc_mode doc for the details
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022
            
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
                error('You need to download the "bricks" and "ml_spikes" toolboxes, and add then to the path')
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
        
        function gui(obj)
            explore_factors(obj);
        end
        
        function weighted_averages = get_dimensionality(obj, cross_validate, n_factors, timepoints, dim_red_type, clust_meth, N_clust)
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
            if nargin < 4 || isempty(timepoints) || (ischar(timepoints) && contains(timepoints, 'full'))
                timepoints                  = true(size(obj.timescale.global_timescale));
                variable                      = 'full_traces'; 
            elseif ischar(timepoints) && strcmpi(timepoints, 'peaks')
            	timepoints     = vertcat(obj.event.peak_time{:});
                 variable       = 'peaks';
            elseif ischar(timepoints)
                [timepoints, beh_sm, active_tp] = obj.get_activity_bout(timepoints, true, 5);
                1
                timepoints = timepoints(2:2:end)
                variable       = timepoints;
            else
                variable                      = 'manual selection'; 
            end   
            if nargin < 5 || isempty(dim_red_type)
                dim_red_type                  = 'pca';
            end
            if nargin < 6 || isempty(clust_meth)
                clust_meth                    = 'none';
            end
            if nargin < 7 || isempty(N_clust)
                N_clust                    = 0;
            end
            
%             t = 1:numel(obj.global_median_raw);
%             t_bap = [obj.event.t_win_no_overlap{:}];
%             t_no_bap = ~ismember(t, t_bap);
% 
%             timepoints = find(t_no_bap);
            
            obj.dimensionality.dim_red_type     = dim_red_type;
            obj.dimensionality.variable         = variable;
            
            data                = obj.rescaled_traces(timepoints,:);
            
            %rescaled_traces     = rescaled_traces - obj.binned_data.global_median;
            
            all_ROIs            = 1:size(data, 2);
            normal_n_NaN        = median(sum(isnan(data(:,~all(isnan(data)))))) * 4; % get an indicative number of NaN in a normal traces, and set acceptable thr at 4 times that
            valid_trace_idx     = sum(isnan(data)) <= normal_n_NaN; % eclude traces with too many NaNs (eg. traces that got masked completely)
            data                = fillmissing(data(:, valid_trace_idx),'spline'); % removed funny traces
            
            %% Need to set it now
            obj.dimensionality.valid_trace_idx = valid_trace_idx;

            
            %data = data - nanmedian(data,2)
            data = data - nanmean(data,2);
            
            
            %% Get single or multiple factor estimate
            if ~cross_validate
                T = [];
                stats = {};
                specVarPM = [];
                
                switch dim_red_type
                    case 'pca'
                        %[LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(rescaled_traces),'NumComponents',obj.dimensionality.n_factors);
                        [LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(data),'NumComponents',obj.dimensionality.n_factors);

                        %[LoadingsPM,F,specVarPM,stats,explained,mu] = pca(double(data) - obj.binned_data.global_median(vertcat(obj.event.peak_time{:})),'NumComponents',obj.dimensionality.n_factors);

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
                
                %% Store results
                obj.dimensionality.LoadingsPM         = LoadingsPM;
                obj.dimensionality.specVarPM          = specVarPM;
                obj.dimensionality.T                  = T;                % Rotation matrix
                obj.dimensionality.stats              = stats;            % Factoran stats
                obj.dimensionality.F                  = F;                % components
                obj.dimensionality.all_ROIs           = all_ROIs;         % first occurences
                obj.dimensionality.mask               = timepoints;
                
                obj.dimensionality.cluster_idx        = [];
                obj.dimensionality.clust_meth         = clust_meth;
                obj.dimensionality.N_clust            = N_clust;

                %% Plot weight-tree for each component
                for comp = 1:obj.dimensionality.n_factors
                    obj.plot_dim_tree(comp);
                end
                
                if obj.rendering
                    obj.plot_factor_tree();      
                end

                %% Plot a map of tree weights by ROI number for each component
                weighted_averages = obj.get_weight_map();
                
                %% If no weighted average, Find clusters from latent variables
                if ~strcmp(clust_meth, 'none')
                    obj.cluster_factors();
                end

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
        
        function cluster_factors(obj)
            if ~obj.dimensionality.N_clust
                obj.dimensionality.cluster_idx = (1:size(obj.dimensionality.LoadingsPM, 1))';
                obj.dimensionality.sorted_idx = (1:size(obj.dimensionality.LoadingsPM, 1))';
                return
            end
            
            if isinf(obj.dimensionality.N_clust)
                R = 1:20;
                if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'kmeans','silhouette','KList',R);
                elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'linkage','silhouette','KList',R);
                else
                    eva = {};
                    eva.OptimalK = inf;
                    fprintf(['No auto auto-determinatiopn of the number of cluster\n'])
                end
                obj.dimensionality.N_clust = eva.OptimalK
                fprintf(['Ideal number of clusters is ',num2str(obj.dimensionality.N_clust), '\n'])
            end
            
            if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                %cluster_idx = kmeans_opt(obj.dimensionality.LoadingsPM, max_n, 0.95, 500);
                cluster_idx = kmeans(obj.dimensionality.LoadingsPM , obj.dimensionality.N_clust);
            elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical') % see https://fr.mathworks.com/help/stats/hierarchical-clustering.html
                cluster_idx = clusterdata(obj.dimensionality.LoadingsPM,'Linkage', 'ward', 'MAXCLUST', obj.dimensionality.N_clust);%, 'Criterion','distance' 'MAXCLUST', 40)
                figure(3);clf(); silhouette(obj.dimensionality.LoadingsPM,cluster_idx)
                %                 Y = pdist(obj.dimensionality.LoadingsPM ,'euclidean');
                %                 Z = linkage(Y,'ward');
                %                 figure();dendrogram(Z);
                %                 W = inconsistent(Z,3)

                %                 pause(0.1); drawnow
                %                 Y = pdist(obj.dimensionality.LoadingsPM ,'cityblock');
                %                 Z = linkage(Y,'average');
                %                 %c = cophenet(Z,Y); % cophenetic correlation coefficient
                %                 %I = inconsistent(Z)
                %                 %T = cluster(Z,'cutoff',obj.dimensionality.N_clust)
                %                 T = cluster(Z,'maxclust',10);
                %                 cluster_idx = T;
            elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                MIN_GP = 3   
%                 
%                 [mIdx,mD] = knnsearch(obj.dimensionality.LoadingsPM,obj.dimensionality.LoadingsPM,'K',2);
%                 dist = nanmean(mD(:,2))
%                 cluster_idx = dbscan(obj.dimensionality.LoadingsPM ,dist*2,MIN_GP, 'Distance', 'euclidean');
                %if isinf(obj.dimensionality.N_clust)                    
%                     mean_v = [];
%                     L = 0.01:0.01:1
%                     [mIdx,mD] = knnsearch(obj.dimensionality.LoadingsPM,obj.dimensionality.LoadingsPM,'K',2);
%                     nanmean(mD(:,2))
%                     
%                     for v = L
%                         %cluster_idx = clusterdata(obj.dimensionality.LoadingsPM, 'Linkage', 'ward', 'MAXCLUST', v);%, 'Criterion','distance' 'MAXCLUST', 40)
%                         cluster_idx = dbscan(obj.dimensionality.LoadingsPM ,v,MIN_GP, 'Distance', 'euclidean');
%                         %labels = labels - min(labels)+ 1
%                         %figure();silhouette(obj.dimensionality.LoadingsPM,cluster_idx);
%                         a = silhouette(obj.dimensionality.LoadingsPM, cluster_idx);
%                         mean_v(end+1) = nanmean(a);
%                     end
%                     [~, best] = nanmax(mean_v);
%                     obj.dimensionality.N_clust = L(best);
%                     cluster_idx = dbscan(obj.dimensionality.LoadingsPM ,obj.dimensionality.N_clust,MIN_GP, 'Distance', 'euclidean');
%                 %end
%                 


                
                kD = pdist2(obj.dimensionality.LoadingsPM ,obj.dimensionality.LoadingsPM ,'euc','Smallest',5)
                figure();
                plot(sort(kD(end,:)));
                title('k-distance graph')
                xlabel('Points sorted with 50th nearest distances')
                ylabel('50th nearest distances')
                grid
                thr = 1;
                
                M = [];
                L = 0.01:0.01:1
                for thr = L
                    thr = thr - thr/20;                           
                    labels = dbscan(obj.dimensionality.LoadingsPM ,thr, MIN_GP, 'Distance', 'euclidean');
                    M(end + 1) = numel(unique(labels))
                end
                [~, maxloc] = max(M);
                L(maxloc)
                labels = dbscan(obj.dimensionality.LoadingsPM ,L(maxloc),MIN_GP, 'Distance', 'euclidean');
                labels = labels - min(labels)+ 1
                col = lines(numel(unique(labels)));
                labels_col = [labels, labels, labels];
                count = 1;
                for v = unique(labels)'
                    subset = labels_col(:,1) == v;                            
                    labels_col(subset, :) = repmat(col(count,:),sum(subset),1);
                    count = count + 1;
                end
                figure();
                scatter3(obj.dimensionality.LoadingsPM (:,1),obj.dimensionality.LoadingsPM (:,2),obj.dimensionality.LoadingsPM (:,3),20,labels_col)
                %gscatter(obj.dimensionality.LoadingsPM (:,2),obj.dimensionality.LoadingsPM (:,3),labels);
                title('epsilon = 2 and minpts = 50')
                grid
                cluster_idx = labels;
            else
                obj.dimensionality.cluster_idx = [];
            end
            
            obj.dimensionality.cluster_idx = ones(size(obj.dimensionality.cluster_idx));
            obj.dimensionality.sorted_idx   = 1:(numel(obj.dimensionality.cluster_idx));
            
            if ~strcmp(obj.dimensionality.clust_meth, 'none')
                %% Sort clusters by number of elements
                [~, gp] = sort(hist(cluster_idx,unique(cluster_idx)), 'descend');
                a = unique(cluster_idx)';
                gp = a(gp);
                idx_sorted = NaN(size(cluster_idx));
                count = 1;
                for gp_idx = gp  
                    idx_sorted(cluster_idx == gp_idx) = count;
                    count = count + 1;
                end
                obj.dimensionality.cluster_idx = idx_sorted; % stored group ids (reordered by number of element)
                [~, obj.dimensionality.sorted_idx] = sort(idx_sorted); % get ROI index to reorder the data by group

                if obj.rendering
                    %% Plot clusters
                    obj.plot_cluster_tree();
                end
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
                    temp = load(obj.updated_data_path{rec});
                    obj.arboreal_scans{rec}.simple_data = temp.obj.simple_data;
                    clear temp;
                end

                obj.arboreal_scans{rec}.update_segment_signals(new_method);
                obj.arboreal_scans{rec}.simple_data = [];
                obj.need_update(rec)                    = true;
            end
            error_box('COMPRESSION MODES WERE UPDATED BUT YOU NEED TO SAVE THE ARBOREAL SCANS TO KEEP THIS CHANGE FOR NEXT RELOADING')
        end

        function find_events(obj, idx_filter, thr_for_detection, method)
            if nargin < 2 || isempty(idx_filter)
                idx_filter = 1:obj.n_ROIs;
            elseif ischar(idx_filter) && strcmp(idx_filter, 'soma')
                idx_filter = obj.ref.indices.somatic_ROIs;
            end
            if nargin < 3 || isempty(thr_for_detection)
                thr_for_detection = 0.2;
            end
            if nargin < 4 || isempty(method)
                method = 'corr';
            end
            
            fprintf(['Now detecting events with global pairwise correlation at least > ',num2str(thr_for_detection),' %% \n'])     

            %% Get original traces (unscaled)
            raw_traces              = obj.extracted_traces_conc(:, idx_filter);

            %% Get original traces
            [obj.event, correlation_res] = detect_events(raw_traces, obj.t, method, thr_for_detection, [], obj.rendering);
            %[obj.event]             = detect_events(raw_traces, obj.t, 'global_amp', 60, corr_window);

            %% Identify and log poorly correlated ROIs
            obj.find_bad_ROIs(correlation_res, idx_filter);
        end

        function [bad_ROIs, mean_corr_with_others] = find_bad_ROIs(obj, correlation_res, ROIs)
            if nargin < 3 || isempty(ROIs)
                ROIs = 1:obj.n_ROIs;
            end            
            if nargin < 2 || isempty(correlation_res)
                warning('BAD ROI IDENTIFICATION RELIES ON EVENT DETECTION. EVENT DETECTION FIELD SEEMS EMPTY. RUNNING DEFAULT EVENT DETECTION')
                raw_data = obj.extracted_traces_conc;
                [events, correlation_res] = detect_events(raw_data(:, ROIs), obj.t, 'corr', 0.2, [], false, []); % no rendering here
            else
                events = obj.event;
            end

            
            THR_FOR_CONNECTION      = 0.2 % Defines what level of minimal pairwise correlation means "these two ROIs are connected"
            fprintf(['* Now detecting ROIs that are either,\n'...
                     '      - so poorly correlated to the rest of the tree that they probably belong to another cell (or have no signal).\n',...
                     '      - are member of batch_params.excluded_branches \n',...
                     '* Threshold for exclusion is  ',num2str(THR_FOR_CONNECTION),' %% \n'])  
            
            %% Show mean correlation with each ROI
            mean_corr_with_others   = [];
            max_corr_for_cell       = max(events.globality_index(2:end)); % QQ 1st point sometimes show some artifacts
            for key = 1:numel(ROIs)
                corr_results_sub    = correlation_res.corr_results(events.t_corr(events.is_global), correlation_res.comb(:,2) == key | correlation_res.comb(:,1) == key);
                corr_results_sub    = [corr_results_sub(:,1:(key-1)), NaN(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];
                mean_corr           = nanmean(corr_results_sub,1);
                mean_corr_with_others(key) = sum(mean_corr > THR_FOR_CONNECTION);
            end

            %% Normalize to 100%
            mean_corr_with_others = mean_corr_with_others / numel(mean_corr_with_others); % renormalize to max possible corr for this cell
            
            if obj.rendering
                figure(88888);cla();hist(100*mean_corr_with_others,0:2:100); hold on; 
                xlabel('%% of correlation with all other ROIs (from the tree)'); ylabel('counts')
            end
            
            %% Renormalize to max possible corr for this cell
            mean_corr_with_others_norm = mean_corr_with_others / max(mean_corr_with_others); 
            if obj.rendering
                figure(88889);cla();hist(100*mean_corr_with_others_norm,0:2:100); hold on; 
                xlabel('%% of correlation with all other ROIs (from the tree)'); ylabel('counts')
            end
            
            %% Update class variables
            if obj.bad_ROI_thr ~= THR_FOR_CONNECTION
                obj.bad_ROI_thr = THR_FOR_CONNECTION;
            end            
            bad_ROIs = mean_corr_with_others_norm < obj.bad_ROI_thr;            
            obj.bad_ROI_list = bad_ROIs; % below threshold% of max correlation
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
                
                %% Plot (if possible) the excluded traces
                try
                    regroup_traces(bad', 80, 'pca') 
                end
                
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
            obj.find_events();

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
            try
                obj.get_dimensionality(false,6,'','factoran',false,'peaks'); % set to true for cross validation
            catch
                warning('Not enough events availabl for factoran on peaks. Using full traces instead')
                try
                    obj.get_dimensionality(false,6,'','pca',false,'peaks'); % set to true for cross validation 
                catch
                    warning('Not enough events available for pca on peaks. Using full traces instead')
                    obj.get_dimensionality(false,6,'','factoran',false,'full_traces'); % set to true for cross validation 
                end
            end

            %% Optionally, if external variables need an update, do it here
            %obj.update_external_metrics(60)
            
            %% Extract behaviours
            obj.get_behaviours();

            %% Spike inference
            try
            %    obj.get_spike_trains();
            end

            %% Save figure and/or analysis
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
