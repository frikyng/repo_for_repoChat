%% Note, go to the process() method for an overview of the available features


%% TO DO :
% - improve update system when changing folder.
% - enable loading if we just have a folder with arboreal_scans objects
% - add warning if source folder is the actual raw experiment
% - add doc to load_Several_Experiment
% - add use_fitting_corrected variable flag, to replace all isfield(obj.event, 'fitting')

%% TO CHECK:
% - loading from raw data
% - Breakpoints not iintorucing artefacts
% - Check if bheviours are fine
% why are camera behaviours downsampled. Do we want that?

%% Suclasses
% signal_manipulation contains tools to rescale and detrend signals, and find uncorrelated ROIs
% arboreal_scan_plotting contains all the rendering tools
% event_detection contains tools to detect and fit events
% behaviours_analysis contains tools to detect and extract behaviour activity bouts
% correlation_analysis contains tools to generate correlation matrices using specific sets of variables
% cluster_analysis can be used to cluster dimensionality reduction results OR correlation matrices

classdef arboreal_scan_experiment < handle & arboreal_scan_plotting & event_detection & behaviours_analysis & correlation_analysis & cluster_analysis & signal_manipulation
    properties
        %% Extraction/Re-extraction Settings
        extraction_method       = 'median';     % the 1D -> 0D compression method. If None, XxT data is kept, but this make the files heavier
        source_folder                           % The folder where individual arboreal_scans were located when you built the object
        update_folder           = '';           % If you need to update the arboral_scans, but the orginal path changed, set the folder containing the new files here
        extracted_data_paths                    % The original location of the individual arboreal_scans when you built the object
        need_update                             % A series of obj.n_rec boolean indicating whether some recordings need to be reextracted

        %% Children
        arboreal_scans                          % A copy of the individual arboreal_scans, but where the uncompressed signal per ROIs had ben deleted. See keep_2D for details

        %% Across expe Timescale info
        timescale                               % time and sampling rate info for the individual recordings

        %% Saving options
        demo                    = 0;            % Defines the level of display info across the different functions
        auto_save_analysis      = false         % If true, the arboreal_scan_experiment object is saved automatically after calling obj.process
        auto_save_figures       = false         % If true, analysis figures from obj.process() are saved automatically

        default_handle
        variability_metric      = 'vmr'         % Defines the mretic that will be used to comput signal or events variability (vmr or cv)        

        %% All the fields computed in
        binned_data                             % Defines how ROIs are grouped
        rescaling_info                          % Defines how each ROi get rescaled to match the cell median
        event                                   % Event detection output (event times, amplitude, correlation etc....)
        variability                             % Signal variability over time
        dimensionality                          % Results of diemensionality reduction
        spiketrains                             % If available, spike inference results
        updatable                               % If arboreal_scan are still available, you could update the arboreal_scan_experiment compression
    end

    properties (Dependent = true, Transient = true)
        updated_data_path                       % If update_folder is used, the updated filpath
        extracted_traces                        % Concatenated version of each obj.arboral_scan.simple_data
        extracted_pop                           % Concatenated version of each obj.arboral_scan.simple_pop_data
        extracted_traces_conc                   % Concatenated version of extracted_traces
        extracted_pop_conc                      % Concatenated version of extracted_pop
        global_median_raw                       % The median of extracted_traces_conc
        global_median_rescaled                  % The median of rescaled traces
        t                                       % Pointer to obj.timescale.global_timescale
        tp                                      % Imaging timpoints
        n_ROIs                                  % Total number of ROIs in the swc, including bad ones
        n_pop_ROIs                              % Total number of population ROIs
        n_rec                                   % Total number of recordings
        ref                                     % A pointer to the first extracted arboreal scan, for conveniency
        batch_params                            % Pointer to obj.ref.batch_params ; the info to rebuild and locate the tree
        logs                                    % Display the logs for all experiments
        usable_swc_ROIs                         % Rois that are neither bad nor excluded
    end

    properties (Dependent = true, Transient = true, Hidden = true)
        external_variables                      % Pointer to behavioural variables of each arboreal scan --> set in obj.behaviours
        
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
            %   source_folder (STR) - Optional - Default is current folder
            %       Path to the raw or extracted recordings. It can be :
            %           * a folder containing extracted arboreal_scans
            %             objects (recommended). If an axisting extracted
            %             arboreal_scan_experiment is available, the code
            %             will as if you want to load it.
            %           * an exisiting arboreal_scan_experiment (.mat file)
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

            if nargin < 1 || (~ischar(source_folder) && isnan(source_folder))
                return % empty object when building whole dataset                
            elseif isempty(source_folder)
                source_folder = pwd;
            end
            if nargin < 2 || isempty(keep_2D)
                keep_2D = false;
            end

            %% Fix paths
            source_folder       = parse_paths(source_folder);
            
            %% Check if there is already an aorboreal_scan_experiment. 
            if contains(source_folder, '.mat')
                obj             = importdata(source_folder, 'obj'); %instead of load. This assign directly the object to obj.
                source_folder   = fileparts(source_folder);
                obj.source_folder = parse_paths(source_folder);  
                return
            else
                [~, tag] = fileparts(fileparts(source_folder));
                if isfile([source_folder, tag, '.mat'])
                    %answ = questdlg('Existing extracted file detected. Do you want to reload existing file or rebuild it?','','Reload','Rebuild','Abort','Reload');
                    answ = 'Reload'
                    if strcmp(answ, 'Abort')
                        return
                    elseif strcmp(answ, 'Reload')
                        source_folder   = [source_folder, tag, '.mat'];
                        obj             = importdata(source_folder, 'obj'); %instead of load. This assign directly the object to obj.
                        source_folder   = fileparts(source_folder);
                        obj.source_folder = parse_paths(source_folder);
                        return
                    end
                end
                obj.source_folder = parse_paths(source_folder);                
            end
   
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
                obj.disp_info(['Original arboreal_scans not found in : ',[obj.source_folder,'/**/*-*-*_exp_*_*-*-*'],' . extraction/re-extraction not available. If you moved the files to a new top folder, change obj.update_folder accordingly'],2);
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
                quest = questdlg('WARNING : UPDATING SOURCES WILL DELETE ALL PROCESSED DATA (unless no folder needs updating). Continue?','Update?','Yes','No','No');
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

                        err = meta_batch_extract_arboreal_scan(obj.source_folder,optional_settings_path, optional_analysis_params, fold);% QQ PLEASE CHECK IF THIS STILL WORKS
                        if isempty(err{1})
                            obj.source_folder       = fold;
                            obj.update(true,keep_2D);
                        else
                            obj.disp_info('Error detected during extraction. Check that the settings.txt file is present in the top_folder or manually indicated, and that it contains the correct paths',4);
                        end
                    else
                        return
                    end
                elseif ~isempty(obj.arboreal_scans)
                    %% Check if all existing arboreal_scans are included
                    lastFolders_existing_as = cellfun(@(x) strsplit(fileparts(x.data_folder),'/'), obj.arboreal_scans, 'UniformOutput', false);
                    lastFolders_existing_as = cellfun(@(x) x{end}, lastFolders_existing_as, 'UniformOutput', false);
                    
                    lastFolders_extracted = cellfun(@(x) strsplit(fileparts(x),'/'), obj.extracted_data_paths, 'UniformOutput', false);
                    lastFolders_extracted = cellfun(@(x) x{end}(end-7:end), lastFolders_extracted, 'UniformOutput', false);
                    
                    if numel(lastFolders_existing_as) ~= numel(lastFolders_extracted) || any(cellfun(@(x,y) any(x ~= y), lastFolders_extracted,lastFolders_existing_as))
                        % at least one mismatch. Will need update
                    else
                        return
                    end
                end
                obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update

                %% Clear all fields
                obj.reset();
                obj.arboreal_scans = {};

                %% Rebuild from extracted arboreal_scans
                for el = fliplr(find(obj.need_update))
                    add_tree(el, keep_2D);
                end
                obj.need_update(:) = false;
            end

            function add_tree(pos, keep_2D)
                %% INTERNAL FUNCTION THAT LOADS THE ARBOREAL_SCAN OBJECT. REMOVE UNCOMPRESSED DATA IF REQUIRED
                obj.arboreal_scans{pos}                             = load(obj.extracted_data_paths{pos});
                if isa(obj.arboreal_scans{pos}.obj, 'arboreal_scan')
                    obj.arboreal_scans{pos}                         = obj.arboreal_scans{pos}.obj;
                    if ~keep_2D
                        obj.arboreal_scans{pos}.full_data               = []; % clear full data to save space.
                        obj.arboreal_scans{pos}.population_data         = []; % clear population_data.
                    end                    
                    obj.need_update(pos)                            = false;
                else
                    obj.arboreal_scans(pos)                         = [];
                    obj.need_update(pos)                            = [];
                end
            end
        end

        function reset(obj, deep_reset)
            %% Reset all analyzed fields
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.reset()
            % -------------------------------------------------------------
            % Inputs:
            %   deep_reset(BOOL) - Optional - default is false
            %       if true, all the trace manipulation options are reset
            %       too
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

            if nargin < 2 || isempty(deep_reset)
                deep_reset = false;
            end
            
            if deep_reset
                obj.disp_info('ALL FIELDS AND OPTIONS HAVE BEEN RESET', 3)
            else
                obj.disp_info('Analysis Fields have been reset. Other settings (smoothing options, thresholds etc...) were conserved', 2)
            end
            
            %% Saving options
            obj.demo                = 0;
            obj.auto_save_analysis  = false;
            obj.auto_save_figures   = false;

            %% Analysis/extraction settings
            if deep_reset
                obj.time_smoothing      = [0, 0];
                obj.filter_type         = 'gaussian';
                obj.peak_thr            = 2;
                obj.bad_ROI_thr         = 0.2;
                obj.cc_mode             = 'groups_peaks'; % or raw
                obj.detrend             = false;
                obj.is_rescaled         = false; % is set to True once you ran the rescaling step.
                obj.rescaling_method    = 'by_trials_on_peaks';
                obj.breakpoints         = []; % if you had a disruptive event during the experiment, a first scaling is done with large blocks
            end
            
            %% All the fields computed in
            obj.binned_data         = [];	% Defines how ROIs are grouped
            obj.rescaling_info      = [];   % Defines how each ROi get rescaled to match the cell median
            obj.event               = [];   % Event detection output (event times, amplitude, correlation etc....)
            obj.variability         = [];   % Signal variability over time
            obj.dimensionality      = [];   % Results of diemensionality reduction
            obj.behaviours          = [];   % List of available behaviours
            obj.spiketrains         = [];   % If available, spike inference results
            obj.bad_ROI_list        = [];   % list of uncorrelated ROIs (following event detection)
            obj.updatable           = [];   % If arboreal_scan are still available, you could update the arboreal_scan_experiment compression
            obj.crosscorr           = [];   % Correlation of peaks/signal across ROIs/groups during/between bAps/activity_bouts
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
            %  * If time_smoothing  > 0, time smoothing is done here
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            if ~obj.use_hd_data
                extracted_traces = cellfun(@(x) x.simple_data, obj.arboreal_scans, 'UniformOutput', false);
            else
                extracted_traces = cellfun(@(x)  squeeze(cat(1, x.full_data{:})), obj.arboreal_scans, 'UniformOutput', false);
                extracted_traces = cellfun(@(x) x(:,:,1)', extracted_traces, 'UniformOutput', false);
            end

            %% If expe was interrupted signal gain changed, we fix it here
            if (~isempty(obj.breakpoints) || ~isempty(obj.detrend))% && ~obj.is_detrended
            	extracted_traces = obj.fix_changes_in_gain(extracted_traces);
                if isfield(obj.binned_data, 'median_traces') && ~obj.is_detrended
                    obj.disp_info('Changing the detrending method affects several preprocessing steps. Analyisis fields need to be reset', 2);
                    obj.reset();
                end                
            end

            %% Time smoothing if required
            if any(obj.time_smoothing)
                extracted_traces = cellfun(@(x) smoothdata(x, 'gaussian', obj.time_smoothing), extracted_traces, 'UniformOutput', false);
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
            %  * If time_smoothing  > 0, time smoothing is done here
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
            extracted_traces_conc(isinf(extracted_traces_conc))    = NaN;
            %figure(666);cla();plot(normalize_sig(smoothdata(extracted_traces_conc,'gaussian',obj.time_smoothing)', '', 'norm_method','dF/F0','percentile',10)')
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
            % figure(666);cla();plot(normalize_sig(smoothdata(extracted_pop_conc,'gaussian',obj.time_smoothing)', '', 'norm_method','dF/F0','percentile',10)')
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
                obj.disp_info('Some timescale are estimated based on scan command, but not measured. This may cause errors in the exact timescale', 2);
            end
            timescale.tp                = cellfun(@(x) x.analysis_params.timepoints, obj.arboreal_scans);
            timescale.durations         = timescale.sr.*timescale.tp; % same as cellfun(@(x) x.analysis_params.duration, obj.arboreal_scans)
            timescale.durations_w_gaps  = cellfun(@(x) x.analysis_params.timescale_w_gaps(end), obj.arboreal_scans);
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
        
        function tp = get.tp(obj)
            %% Return the number of imaging timepoint
            % -------------------------------------------------------------
            % Syntax:
            %   timescale = obj.timescale;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   t (1xT DOUBLE)
            %       equivalent to sum(obj.timescale.tp). Use it to
            %       simplify your code
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   28/02/2023

            tp = sum(cellfun(@(x) x.analysis_params.timepoints, obj.arboreal_scans));
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
                obj.disp_info('breakpoint field added post hoc. plase re-extract',3)
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
            %   This should match obj.ref.indices.n_tree_ROIs
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            if ~obj.use_hd_data
                n_ROIs = size(obj.ref.simple_data,2);
            else
                n_ROIs = sum(obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1));
            end
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
            %   This should match obj.ref.indices.n_pop_ROIs
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            n_pop_ROIs = size(obj.ref.simple_pop_data,2);
        end
        function n_rec = get.n_rec(obj)
            %% Total Number of records in the experiment
            % -------------------------------------------------------------
            % Syntax:
            %   n_rec = obj.n_rec;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   n_ROIs (INT)
            %       Number of records in the experiment
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   30/03/2023

            n_rec = numel(obj.arboreal_scans);
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

        function global_median_rescaled = get.global_median_rescaled(obj)
            %% Median trace of all non-excluded data, after rescaling
            % -------------------------------------------------------------
            % Syntax:
            %   global_median_rescaled = obj.global_median_rescaled;
            % -------------------------------------------------------------
            % Inputs:
            % -------------------------------------------------------------
            % Outputs:
            %   global_median_rescaled (Tx1 SINGLE)
            %       median rescale trace
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022

            global_median_rescaled = nanmedian(obj.rescaled_traces, 2);
            if isempty(global_median_rescaled)
                obj.disp_info('Global_median_rescaled not available until you called obj.rescale_traces',2);
            end
        end
        
        function binned_data = get.binned_data(obj)
            binned_data = obj.binned_data;
            if isempty(binned_data)
                obj.disp_info('Binned data not available because no groups were set. Using unique group will all segments instead',1);
                binned_data = {};
                binned_data.condition   = 'single group';
                binned_data.groups      = {1:obj.n_ROIs};
                binned_data.metrics     = 1;
                binned_data.bin_legend  = {'all ROIs'};
                binned_data.readmap     = sort(unique([binned_data.groups{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
            end
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
            external_variables = cellfun(@(x) x.external_var, obj.arboreal_scans, 'UniformOutput', false);
        end


        %% ########### ... ###################
        function disp_info(~, message, level)
            %% Print readable Info, warning and error messages
            % -------------------------------------------------------------
            % Syntax:
            %   obj.disp_info(message, level)
            % -------------------------------------------------------------
            % Inputs:
            %   message (char or cell array)
            %       The message to print. a tab heading is automatically
            %       added. If you want to print multiple lines, pass a cell
            %       array. Tabulation and new lines are automatically
            %       added.
            %   level (INT)
            %       the level controls the appearence of the message. 0 is
            %       regular printing, 1 is black bold, 2 is orange bold, 3
            %       is red bold and 4 is red bold and interrupts the code
            %       by returning an error message
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            % * for a simple message use
            %       obj.disp_info('blabla',level).
            % * for several lines use a cell array
            %             obj.disp_info({'line 1,'...
            %                            '\t- line 2 with indent',...
            %                            ['line ',num2str(3),' with var']},1);
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   14/04/2022
            
            if nargin < 3 || isempty(level)
                level = 0;
            elseif level < 0 || level > 4                
                error('message level must be between 0 and 4. See obj.disp_info doc')
            else
                level = round(level);
            end
            
            if iscell(message)
               message  = strcat('\t',message,'\n');
               if numel(message) > 1
                   message(2:end)  = strcat(repmat('\t',1,level+1),message(2:end));
               end
               message  = [message{:}];
            end
            
            if level == 0
                fprintf([message,'\n']);
            elseif level == 1
                fprintf(['\t<strong>INFO : ',message,'\n</strong>']);
            elseif level == 2
                fprintf(1, ['\t[\b<strong>WARNING : ',message,'</strong>\n]\b']);
            elseif level == 3                
                fprintf(2, ['\t<strong>MAJOR WARNING : ',message,' \n</strong>\n']);
            elseif level == 4                
                fprintf(2, ['\t<strong>!!PROCESSING ERROR!! : ',message,' \n</strong>\n']);
                error(message)
            end
        end
        
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
            obj.is_detrended = true;
        end

        %% #############################

        function prepare_binning(obj, condition, demo)
            %% Creates bins of ROIs based on morphometric criterions
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.prepare_binning(condition, demo)
            % -------------------------------------------------------------
            % Inputs:
            %   condition (STR or {STR, INT} Cell or {1xN INT MAtrix}) - 
            %   Optional - Default is 'single group' (no binning)
            %       Defines the type of binning required, and for some
            %       binning, the size of the bins. 
            %       - if '' or 'single_group' or 'none' or 'all', no
            %       binning is done and there is a single global average
            %      - If the input is a cell array, then see extra notes
            %      arboreal_scan.get_ROI_groups documentation. 
            %       - If input is a matrix, or a cell array with a matrix,
            %       then we use a customized binning.
            %   demo (BOOL) - Optional - Default obj.demo
            %       If > 0, additioan lfiures are genrated
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            %   * The condition chosen is set in Value is set in
            %       EXPE.binned_data.condition
            %   * If obj.demo is > 0, additional figures are genrated
            %   showing the binning
            %
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   23/08/2022
            
            %% todo : move 'custom' into get_roi_groups, and add doc there

            if nargin < 3 || isempty(demo)
                demo = obj.demo;
            end
            
            %% Clear fields depending on a different scaling
            obj.rescaling_info              = {};
            obj.need_update(:)              = true;
            if nargin < 2 || (ischar(condition) && any(strcmp(condition, {'','single_group','none','all'})))
                obj.binned_data.condition   = 'single group';
                obj.binned_data.groups      = {1:obj.n_ROIs};
                obj.binned_data.metrics     = 1;
                obj.binned_data.bin_legend  = {'all ROIs'};
            elseif iscell(condition) && ischar(condition{1})
                obj.binned_data.condition   = condition;

                %% Define current binning rule. See arboreal_scan.get_ROI_groups for more info % CC and peak extractions are based on this binning
                [obj.binned_data.groups, obj.binned_data.metrics, obj.binned_data.bin_legend] = obj.ref.get_ROI_groups(obj.binned_data.condition, demo);
            elseif iscell(condition) && ismatrix(condition{1})
                obj.binned_data.condition   = 'custom';
                obj.binned_data.groups      = condition;
                obj.binned_data.metrics     = 1:numel(condition);
                legends = strcat('group ', num2str(1:numel(condition))');
                obj.binned_data.bin_legend  = cellstr(legends(1:3:end,:))';
            else
                obj.disp_info('Binning condition not identified. Condition must be a manual cell array of ROIs, or a valid method as defined in arboreal_scans.get_ROI_groups()',4)
            end
            obj.binned_data.readmap         = sort(unique([obj.binned_data.groups{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
            obj.set_median_traces(false);
            obj.disp_info('Bins were created. obj.binned_data is now available including group medians',1);
        end



        function pxl_list = get_voxel_for_ROI(obj, ROIs)
            pxl_list = [];
            res_list = obj.ref.header.res_list(:,1);
            pxls = [cumsum(res_list) - res_list(1) + 1;  sum(res_list)];
            for el = 1:numel(ROIs)
                pxl_list = [pxl_list, pxls(ROIs(el)):(pxls(ROIs(el)+1)-1)];
            end
        end

        function norm_cumsum = get_events_statistics(obj)
            if isfield(obj.event, 'fitting')
                peaks = obj.event.fitting.post_correction_peaks;
            else
                peaks = obj.binned_data.median_traces(vertcat(obj.event.peak_time{:}), :);
            end

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
                    h = histogram(peaks(:,gp),0:bin_size:max_peaks, 'FaceAlpha', 0.8,'EdgeColor','none', 'FaceColor', cmap(gp, :));hold on;
                end
                obj.binned_data.event_ditribution{gp} = [h.BinEdges', [h.BinCounts,NaN]'];
                title(obj.binned_data.bin_legend(gp));
            end

            %% Plot all individual peaks to see if they covary
            figure(1006);cla();plot(peaks);legend(obj.binned_data.bin_legend);title('peaks values per subgroups'); xlabel('event #'); ylabel('Amplitude');set(gcf,'Color','w');

            %% Plot the mean amplitude per subgroup of peaks, per distance
            %bin_step = 10;
            %norm_cumsum = cumsum(obj.binned_data.median_traces) ./ nanmax(cumsum(obj.binned_data.median_traces));
            norm_cumsum = cumsum(peaks) ./ nanmax(cumsum(peaks));
            obj.binned_data.norm_cumsum = norm_cumsum;
            figure(1007);cla();plot(norm_cumsum); hold on;set(gcf,'Color','w');xlabel('Event #');ylabel('normalized cumulative amplitude')
            title('cumulative sum of peaks'); ylim([0, 1]);legend(obj.binned_data.bin_legend,'Location','southeast');
            obj.disp_info('Detected events statistics computed. See obj.binned_data for the group-specific values',1);
        end

        function norm_vmr = compute_variability(obj, metric)   
            if nargin < 2 || isempty(metric)
            	metric = obj.variability_metric;
            else
                obj.variability_metric = metric;
            end
            
            %% Get summary covariance plots of the raw traces
            obj.compute_similarity();
            
            %% Get peaks time or, if available decay-corrected peak times
            if isfield(obj.event, 'fitting')
                peaks = obj.event.fitting.post_correction_peaks;
                times = obj.event.fitting.peak_times;
            else
                peaks = obj.binned_data.median_traces(vertcat(obj.event.peak_time{:}), :);
                times = obj.t(vertcat(obj.event.peak_time{:}));
            end
            
            %% Comput variability on peaks
            if strcmpi(metric, 'vmr')            
                vmr = nanvar(peaks,[],2)./nanmean(peaks, 2);
            elseif strcmpi(metric, 'CV')   
                cv  = nanstd(peaks,[],2)./nanmean(peaks, 2); 
            else
                obj.disp_info('error, valid variability metrics are vmr and cv',4)
            end

            [~, idx] = sort(vmr,'descend');

            %figure(123); cla();ylim([0, nanmax(obj.event.fitting.post_correction_peaks(:))]); hold on;
            %     for event = idx'
            %         show_event(obj.binned_data.median_traces, round(obj.event.fitting.peak_times/sr), event);
            %         drawnow;%pause(0.1)
            %     end

            figure(1009);cla();plot(peaks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')

            R = max(range(obj.binned_data.median_traces));
            norm_vmr = vmr/mean(vmr);
            obj.variability.index_of_disp = norm_vmr;
            figure(1010);clf();
            ax1 = subplot(2,1,1);plot(obj.t, obj.binned_data.median_traces); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
            ax2 = subplot(2,1,2);plot(times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
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

            figure(1011);cla();scatter(nanmedian(peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w');
            obj.disp_info('General variability statistics have been computed. See obj.variability',1)
        end
        
        function precision = compute_similarity(obj, use_bins)
            if nargin < 2 || isempty(use_bins)
                use_bins = true;
            end
            
            win                         = ceil(1./median(diff(obj.t)));
            if use_bins
                data                    = obj.binned_data.median_traces;
                obj.variability.source  = 'binned data';
            else
                data                    = obj.rescaled_traces;
                obj.variability.source  = 'ROIs';
            end
            if size(data,2) > 1
                comb                        = nchoosek(1:size(data,2),2);
                corr_results                = {};
                for pair = 1:size(comb,1)
                    corr_results{pair}      = movcorr(data(:,comb(pair, 1)),data(:,comb(pair, 2)),[win, 0]);
                end

                corr_results                = cell2mat(corr_results);
                precision                   = 1./nanvar(corr_results,[], 2);
                out                         = nanmax(precision(:)) / 10;
                precision(precision > out)  = NaN;
                obj.variability.corr_results= corr_results;
                obj.variability.precision   = precision;
            else
                obj.variability.corr_results= NaN(size(data));
                obj.variability.precision   = NaN(size(data));
            end
            if obj.rendering
                obj.plot_similarity();arrangefigures([1,2]);
            end
            obj.disp_info(['Event variability statistics have been computed using ',obj.variability.source,'. See obj.variability'],1)
        end
        
        function set.variability_metric(obj, variability_metric)
            % add a check for valid methods 
            if ~strcmpi(variability_metric, obj.variability_metric)
                obj.variability = {};
            end
            obj.variability_metric = variability_metric;
        end
        
        function usable_swc_ROIs = get.usable_swc_ROIs(obj)
            usable_swc_ROIs = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list));
        end
  
        function [tp, beh, analysis_mode] = get_tp_for_condition(obj, analysis_mode)
            beh             = {};  
            
            %% If you pass manual set of inputs, use that
            if isnumeric(analysis_mode)
                tp          = false(1,obj.tp);
                tp(analysis_mode) = true;
                beh = {};
                analysis_mode = 'manual';
                return
            end
            
            %% Get signal time range     
            tp              = true(1,obj.tp);
            using_peaks     = false;
            analysis_mode   = strrep(analysis_mode, '~active', 'quiet');
            if contains(analysis_mode, 'peaks') %% event time
                tp              = ~tp; %set all timepoints to false
                %tp(obj.event.peak_time{:}) = true
                %tp(obj.event.fitting.pre_correction_peaks) = true;
                tp(vertcat(obj.event.peak_time{:})) = true;% could be using obj.event.fitting.pre_correction_peaks
                using_peaks     = true;
            elseif contains(analysis_mode, 'quiet') %% low corr window
                tp_of_events    = sort(unique([obj.event.t_win{:}]));                
                tp(tp_of_events)= false;
            elseif contains(analysis_mode, 'active') %% high corr window
                tp_of_events    = sort(unique([obj.event.t_win{:}]));
                tp              = ~tp;
                tp(tp_of_events)= true;
            else
                %% keep all tp
            end  
            
            invert          = contains(analysis_mode, '~');
            analysis_mode   = erase(analysis_mode, {'peaks','quiet','active','~','_'});  
            
            analysis_mode_no_win = strsplit(analysis_mode,{'[',']'});
            analysis_mode_no_win = analysis_mode_no_win(1:3:end);
                  
            if ~all(cellfun(@isempty, analysis_mode_no_win)) && any(contains(obj.behaviours.types, analysis_mode_no_win, 'IgnoreCase',true))
                to_test = [];
                for el = obj.behaviours.types
                    if contains(analysis_mode_no_win, el{1},'IgnoreCase',true) || contains(el{1},analysis_mode_no_win,'IgnoreCase',true)
                        to_test(end+1) = 1;
                    else
                        to_test(end+1) = 0;
                    end
                end
                beh_name = {obj.behaviours.types{find(to_test)}};
 
                %% Check if there is window suffix
                beh_end_loc = strfind(analysis_mode, beh_name{1}) + numel(beh_name{1});
                if ~isempty(beh_end_loc) && beh_end_loc < numel(analysis_mode) && strcmp(analysis_mode(beh_end_loc), '[')
                    loc             = strfind(analysis_mode,']');
                    loc             = loc(find(loc > beh_end_loc, 1, 'first'));
                    range           = analysis_mode(beh_end_loc:loc);
                    bout_extra_win  = str2num(range);
                else
                    bout_extra_win  = obj.bout_extra_win;
                end

                %% Get behaviour bouts
                [~, ~, beh]         = obj.get_behaviours(beh_name); 
                [~, ~, active_tp]   = obj.get_activity_bout(beh_name, true, invert, '', bout_extra_win);
                
                %% Valid tp are either (in)active behaviour, or peaks during behaviours
                if using_peaks
                    tp = tp & active_tp{1};
                else
                    tp = active_tp{1};
                end
            elseif invert
                tp = ~tp;
            end
        end

        function get_spike_trains(obj)
            %% Use ML spike to get a spike train estimate
            if ~exist('spk_autocalibration.m','file') || ~exist('fn_getfile.m','file')
                obj.disp_info('You need to download the "bricks" and "ml_spikes" toolboxes, and add then to the path',4)
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
        
        function cross_validate(obj)
            if nargin < 2 || isempty(cross_validate)
                cross_validate                  = false;
            end
            
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
            %obj.get_dimensionality(false, n_factor, timepoints, dim_red_type, weigthed_average)
        end
        
        function [weighted_averages, data] = get_dimensionality(obj, data, dim_red_type, timepoints, n_components, varargin)
            
            % obj.get_dimensionality() always removes data rows that have too many NaNs.  They end up flagged as ~valid_trace_idx. This includes
            % rows that are fully NaN's 
            % rows that have more than 4 times the median number of NaNs (basically, low quality rows that would flicker on and off because of posthoc MC for example)
            % 
            % The rows being deleted depends on your usage of the function.
            % If you directly provide the data, this is the only filter applied, so it is up to you to discard the bad rows, either from the start, or by NaNing them (this will impact valid_trace_idx sequence, but not the final result)
            % If you do not provide the data, data is obtained from obj.rescaled_traces, which always automatically NaN obj.bad_ROI_list
            if nargin < 2 || isempty(data)
                data                = obj.rescaled_traces;
            end            
            if nargin < 3 || isempty(dim_red_type)
                dim_red_type        = 'pca';
            end
            if nargin < 4 || isempty(timepoints)
                timepoints          = 1:size(data,1);
            end
            
            %% Adjust timepoints
            median_subtracted = false;            
            if ischar(timepoints)
                original_variable = timepoints;
                variable    = timepoints;
                if contains(variable, 'subtracted')
                    median_subtracted = true;
                    variable = erase(variable, 'subtracted');
                end
                timepoints  = obj.get_tp_for_condition(variable);
                variable  	= original_variable;
            else
                variable = 'manual_input';
            end                
           
            %% Filter timepoints
            data                    = data(timepoints,:);

            %% Reset dimensionality field
            if nargin < 5 || isempty(n_components)
                [~,~,~,~,explained]     = pca(obj.rescaled_traces(:, ~all(isnan(obj.rescaled_traces), 1)));
                n_components            = find(cumsum(explained)  > 90, 1, 'first');
                obj.disp_info(['n_components was not specified. Number of factors was estimated as the number of PCA components required to explain 90 %% of the variance (',num2str(n_components),') components)'],2)
            end
            obj.dimensionality                  = {};
            obj.dimensionality.n_factors        = n_components;
            obj.dimensionality.dim_red_type     = dim_red_type;
            obj.dimensionality.variable         = variable;

            %% Remove NaN
            data(isinf(data))   = NaN;
            all_ROIs            = 1:size(data, 2);
            normal_n_NaN        = median(sum(isnan(data(:,~all(isnan(data)))))) * 4; % get an indicative number of NaN in a normal traces, and set acceptable thr at 4 times that
            valid_trace_idx     = sum(isnan(data)) <= normal_n_NaN; % exclude traces with too many NaNs (eg. traces that got masked completely)
            data                = fillmissing(data(:, valid_trace_idx),'spline'); % removed funny traces

            %% Subtract signal median if required
            if median_subtracted
                type_of_trace = 'median subtracted';
                data = data - nanmedian(data,2);
            else
                type_of_trace = 'original signal';
            end
            
            %% Need to set it now
            obj.dimensionality.valid_trace_idx = valid_trace_idx;

            %% Get single or multiple factor estimate
            T = [];
            stats = {};
            specVarPM = [];

            switch dim_red_type
                case 'pca'
                    [LoadingsPM, F, specVarPM,stats,explained,mu] = pca(double(data),'NumComponents',obj.dimensionality.n_factors);
                    T                                           = [];
                case 'nnmf'
                    N_REPLICATES  = 20
                    [F, LoadingsPM, D]  = nnmf(double(data),obj.dimensionality.n_factors,'replicates',N_REPLICATES,'algorithm','als');
                    LoadingsPM          = LoadingsPM';
                    T                   = [];
                    stats               = {};
                    specVarPM           = [];
                case 'factoran'
                    [LoadingsPM, specVarPM, T, stats, F] = factoran(double(data), obj.dimensionality.n_factors,'rotate','varimax'); % varimax, quartimax
                case 'phate'
                    varargin = {'ndim', obj.dimensionality.n_factors, varargin{:}};
                    [LoadingsPM, P, K]  = phate(double(data'),varargin{:});
                    T                   = [];
                    stats               = {};
                    specVarPM           = [];
                    F                   = [];
            end

            %% Store results
            obj.dimensionality.LoadingsPM           = LoadingsPM;       % loadings / PC / modes / Factors / PHATE modes
            obj.dimensionality.specVarPM            = specVarPM;
            obj.dimensionality.T                    = T;                % Rotation matrix
            obj.dimensionality.stats                = stats;            % Factoran stats
            obj.dimensionality.F                    = F;                % components
            obj.dimensionality.all_ROIs             = all_ROIs;         % first occurences
            obj.dimensionality.mask                 = timepoints;

            obj.dimensionality.cluster_idx          = [];
            obj.dimensionality.clust_meth           = [];
            obj.dimensionality.N_clust              = [];
            obj.dimensionality.clust_groups         = {};

            %% Plot location of strongest component
            if obj.rendering
                obj.plot_factor_tree();
            end
            
            %% Plot a map of tree weights by ROI number for each component
            weighted_averages = obj.get_weight_map();
            
            %% Plot weight-tree for each component
            if obj.rendering
                if obj.dimensionality.n_factors > 9
                    obj.disp_info('Only the first 9 dimensions were displayed. To see more type obj.plot_dim_tree(1:obj.dimensionality.n_factors)',1)
                end
                obj.plot_dim_tree(1:min(obj.dimensionality.n_factors, 9));
            end
            obj.disp_info({['Dimensionality reduction done using ',obj.dimensionality.dim_red_type],...
                            ['analysis done on ',num2str(obj.dimensionality.n_factors),' factors'],...
                            ['Data obtained from : ',type_of_trace,' traces'],...
                            ['using the condition filter : ',obj.dimensionality.variable]},1)
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

%         function [tree, soma_location, tree_values, values] = plot_strongest_comp_tree(obj, n_dim, fig_handle)
%             if nargin < 2 || isempty(n_dim)
%                 n_dim = size(obj.dimensionality.LoadingsPM, 2);
%             end
%             if nargin < 3 || isempty(fig_handle)
%                 fig_handle = 10200;
%             end
% 
% 
%             [~, loc]    = nanmax(obj.dimensionality.LoadingsPM(:,1:n_dim),[],2);
%             Valid_ROIs  = find(obj.dimensionality.valid_trace_idx);
%             values      = NaN(size(Valid_ROIs));
%             for ROI = 1:numel(Valid_ROIs)
%                 values(ROI)     = loc(ROI);
%             end
%             values = values(~isnan(values));
% 
%             %% Map dimension weights on the tree
%             if obj.rendering
%                 [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, 'Location of strongest component','',fig_handle, 'regular', 'jet');
%                 colorbar('Ticks',1:nanmax(values));colormap(jet(nanmax(values)))
%             end
%         end

        %% #############################################

        function update_all_signals(obj, new_method)
            if nargin < 2 || isempty(new_method)
                new_method = obj.extraction_method;
            elseif ~any(strcmp(new_method, {'max', 'mean', 'median','min'}))
                obj.disp_info('Only max, min, mean and median method are supported',4)
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
                obj.need_update(rec)                = true;
            end
            error_box('COMPRESSION MODES WERE UPDATED BUT YOU NEED TO SAVE THE ARBOREAL SCANS TO KEEP THIS CHANGE FOR NEXT RELOADING')
        end

        function process(obj, condition, time_smoothing, rendering)
            %% Call all processing steps
            if nargin < 2 || isempty(condition)
                condition = {'distance',Inf};
            end
            if nargin > 3 && ~isempty(time_smoothing)
                obj.time_smoothing = time_smoothing;
            end
            if nargin >= 4 && ~isempty(rendering)
                obj.rendering = rendering;
            end

            %% Remove previous plots
            obj.clear_plots();

            %% Load and concatenate traces for the selected experiment
            % obj.load_extracted_data();   % this also sets the current expe #

            %% Prepare binning of ROIs based on specific grouping condition
            obj.prepare_binning(condition); % eet obj.demo = 1 or obj.prepare_binning(condition, true); to display the figure wirh the bins

            %% Find peaks based on amplitude AND correlation
            obj.find_events();

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces(); % note that signal rescaling is computed on peaks, not on the entire trace

            %% Create median trace per bins
            obj.set_median_traces();

            %% Correct for decay to avoid overestimating peak amplitude
            %obj.fit_events();

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.get_events_statistics();

            %% Compute signal variability
            obj.compute_variability();

            %% Check how correlated are the peaks between different parts of the tree
            obj.get_correlations();

            %% Dimensionality reduction
            obj.get_dimensionality([],'phate','peaks_subtracted',[]); 
            
            %% Clustering
            obj.cluster_factors('dbscan', [])           

            %% Optionally, if external variables need an update, do it here
            %obj.update_external_metrics(60)

            %% Extract behaviours
            obj.get_behaviours();

            %% Spike inference
            try
            %    obj.get_spike_trains();
            end

            %% Save figure and/or analysis
            obj.auto_save_analysis = false;
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
