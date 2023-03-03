%% Note, go to the process() method for an overview of the available features


%% TO DO :
% - improve update system when changing folder.
% - enable loading if we just have a folder with arboreal_scans objects
% - add warning if source folder is the actual raw experiment
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
        is_detrended    = false; % is set to True once you ran the detrending once.
        is_rescaled     = false; % is set to True once you ran the rescaling step.
        default_handle
        rescaling_method= 'by_trials_on_peaks'; 
        breakpoints     = []; % if you had a disruptive event during th experiment, a first scaling is done with large blocks
        use_hd_data     = false;

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
        updated_data_path       % If update_folder is used, the updated filpath
        extracted_traces        % Concatenated version of each obj.arboral_scan.simple_data
        extracted_pop           % Concatenated version of each obj.arboral_scan.simple_pop_data
        extracted_traces_conc   % Concatenated version of extracted_traces
        rescaled_traces         % Rescaled Traces according to rescaling_info
        extracted_pop_conc      % Concatenated version of extracted_pop
        global_median_raw       % The median of extracted_traces_conc
        global_median_rescaled  % The median of rescaled traces
        t                       % Pointer to obj.timescale.global_timescale
        tp                      % Imaging timpoints
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
            %   source_folder (STR) - Optional - Default is current folder
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

            if nargin < 1 || (~ischar(source_folder) && isnan(source_folder))
                return % empty object when building whole dataset                
            elseif isempty(source_folder)
                source_folder = pwd;
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
                            error('Error detected during extraction. Check that the settings.txt file is present in the top_folder or manually indicated, and that it contains the correct paths')
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
                for field = {'arboreal_scans','binned_data','rescaling_info','event','variability','dimensionality'}
                    obj.(field{1}) = {};
                end

                %% Rebuild from extracted arboreal_scans
                for el = fliplr(find(obj.need_update))
                    add_tree(el, keep_2D);
                end
                obj.need_update(:) = false;
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
            obj.rescaling_method= 'by_trials_on_peaks';
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
                error('Cutoff must be between 0 and 1')
            else
                obj.bad_ROI_thr = value;
                a = dbstack();
                if ~contains([a.name], 'arboreal_scan_experiment.reset') % not useful when resetting or initializing
                  %  obj.find_bad_ROIs();
                end
            end
        end
        
        function set.filter_win(obj, filter_win)
            %% Defines the gaussian filtering window to use across analyses
            % -------------------------------------------------------------
            % Syntax:
            %   obj.filter_win = filter_win;
            % -------------------------------------------------------------
            % Inputs:
            %   filter_win (FLOAT OR 2x1 FLOAT)
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
            
            
            if all(isnumeric(filter_win)) && numel(filter_win) == 2 && all(filter_win == obj.filter_win)
                %% no change, pass                
            elseif all(isnumeric(filter_win)) && numel(filter_win) == 1 || numel(filter_win) == 2 
                try
                    filter_win(filter_win < 0)  = filter_win(filter_win < 0) * nanmedian(1./obj.timescale.sr);
                end
                filter_win                  = abs(round(filter_win));                
                if numel(filter_win)    == 1                 
                    filter_win = [filter_win, filter_win];
                end
                obj.filter_win              = filter_win;
                if isfield(obj.binned_data, 'median_traces') && ~isempty(obj.binned_data.median_traces) %empty upon initialization, in which case we do not care
                    warning('Changing the filter window affects several preprocessing steps. Metaanalysises fields were reset')
                    obj.reset();
                end
            else
                error('filter window must be a set of one (for symmetrical gaussian kernel) or 2 (for asymetrical gaussian kernel) values. Values are rounded. If values are < 11, window is converted in seconds')
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
                    warning('Changing the detrending method affects several preprocessing steps. Meta analysises fields were reset')
                    obj.reset();
                end                
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
            extracted_traces_conc(isinf(extracted_traces_conc))    = NaN;
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
        end
        
        
        function binned_data = get.binned_data(obj)
            binned_data = obj.binned_data;
            if isempty(binned_data)
                warning('binned data not available because no groups were set. using unique group instead')
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

        function bad_ROI_list = get.bad_ROI_list(obj)
            bad_ROI_list = obj.bad_ROI_list;
            if ischar(bad_ROI_list) && strcmpi(bad_ROI_list, 'unset')
                warning('Uncorrelated ROIs were never set. Detecting now. To ignore this set value to []')
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
                warning('unuable to use HD data as it is not embedded in the current arboreal_scan_experiment object. Rebuild the object with HD data using expe = arboreal_scan_experiment([arboreal_scans folder PATH],true).')
            else
                obj.use_hd_data = false;
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
            obj.is_detrended = true;
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
                timepoints = 1:size(values,1);
            end
            if nargin < 4 || isempty(ROIs)
            	ROIs =  obj.ref.indices.complete_ROIs_list(:,1);
            end
            if nargin < 5 || isempty(color_range)
            	color_range = [nanmin(values(:)),nanmax(values(:))];
            end

            obj.ref.animate_tree(values, timepoints, ROIs, color_range)
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
            if  nargin < 2 || (ischar(condition) && any(strcmp(condition, {'','single_group','none','all'})))
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
                error('binning condition not identified')
            end
            obj.binned_data.readmap         = sort(unique([obj.binned_data.groups{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
            obj.set_median_traces(false);
        end

        function rescale_traces(obj, method, smoothing)
            if nargin >= 2 && ~isempty(method)
                obj.rescaling_method = method;
            end
            smoothing = 0;
            
            if isempty(obj.event)
                warning('LD RESCALING REQUIRES DETECTED EVENTS. You must run obj.find_events() first');
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
                error('rescaling method not valid, It must specify "global" or "trials"')
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
                warning('TRACES HAVE NOT BEEN RESCALED YET - CALL obj.rescale_traces() first.\n');
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

        function pxl_list = get_voxel_for_ROI(obj, ROIs)
            pxl_list = [];
            res_list = obj.ref.header.res_list(:,1);
            pxls = [cumsum(res_list) - res_list(1) + 1;  sum(res_list)];
            for el = 1:numel(ROIs)
                pxl_list = [pxl_list, pxls(ROIs(el)):(pxls(ROIs(el)+1)-1)];
            end
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
            if isfield(obj.event, 'fitting')
                peaks = obj.event.fitting.post_correction_peaks;
                times = obj.event.fitting.peak_times;
            else
                peaks = obj.binned_data.median_traces(vertcat(obj.event.peak_time{:}), :);
                times = obj.t(vertcat(obj.event.peak_time{:}));
            end
            
            %sr = nanmedian(diff(obj.timescale{expe}.global_timescale));
            vmr = nanvar(peaks,[],2)./nanmean(peaks, 2);
            %cv  = nanstd(obj.event.fitting.post_correction_peaks,[],2)./nanmean(obj.event.fitting.post_correction_peaks, 2); % (maybe chack snr at one point?  mu / sigma)
            %fano = []; % windowed VMR. usually for spike trains
            [~, idx] = sort(vmr,'descend');

            %figure(123); cla();ylim([0, nanmax(obj.event.fitting.post_correction_peaks(:))]); hold on;
            %     for event = idx'
            %         show_event(obj.binned_data.median_traces, round(obj.event.fitting.peak_times/sr), event);
            %         drawnow;%pause(0.1)
            %     end

            figure(1009);cla();plot(peaks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')

            R = max(range(obj.binned_data.median_traces));
            %norm_vmr = vmr/range(vmr);
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

            figure(1011);cla();scatter(nanmedian(peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
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
            %        timepoints when cell is active (i.e. during bAPs)
            %      - 'peaks' : correlation between all ROIS, only
            %        at peak time.
            %      - 'encoder_peaks' : correlation between all
            %        ROIS, only at peak time, AND only when running speed
            %        is high
            %      - '~encoder_peaks' : correlation between all
            %        ROIS, only at peak time, AND only when running speed
            %        is low
            %   * The reference (located at index 1 of the matrix) is
            %    either
            %       - The averaged somatic data when available (as defined
            %         by ref.indices.somatic_ROIs), or the nexus signal
            %         when somatic ROIS are not available
            %       - The behaviour data if you specify 'behref', AND IF
            %         you pass a valid behaviour in obj.cc_mode (list
            %         available when typing obj.behaviours.types). eg :
            %         '~encoder_peaks_behref'
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
            %        binned_data.median_traces).
            %   * To include population data, include 'pop' in cc_mode
            % -------------------------------------------------------------
            % Example - How To
            %
            % * Get correlations across all ROIs
            %   EXPE.get_correlations(''); % Rejected ROIs appear in gray
            %
            % * Get correlations across all ROIs
            %   EXPE.get_correlations('groups');
            %
            % * Get correlations across all ROIs, only at peak time
            %   EXPE.get_correlations('peaks');
            %
            % * Get correlations across binned traces, only at peak time
            %   EXPE.get_correlations('peaks_groups');
            %
            % * Get correlations only when the cell is active
            %   EXPE.get_correlations('active');
            %
            % * Get correlations only when the mouse is running (all
            %   timepoints when running)
            %   EXPE.get_correlations('encoder');
            %
            % * Get correlations only when the mouse is running (only peak
            %   times when running)
            %   EXPE.get_correlations('encoder_peaks');
            %
            % * Get correlations only when the mouse is NOT running (only
            %   peak  times when not running)
            %   EXPE.get_correlations('~encoder_peaks');
            %
            % * Get correlations only when the mouse is running (all
            %   timepoints when running), but use the behaviour as a ref
            %   for correlation (instead of the somatic data)
            %   EXPE.get_correlations('encoder_behref');
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022

        	if ~strcmp(obj.cc_mode, cc_mode)
                obj.crosscorr = []; %if you change the cc mode, clear the previous correlation results
            end

            original_cc_mode = cc_mode;
            msg = '\n\t';
            if contains(cc_mode, 'groups')
                msg = [msg, 'Correlation done by group median, using your current binning, '];
            else
                msg = [msg, 'Correlation done between all ROIs, '];
            end
            cc_mode = erase(cc_mode, {'groups', 'ROIs'});
            if contains(obj.cc_mode, 'pop')
                msg = [msg, 'including population data.'];
            end
            cc_mode = erase(cc_mode, 'pop');
            if contains(cc_mode, 'behref')
                msg = [msg, '\n\tReference trace is indicated behaviour, at the indicated timpoints. '];
            else
                msg = [msg, '\n\tReference trace is the somatic signal average, at the indicated timpoints. '];
            end
            cc_mode = erase(cc_mode, 'behref');
            if contains(cc_mode, 'peaks')
                msg = [msg, '\n\tCorrelation computed using data at peak time only, '];
            elseif contains(cc_mode, 'active')
                msg = [msg, '\n\tCorrelation computed using timepoints when the cell is active (i.e. during bAPs), '];
            elseif contains(cc_mode, 'quiet')
                msg = [msg, '\n\tCorrelation computed using timepoints when the cell is quiet (i.e. between bAPs), '];
            else
                msg = [msg, '\n\tCorrelation done using all timepoints, '];
            end
            cc_mode = erase(cc_mode, {'peaks', 'active', 'quiet', '_',' '});
            if ~isempty(cc_mode)
                if contains(cc_mode, '~')
                    msg = [msg, ' WHEN "',erase(cc_mode, '~'),'" behaviour IS NOT ongoing.'];
                else
                    msg = [msg, ' WHEN "',cc_mode,'" behaviour IS ongoing.'];
                end
            end

            fprintf([msg ,  '\n'])
            obj.cc_mode = original_cc_mode;
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
            else
                fprintf('Using current obj.cc_mode\n')
            end

            %% Get CC (and build the correlation matrix if the settings changed)
            cross_corr = obj.crosscorr;

            %% Show tree if required
            if obj.rendering
                obj.plot_correlation_results(cross_corr);
            end
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
            tp          = true(1,obj.tp);
            using_peaks = false;
            if contains(analysis_mode, 'peaks') %% event time
                tp              = ~tp;
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
                [~, ~, active_tp]   = obj.get_activity_bout(beh_name, true, [], invert, '', bout_extra_win);
                
                %% Valid tp are either (in)active behaviour, or peaks during behaviours
                if using_peaks
                    tp = tp & active_tp{1};
                else
                    tp = active_tp{1};
                end
            end
        end

        function tp = set_crosscorr(obj, cc_mode)
            %% Compute cross correlation
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.set_crosscorr(cc_mode)
            % -------------------------------------------------------------
            % Inputs:
            %   cc_mode (STR) - See Description for details - Optional -
            %       Default is [];
            %       If provided, update EXPE.cc_mode. see set.cc_mode for
            %       more details
            % -------------------------------------------------------------
            % Outputs:
            %   tp(1 x N INT)
            %       The timepoints used based on your fitlering criteria
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

            if nargin > 1
                obj.cc_mode = cc_mode; % change cc mode
            else
                cc_mode = obj.cc_mode;
                fprintf(['Using current obj.cc_mode : ',cc_mode,'\n'])
            end

            %% Get signal to use
            if contains(obj.cc_mode, 'groups')
                signal = obj.binned_data.median_traces;
            else% if contains(obj.cc_mode, 'ROIs')
                try
                    if ~obj.use_hd_data
                        signal = obj.rescaled_traces(:,obj.ref.indices.valid_swc_rois);
                    else
                        if contains(obj.cc_mode, 'ROIs')
                            error('to use ROIs, set obj.use_hd_data to false');
                        end
                        signal = obj.rescaled_traces;
                        warning('full hd not filtering excluded ROIS / voxels for now')
                    end
                catch
                    error('unable to get rescaled data. try to run obj.rescale_traces()')
                end
            end           
            cc_mode = erase(cc_mode, {'ROIs','groups'});
            
            %% Subtract median if required
            if contains(obj.cc_mode, 'subtracted')
                signal = signal - nanmedian(signal,2);
            end
            cc_mode = erase(cc_mode, {'subtracted'});

            %% Get ref ROIs and trace
            if ~obj.use_hd_data
                somatic_ROIs= obj.ref.indices.somatic_ROIs;
            else
                somatic_ROIs= obj.ref.indices.HD_somatic_ROIs;
            end
            ref         = nanmean(obj.rescaled_traces(:, somatic_ROIs),2); %always ref, unless you pass 'behref'

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
            cc_mode = erase(cc_mode, {'pop'});

            %% Get timepoints base on filter
            [tp, beh, cc_mode] = obj.get_tp_for_condition(cc_mode);
            if iscell(tp)
                tp = any([vertcat(tp{:})]); % if multipe behaviours, using BITWISE OR
            end

            %% Update ref if we want directly the behaviour data instead of the somatic ROIs
            if ~isempty(beh) && contains(cc_mode, 'behref')
                ref         = beh.value';
            end

            %% Build the arrays used for the correlation matrix
            variable        = [ref(tp, :), signal(tp, :)];
            if ~isempty(pop)
                variable = [variable, pop(tp,:)];
            end

            cc   = corrcoef(variable,'Rows','Pairwise')';
            obj.crosscorr = cc;
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

            if isempty(obj.binned_data) && contains(obj.cc_mode, 'groups')
                fprintf('\tcrossscorr cannot be calculated using the "groups" flag if no groups were defined. Use obj.prepare_binning first\n')
                crosscorr = [];
                return
            elseif isempty(obj.event) && contains(obj.cc_mode, 'peaks')
                fprintf('\tcrossscorr cannot be calculated using the "peaks" flag if you have not run event detection. Use obj.detect_events first\n')
                crosscorr = [];
                return
            elseif isempty(obj.crosscorr) || size(obj.crosscorr, 1) ~= numel(obj.binned_data)
                obj.set_crosscorr();
            end
            crosscorr = obj.crosscorr;
        end

        function [tree, soma_location, tree_values, mean_bin_cc] = plot_corr_tree(obj, cc, ref_column)
            if nargin < 2 || isempty(cc)
                cc = obj.crosscorr;
                if contains(obj.cc_mode,'pop')
                    pop_sz = size(obj.extracted_pop_conc,2);
                    cc = cc(1:(end-pop_sz),1:(end-pop_sz));
                end
            end
            if nargin < 3
                ref_column = [];
            end

            %% Erase diagonal
            cc(1:size(cc,1)+1:end)= NaN;

            %% Remove "ref" row/column when you pass a matrix with one row per ROI
            if size(cc,1) > 2 && (size(cc,1) == (obj.ref.indices.n_tree_ROIs+1) || (~isempty(obj.binned_data) && size(cc,1) == (numel(obj.binned_data.groups)+1)) && (isempty(ref_column) || ~all(ref_column == 1)))
                cc = cc(2:end,2:end);
            end

            %% If no ref were provided, we will use the average correlation
            if isempty(ref_column)
                ref_column = 1:size(cc,1);
            end

            %%
            if contains(obj.cc_mode,'groups')
                %% Identify valid set of traces
                valid_gp            = find(~all(isnan(obj.binned_data.median_traces))); % You get NaN'ed bins if the soma location is not scanned (eg a big pyramidal cell)

                %% Build tree values per bin
                mean_bin_cc     = [];
                ROIs_list       = [];
                if ~isempty(cc)
                    for gp = 1:numel(obj.binned_data.groups)
                        roi_of_gp           = obj.binned_data.groups{gp};
                        %roi_of_gp           = roi_of_gp(~ismember(roi_of_gp, obj.bad_ROI_list)); %% COMMENT OUT TO INCLUDE BAD ROIS
                        v_of_gp             = cc(ref_column,gp);
                        ROIs_list           = [ROIs_list, roi_of_gp];
                        mean_bin_cc         = [mean_bin_cc, repmat(nanmean(v_of_gp), 1, numel(roi_of_gp))];
                    end
                end
            else
                ROIs_list   = obj.ref.indices.valid_swc_rois;
                mean_bin_cc = cc(:,find(obj.dimensionality.valid_trace_idx, 1, 'first'));
            end

            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(mean_bin_cc, ROIs_list, obj.default_handle, 'Correlation with most proximal segment','',1018);
            caxis([0,1]);
            col = colorbar; col.Label.String = 'Spatial correlation between ROIs/groups with soma';

%
%             [S,Q] = genlouvain(double(cc),[],[],1);
%             [a,b] = sort(S);
%             figure(1008);clf();imagesc(cc(b,b))
%             obj.ref.plot_value_tree(S,'','','','',124,'ddd','lines');
%             R = unique(S);
%             colorbar('Ticks',R);
%             caxis([nanmin(R)-0.5, nanmax(R)+0.5])
%             colormap(lines(numel(unique(S))));
%
%             figure(125);clf();
%             for community = R'
%                 plot(nanmean(obj.extracted_traces_conc(:,S == community),2)); hold on;
%             end
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
                variable = timepoints;
                if contains(variable, 'subtracted')
                    median_subtracted = true;
                    variable = erase(variable, 'subtracted');
                end
                timepoints = obj.get_tp_for_condition(variable);
            else
                variable = 'manual_input';
            end                
           
            %% Filter timepoints
            data                    = data(timepoints,:);

            %% Reset dimensionality field
            if nargin < 5 || isempty(n_components)
                [~,~,~,~,explained]     = pca(obj.rescaled_traces(:, ~all(isnan(obj.rescaled_traces), 1)));
                n_components            = find(cumsum(explained)  > 90, 1, 'first');
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
                data = data - nanmedian(data,2);
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
                    warning('only the first 9 dimensions were displayed. To see more type obj.plot_dim_tree(1:obj.dimensionality.n_factors)')
                end
                obj.plot_dim_tree(1:min(obj.dimensionality.n_factors, 9));
            end
        end

        function cluster_factors(obj, clust_meth, N_clust)
            if nargin < 2 || isempty(clust_meth)
                clust_meth                    = 'hierarchical';
            end
            if nargin < 3 || isempty(N_clust)
                N_clust                    = [];
            end
            
            %% Clear previous results
            obj.dimensionality.cluster_idx      = [];
            obj.dimensionality.clust_meth       = clust_meth;
            obj.dimensionality.N_clust        	= N_clust;
            obj.dimensionality.clust_groups     = {};
            obj.dimensionality.epsilon       	= [];
            obj.dimensionality.labels           = [];

            %% If N cluster is 0, skip clustering
            if obj.dimensionality.N_clust == 0
                obj.dimensionality.cluster_idx  = (1:size(obj.dimensionality.LoadingsPM, 1))';
                obj.dimensionality.sorted_idx   = (1:size(obj.dimensionality.LoadingsPM, 1))';
                obj.dimensionality.clust_groups = {1:size(obj.dimensionality.LoadingsPM, 1)};
                return
            end

            %% If N cluster is unknow, try to guess
            if isempty(obj.dimensionality.N_clust) || isnan(obj.dimensionality.N_clust)
                R = 1:20;
                if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'kmeans','silhouette','KList',R);
                    obj.dimensionality.N_clust = eva.OptimalK;
                    fprintf(['Optimal number of clusters is ',num2str(obj.dimensionality.N_clust), '\n'])
                elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical')
                    eva = evalclusters(obj.dimensionality.LoadingsPM,'linkage','silhouette','KList',R);
                    obj.dimensionality.N_clust = eva.OptimalK;
                    fprintf(['Optimal number of clusters is ',num2str(obj.dimensionality.N_clust), '\n'])
                elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                    obj.dimensionality.epsilon = test_epsilon(obj, obj.dimensionality.LoadingsPM);
                    fprintf(['Optimal espilon s ',num2str(obj.dimensionality.epsilon), '\n'])
                elseif strcmp(obj.dimensionality.clust_meth, 'strongest')
                    obj.dimensionality.N_clust = NaN;
                else
                    error(['No auto auto-determination of the number of cluster for this method\n'])
                end
            elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                if obj.dimensionality.N_clust > 0
                    obj.dimensionality.epsilon = obj.dimensionality.N_clust;
                    obj.dimensionality.N_clust = [];
                else
                    [obj.dimensionality.epsilon, obj.dimensionality.N_clust] = test_epsilon(obj, obj.dimensionality.LoadingsPM,[],[],[],obj.dimensionality.N_clust);
                end
            end
            
            if strcmp(obj.dimensionality.clust_meth, 'kmeans')
                cluster_idx = kmeans(obj.dimensionality.LoadingsPM , obj.dimensionality.N_clust);
                figure(3);clf(); silhouette(obj.dimensionality.LoadingsPM,cluster_idx)
            elseif strcmp(obj.dimensionality.clust_meth, 'hierarchical') % see https://fr.mathworks.com/help/stats/hierarchical-clustering.html
                cluster_idx = clusterdata(obj.dimensionality.LoadingsPM,'Linkage', 'ward', 'MAXCLUST', obj.dimensionality.N_clust);%, 'Criterion','distance' 'MAXCLUST', 40)
                figure(3);clf(); silhouette(obj.dimensionality.LoadingsPM,cluster_idx)
                %Y = pdist(obj.dimensionality.LoadingsPM ,'euclidean');Z = linkage(Y,'ward');figure();dendrogram(Z);                                
            elseif strcmp(obj.dimensionality.clust_meth, 'dbscan')
                MIN_GP = 4
                cluster_idx                 = dbscan(obj.dimensionality.LoadingsPM, obj.dimensionality.epsilon, MIN_GP, 'Distance', 'euclidean');
                obj.dimensionality.N_clust  = numel(unique(cluster_idx(cluster_idx > 0)));
                if any(cluster_idx <= 0)
                    col = UNASSIGNED_ROI_COLOR;
                else
                    col = [];
                end
                col                         = [col ; jet(obj.dimensionality.N_clust)];
                obj.dimensionality.labels   = [cluster_idx, cluster_idx, cluster_idx];
                count                       = 1;
                for v = unique(cluster_idx)'
                    subset                                  = obj.dimensionality.labels(:,1) == v;
                    obj.dimensionality.labels(subset, :)    = repmat(col(count,:),sum(subset),1);
                    count                                   = count + 1;
                end
                figure(1111);clf();
                for offset = 0:(min(size(obj.dimensionality.LoadingsPM, 2)-3, 2))
                    subplot(2,2,offset+1); hold on; grid; hold on
                    scatter3(obj.dimensionality.LoadingsPM (:,1+offset),obj.dimensionality.LoadingsPM (:,2+offset),obj.dimensionality.LoadingsPM (:,3+offset),20,obj.dimensionality.labels, 'filled');
                end               
            elseif strcmp(obj.dimensionality.clust_meth, 'strongest')
                %% To assign to strongest component
                for row = 1:size(LoadingsPM,1)
                    [~, maxloc]                     = max(LoadingsPM(row, :));
                    LoadingsPM(row, :)              = 0;
                    LoadingsPM(row, maxloc)         = 1;
                end
            else
                obj.dimensionality.cluster_idx      = [];
                obj.dimensionality.clust_groups     = {};
            end
            
            %% Sort clusters by number of elements
            if obj.dimensionality.N_clust > 0
                [~, gp] = sort(hist(cluster_idx,unique(cluster_idx)), 'descend');
            else
                gp = 1;
            end
            a = unique(cluster_idx)';
            gp = a(gp);
            idx_sorted = NaN(size(cluster_idx));
            count = 1;
            for gp_idx = gp(gp > 0)
                idx_sorted(cluster_idx == gp_idx) = count;
                count = count + 1;
            end
            idx_sorted(isnan(idx_sorted))       = 0;
            obj.dimensionality.cluster_idx      = idx_sorted;       % stored group ids (reordered by number of element)
            [~, obj.dimensionality.sorted_idx]  = sort(idx_sorted); % get ROI index to reorder the data by group

            %% Now, build groups
            obj.dimensionality.clust_groups     = {};
            for el = 1:obj.dimensionality.N_clust
                obj.dimensionality.clust_groups{el} = find(obj.dimensionality.cluster_idx == el);
            end
            
            %% Plot clusters
            if obj.rendering                
                obj.plot_cluster_tree();
            end
        end


        function [tree, soma_location, tree_values, values] = plot_dim_tree(obj, comp, fig_handle, cmap, tree_type)
            % check 58, % noise issue 64 'D:/Curated Data/2019-09-24/experiment_1/18-13-20/'
            if nargin < 2 || isempty(comp)
                comp        = 0;
            end
            if nargin < 3 || isempty(fig_handle)
                fig_handle = 2000;% fig number or fig hande
            end
            if nargin < 4 || isempty(cmap)
                cmap       = 'redbluesymmetrical';
            end
            if nargin < 5 || isempty(tree_type)
                tree_type  = 'simple';
            end
            figure(fig_handle);clf();

            %% If you asked more dimensions than tavailable diemsnions, clip the list
            if any(comp)
                comp(comp > size(obj.dimensionality.LoadingsPM, 2)) = [];
            end
            
            %% Recover Factors / Loadings
            LoadingsPM  = obj.dimensionality.LoadingsPM;
            Valid_ROIs  = obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx);
            values      = NaN(numel(comp), numel(Valid_ROIs));
            [~, loc]    = nanmax(obj.dimensionality.LoadingsPM,[],2);
            
            %% Prepare rendering
            N_Dim = numel(comp);
            n_row = floor(sqrt(N_Dim));
            n_col = ceil(numel(comp) / n_row);

            %% Plot components (or strongest factor location if comp == 0)
            tiledlayout(n_row, n_col, 'Padding', 'none', 'TileSpacing', 'none'); 
            for comp_idx = 1:numel(comp)
                ax = nexttile; % a bit beter than subplot, but if you have matlab < 2019b, you can use the line below
                %ax = subplot(n_row, n_col, comp_idx);
                dim = comp(comp_idx);
                for roi = 1:numel(Valid_ROIs)
                    if dim
                        values(comp_idx, roi) = LoadingsPM(roi, dim);
                    else                        
                        values(comp_idx, roi) = loc(roi);
                    end
                end

                %% Map dimension weights on the tree
                if obj.rendering || ishandle(fig_handle)
                    if dim
                        titl = ['Component ',num2str(dim), ' weights'];
                    else
                        titl = 'Location of strongest component';
                        cmap = jet(nanmax(values));
                        cmap = cmap(values, :);
                    end  
                    
                    if obj.use_hd_data
                        [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(split_values_per_voxel(values(comp_idx, :), obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1), signal_indices), '','',['phate #',num2str(dim),' Loadings (per voxel)'],'',ax,tree_type,cmap);
                    else
                    	[f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values(comp_idx, :), Valid_ROIs,'',titl,'',ax,tree_type,cmap);
                    end
                end
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
                obj.need_update(rec)                = true;
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
            
%             sz = vertcat(obj.event.peak_value{:});
%             thr1 = max(sz) * 0.8;
%             thr2 = max(sz) * 0.2;
%             mask = cellfun(@(x) mean(x) < thr1 & mean(x) > thr2, obj.event.peak_value);            
%             for fn = fieldnames(obj.event)'                
%                  if numel(obj.event.(fn{1})) == numel(mask)
%                      obj.event.(fn{1}) = obj.event.(fn{1})(mask);
%                  end
%             end

            %% Identify and log poorly correlated ROIs
            obj.find_bad_ROIs(correlation_res, idx_filter);
        end

        function [bad_ROIs, mean_corr_with_others] = find_bad_ROIs(obj, correlation_res, ROIs)
            if nargin < 3 || isempty(ROIs)
                ROIs = 1:obj.n_ROIs;
            end      
            mean_corr_with_others   = [];
            bad_ROIs                = [];
            if nargin < 2 || isempty(correlation_res)
                warning('BAD ROI IDENTIFICATION RELIES ON EVENT DETECTION. EVENT DETECTION FIELD SEEMS EMPTY.  run obj.detect_events() first');                
                return
            else
                events = obj.event;
            end
            THR_FOR_CONNECTION      = 0.2 % Defines what level of minimal pairwise correlation means "these two ROIs are connected"
            fprintf(['* Now detecting ROIs that are either,\n'...
                     '      - so poorly correlated to the rest of the tree that they probably belong to another cell (or have no signal).\n',...
                     '      - are member of batch_params.excluded_branches \n',...
                     '* Threshold for exclusion is  ',num2str(THR_FOR_CONNECTION),' %% \n'])

            %% Show mean correlation with each ROI
            max_corr_for_cell       = max(events.globality_index(2:end)); % QQ 1st point sometimes show some artifacts
            for key = 1:numel(ROIs)
                corr_results_sub    = correlation_res.corr_results(events.t_corr(events.is_global), correlation_res.comb(:,2) == key | correlation_res.comb(:,1) == key);
                corr_results_sub    = [corr_results_sub(:,1:(key-1)), NaN(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];
                mean_corr           = nanmean(corr_results_sub,1);
                mean_corr_with_others(key) = sum(mean_corr > THR_FOR_CONNECTION);
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
            RECOVERY_THR        = 1 - THR_FOR_CONNECTION
            excl_but_not_bad    = was_excluded & ~obj.bad_ROI_list';
            excl_but_good       = was_excluded & (mean_corr_with_others_norm > RECOVERY_THR)';
            if any(excl_but_good)
                 fprintf(['!!! ROIs ',num2str(find(excl_but_good')),' was/were excluded but seem highly correlated\n'])
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

            %% Prepare binning of ROIs based on specific grouping condition
            obj.prepare_binning(condition); % eet obj.demo = 1 or obj.prepare_binning(condition, true); to display the figure wirh the bins

            %% Find peaks based on amplitude AND correlation
            obj.find_events();

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces(); % note that signal rescaling is computed on peaks, not on the entire trace

            %% Create median trace per bins
            obj.set_median_traces()

            %% Get summary covariance plots of the raw traces
            obj.compute_similarity();

            %% Correct for decay to avoid overestimating peak amplitude
            %obj.fit_events();

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.get_events_statistics();

            %% Study bAP heterogeneity
            obj.assess_variability()

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
