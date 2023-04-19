classdef correlation_analysis < handle
    %% Subclass of arboreal_scan_experiment
    % See Demo_CrossCorrelation.mlx for examples    
    properties
    	cc_mode                 = 'groups_peaks';% Type of data used to compute cross correlation. 
        crosscorr                                % Correlation matrix of peaks/signal across ROIs/groups during/between bAps/activity_bouts
        crosscorr_ref           = 'soma'         % Defines what is used to compute the reference column of the crosscorr
        cc_mode_label     
    end

    methods
        function set.cc_mode(obj, cc_mode)
            %% Set cross correlation mode
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.get_correlations(cc_mode)
            % -------------------------------------------------------------
            % Inputs:
            %   cc_mode (STR) - See Description for details - Optional -
            %           Default is your current obj.cc_mode value
            %       update EXPE.cc_mode with existing value. 
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
            %   * An ampty or invalid option will use the entire signal for
            %       all ROIs
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

            if numel(obj.arboreal_scans)
                obj.disp_info(obj.cc_mode_label, 1);
                obj.cc_mode = cc_mode;
            end
        end
        
        function msg = get.cc_mode_label(obj)
            temp_cc_mode = obj.cc_mode;
            msg = '';
            if contains(temp_cc_mode, 'groups')
                msg = [msg, 'Correlation done by group median, using your current binning, '];
            else
                msg = [msg, 'Correlation done between all valid ROIs, '];
            end
            temp_cc_mode = erase(temp_cc_mode, {'groups', 'ROIs'});
            if contains(obj.cc_mode, 'pop')
                msg = [msg, 'including population data.'];
            end
            temp_cc_mode = erase(temp_cc_mode, 'pop');
            if contains(temp_cc_mode, 'behref')
                msg = [msg, '\n\tReference trace is indicated behaviour, at the indicated timpoints. '];
            else
                if isequal(obj.crosscorr_ref, 1:obj.n_ROIs)
                    ref_name    = 'the Cell-wide averaged signal';
                elseif isequal(obj.crosscorr_ref, obj.ref.indices.somatic_ROIs) || isequal(obj.crosscorr_ref, 'soma')
                    ref_name    = 'the perisomatic averaged signal';
                elseif isnumeric(obj.crosscorr_ref) && numel(obj.crosscorr_ref(obj.crosscorr_ref ~= 0)) == 1
                    ref_name    = ['ROI # ',num2str(obj.crosscorr_ref)];
                end                
                msg = [msg, '\n\tReference trace is ',ref_name,', at the indicated timpoints. '];
            end
            temp_cc_mode = erase(temp_cc_mode, 'behref');
            if contains(temp_cc_mode, 'peaks')
                msg = [msg, '\n\tCorrelation computed using data at peak time only, '];
            elseif contains(temp_cc_mode, 'active')
                msg = [msg, '\n\tCorrelation computed using timepoints when the cell is active (i.e. during bAPs), '];
            elseif contains(temp_cc_mode, 'quiet')
                msg = [msg, '\n\tCorrelation computed using timepoints when the cell is quiet (i.e. between bAPs), '];
            else
                msg = [msg, '\n\tCorrelation done using all timepoints, '];
            end
            temp_cc_mode = erase(temp_cc_mode, {'peaks', 'active', 'quiet', '_',' '});
            if contains(temp_cc_mode, 'subtracted')
                msg = [msg, 'on median-subtracted traces'];
            end
            temp_cc_mode = erase(temp_cc_mode, {'subtracted', '_',' '});
            if ~isempty(temp_cc_mode)
                if contains(temp_cc_mode, '~')
                    msg = [msg, ' WHEN "',erase(temp_cc_mode, '~'),'" behaviour IS NOT ongoing.'];
                else
                    msg = [msg, ' WHEN "',temp_cc_mode,'" behaviour IS ongoing.'];
                end
            end
        end
        
        function crosscorr_ref = get.crosscorr_ref(obj)
            crosscorr_ref = obj.crosscorr_ref;
            if ischar(crosscorr_ref)
                obj.crosscorr_ref = crosscorr_ref; % will force convertion to numbers
            end
        end
        
        function set.crosscorr_ref(obj, crosscorr_ref)
            %% Defines the ROIs to use for the reference comun of the
            % correlation matrix
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.crosscorr_ref = crosscorr_ref
            % -------------------------------------------------------------
            % Inputs:
            %   crosscorr_ref ([] or INT or INT ARRAY or CHAR)
            %       If a a list of integer is provided, the ref will be
            %       calculated with the average signal from these ROIs.
            %       If a single value is passed, the reference comlumn wil
            %       lcorrespond to this ROi (in which case there will be 2
            %       identical columns in the crosscorr matrix)
            %       If the value is 0 or [], all ROIs will be used to
            %       computed the average
            %       If the value is 'soma', the average periosomatic signal
            %       will be used as a reference. 
            %       If the value is 'average', this is equivalent to 0 or
            %       [].
            % -------------------------------------------------------------
            % Outputs:
            % -------------------------------------------------------------
            % Extra Notes:
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   06/04/2023
            
            if isnumeric(crosscorr_ref)
                crosscorr_ref = round(crosscorr_ref);
                if ~any(crosscorr_ref  <= 0) && ~any(crosscorr_ref > obj.n_ROIs)
                    obj.crosscorr_ref = crosscorr_ref;
                elseif ~any(crosscorr_ref) || isempty(crosscorr_ref)
                    obj.crosscorr_ref = 1:obj.n_ROIs;
                else
                    obj.disp_info('Error : when passing manual values for crosscorr_ref, values must be >= 0 and <= obj.n_ROIs',4);
                end
            elseif ischar(crosscorr_ref) && strcmpi(crosscorr_ref, 'soma')
                if ~obj.use_hd_data
                    obj.crosscorr_ref = obj.ref.indices.somatic_ROIs;
                else
                    obj.crosscorr_ref = obj.ref.indices.HD_somatic_ROIs;
                end
            elseif ischar(crosscorr_ref) && strcmpi(crosscorr_ref, 'average') || isempty(crosscorr_ref)
                obj.crosscorr_ref = 1:obj.n_ROIs;
            else
                obj.disp_info('Reference for the crosscorr must be a ROI or a list of ROI to average, "average", 0 or [] for a cell wide average ROIs, or "soma" to use the average of perisomatic ROIs' ,4)
            end
        end

        function cross_corr = get_correlations(obj, cc_mode, crosscorr_ref)
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
            %   crosscorr_ref ([] or INT or INT ARRAY or CHAR)
            %       If a a list of integer is provided, the ref will be
            %       calculated with the average signal from these ROIs.
            %       If a single value is passed, the reference comlumn wil
            %       lcorrespond to this ROI (in which case there will be 2
            %       identical columns in the crosscorr matrix)
            %       If the value is 0 or [], all ROIs will be used to
            %       computed the average
            %       If the value is 'soma', the average periosomatic signal
            %       will be used as a reference. 
            %       If the value is 'average', this is equivalent to 0 or
            %       [].
            % -------------------------------------------------------------
            % Outputs:
            %   crosscorr (NxN DOUBLE)
            %       Cross correlation between ROIs OR groups based on
            %       traces OR Peaks, depending on the settings. see
            %       set.cc_mode doc for the deails
            % -------------------------------------------------------------
            % Extra Notes:
            %   * obj.bout_extra_win is updated whenever 
            % -------------------------------------------------------------
            % Author(s):
            %   Antoine Valera.
            %--------------------------------------------------------------
            % Revision Date:
            %   13/05/2022

            if nargin >= 2
                obj.cc_mode = cc_mode; % change cc mode                
            end
            obj.disp_info(['Current obj.cc_mode is : ',obj.cc_mode],1)
            if nargin >= 3
                obj.crosscorr_ref = crosscorr_ref; % change cc reference
            end
            obj.disp_info(['Current obj.crosscorr_ref is : ',obj.crosscorr_ref],1)

            %% Get CC (and build the correlation matrix if the settings changed)
            cross_corr = obj.crosscorr;

            %% Show tree if required
            if obj.rendering
                obj.plot_correlation_results(cross_corr);
            end
        end
        
        function set.crosscorr(obj, cc_mode)
            %% Set method that compute cross correlation
            % -------------------------------------------------------------
            % Syntax:
            %   EXPE.set.crosscorr(cc_mode)
            % -------------------------------------------------------------
            % Inputs:
            %   cc_mode (STR) - See Description for details - Optional -
            %       Default is [];
            %       If provided, update EXPE.cc_mode. see set.cc_mode for
            %       more details. Empty option will use the whole signal on 
            %       ROIs all
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

            if ~numel(obj.arboreal_scans)
                return
            end
            if ~ischar(cc_mode)
                cc_mode = '';
            end
            obj.cc_mode = cc_mode;
            if strcmp(cc_mode, 'none')
                obj.crosscorr = [];
            end
            obj.disp_info(['Correlation based on ',obj.cc_mode],1)

            %% Get signal to use
            if contains(obj.cc_mode, 'groups')
                signal = obj.binned_data.median_traces;
            else % 'ROIs', or no specific input
                %% Do some quick check to see if we can do it
                if obj.use_hd_data
                    if contains(obj.cc_mode, 'ROIs')
                        obj.disp_info('To use ROIs in correlation analysis,set obj.use_hd_data to false',4)
                    end                        
                    obj.disp_info('Full HD not filtering excluded ROIS / voxels for now',3)
                end
                
                %% get ref signal
                if contains(obj.cc_mode, 'subtracted')
                    signal = obj.rescaled_traces;
                    if isempty(signal)
                        obj.disp_info('unable to get rescaled data, which is required for median subtraction. try to run obj.rescale_traces()',4)
                    end
                else
                    signal = obj.extracted_traces_conc;
                end
                
            end           
            cc_mode = erase(cc_mode, {'ROIs','groups'});
            
            %% Subtract median if required
            if contains(obj.cc_mode, 'subtracted')
                signal = signal - obj.global_median_rescaled;
            end
            cc_mode = erase(cc_mode, {'subtracted'});

            %% Get ref ROIs and trace (unless you passed 'behref')
            ref         = nanmean(obj.extracted_traces_conc(:, obj.crosscorr_ref),2); % since correlation are scale insensitive, we do not need obj.rescaled_traces

            %% Add population signal if needed
            if contains(obj.cc_mode, 'pop')
                if isempty(obj.extracted_pop_conc)
                    obj.disp_info('You asked for population correlation, but there is population data for this recording (or it has not been extracted)',2)
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
                obj.disp_info('crossscorr cannot be calculated using the "groups" flag if no groups were defined. Use obj.prepare_binning first',4)
                obj.crosscorr = 'none';
                return
            elseif isempty(obj.event) && contains(obj.cc_mode, 'peaks')
                obj.disp_info('crossscorr cannot be calculated using the "peaks" flag if you have not run event detection. Use obj.find_events first',4)
                obj.crosscorr = 'none';
                return
            else
                obj.crosscorr = obj.cc_mode;
%             elseif isempty(obj.crosscorr) || size(obj.crosscorr, 1) ~= numel(obj.binned_data)
%                 obj.crosscorr = obj.cc_mode;
            end
            crosscorr = obj.crosscorr;
        end
    end
end
