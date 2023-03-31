classdef correlation_analysis < handle
    %% Subclass of arboreal_scan_experiment
    properties
    	cc_mode                 = 'groups_peaks';% Type of data used to compute cross correlation. 
        crosscorr                               % Correlation matrix of peaks/signal across ROIs/groups during/between bAps/activity_bouts
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

            obj.disp_info(msg,1);
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
                obj.disp_info(['Correlation based on ',obj.cc_mode],1)
            end

            %% Get CC (and build the correlation matrix if the settings changed)
            cross_corr = obj.crosscorr;

            %% Show tree if required
            if obj.rendering
                obj.plot_correlation_results(cross_corr);
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
                obj.disp_info(['Correlation based on ',obj.cc_mode],1)
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
                            obj.disp_info('To use ROIs in correlation analysis,set obj.use_hd_data to false',4)
                        end
                        signal = obj.rescaled_traces;
                        obj.disp_info('Full HD not filtering excluded ROIS / voxels for now',3)
                    end
                catch
                    obj.disp_info('unable to get rescaled data. try to run obj.rescale_traces()',4)
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
                obj.disp_info('crossscorr cannot be calculated using the "groups" flag if no groups were defined. Use obj.prepare_binning first',3)
                crosscorr = [];
                return
            elseif isempty(obj.event) && contains(obj.cc_mode, 'peaks')
                obj.disp_info('crossscorr cannot be calculated using the "peaks" flag if you have not run event detection. Use obj.detect_events first',3)
                crosscorr = [];
                return
            elseif isempty(obj.crosscorr) || size(obj.crosscorr, 1) ~= numel(obj.binned_data)
                obj.set_crosscorr();
            end
            crosscorr = obj.crosscorr;
        end
    end
end
