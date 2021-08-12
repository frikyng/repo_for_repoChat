% 1001 : Median scaled responses per group
% 1002 : Event detection on median
% 1003 (NO TITLE) : global distribution of event amplitudes
% 1004 (NO TITLE) : distribution of event amplitudes per region
% 1005 : Median Events and fitted decay
    % 10051-1005n : Median Events and fitted decay for group n
% 1006 : peaks values per subgroups
% 1007 : cumulative sum of peaks
% 1008 : CC matrix
% 1009 : event per index of dispersion
% 1010 : index of dispersion vs traces
% 1011 : index of dispersion vs amplitude
% 1012 : Scaling factor per ROI
% 1013 : Global scaling factor and offset per ROI
% 1014 : median tau per group
% 1015 : median tau 1 per event
% 1016 : tau vs event amplitude
% 1017 : components weight per ROI (matrix)
% 1018 : corr_tree
% 1019 – no title  individual events?
% 1020 : strongest component
    % 10201 – 1020n : individual components
% 1021 : Weighted signal average per component
% 1022 : variability assessment 
% 1023 : median traces  - overall median
% 1024 : …
% 1025 : Spike inference (pr bin)
% 1026 : Behaviour
% 1027 : Behaviour activity bouts
% 1028 : cross validation result - training
% 1029 : cross validation result - testing
% to find un-numbered figures look for "% FF"

classdef arboreal_scan_plotting < handle
    %% subclass of arboreal_scan_experiment 
    properties
        rendering = true;
    end
    
    methods
        function clear_plots(obj)
            close_list = [obj.get_fig_list(50,50), 124]; % adding temporary figures # 124
            figHandles = get(groot, 'Children');
            for f = 1:numel(figHandles)
                if ismember(figHandles(f).Number, close_list)
                	close(figHandles(f));
                end
            end      
        end
        
        function figure_list = get_fig_list(obj, n_groups, n_dims)   
            if nargin < 2 || isempty(n_groups)
                n_groups = numel(obj.binned_data.groups);
            end
            if nargin < 3 || isempty(n_dims)
                n_dims = size(obj.dimensionality.LoadingsPM, 2);
            end
            figure_list = [1001:1035, 10051:(10050 + n_groups), 10021:(10020+n_dims), 10200:(10200+n_dims),10830];            
        end
        
        function plot_rescaling_info(obj)
            n_gp = numel(obj.rescaling_info.individual_scaling{1});
            n_rec = numel(obj.rescaling_info.individual_scaling);
            figure(1012);clf();hold on;
            col = mat2cell(viridis(n_gp), ones(1,n_gp), 3);
            p = plot(1:n_rec, vertcat(obj.rescaling_info.individual_scaling{:}),'o-');set(gcf,'Color','w');
            arrayfun(@(x, y) set(x, 'Color', y{1}), p, col);
            set(gca, 'YScale', 'log');
            %legend(bin_legend);
            xticks(1:n_rec);
            xlim([0.8,n_rec+0.2]);
            ylabel('Scaling Factor');
            xlabel('Recording #');
            title('Scaling factor per ROI, per trial');

            figure(1013);clf();hold on;set(gcf,'Color','w');
            xlim([1,n_gp]);hold on;
            title('global scaling factor and offset per ROI');
            xlabel('ROI');
            plot(obj.rescaling_info.scaling,'ko-'); hold on;ylabel('Consensus Scaling Factor');    
            yyaxis right;plot(obj.rescaling_info.offset,'ro-');ylabel('Consensus Baseline Percentile'); 
            ax = gca;
            ax.YAxis(1).Color = 'r';
            ax.YAxis(2).Color = 'k';
        end

        function plot_median_traces(obj)            
            %% Plot the mean trace for each bin  
            if isempty(obj.binned_data)
               warning('Cannot plot binned data before having formed some groups. use obj.prepare_binning(condition)')
               return 
            end            
            traces = obj.binned_data.median_traces;
            figure(1001);cla();plot(obj.t, traces); hold on;
            legend(obj.binned_data.bin_legend); hold on;
            if obj.is_rescaled
                suffix = ' (rescaled)';
            else
                suffix = ' (raw)';
            end
            title(['median traces trace per group',suffix]);xlabel('time (s)');set(gcf,'Color','w');
            if obj.is_rescaled
                figure(1023);cla();plot(obj.t, traces - nanmean(traces,2));hold on;title('group traces median - overall median');xlabel('time (s)');set(gcf,'Color','w');
            end
        end
        
        function plot_detected_events(obj)
            raw_med = nanmedian(obj.extracted_traces_conc, 2);raw_med = raw_med - prctile(raw_med,1); %% QQ CHCK WHAT@S THE THING USED IN DETTECT_EVETS
            figure(1035);clf(); subplot(2,1,1);plot(obj.t,raw_med , 'r'); hold on; plot(obj.t, obj.binned_data.global_median); hold on;scatter(obj.t(vertcat(obj.event.peak_time{:})), vertcat(obj.event.peak_value{:}), 'filled');xlabel('time (s)');ylabel('median signal (A.U)');
            subplot(2,1,2); plot(obj.t, obj.event.mean_pairwise_corr); hold on;scatter(obj.t(obj.event.t_corr), obj.event.mean_pairwise_corr(obj.event.t_corr), [], obj.event.globality_index, 'filled');xlabel('time (s)');ylabel('Mean pairwise correlation');set(gcf,'Color','w');
        end
        
        function plot_rescaled_traces(obj)
            M = prctile(obj.rescaled_traces(:),90);
            valid_ROIS = ~ismember(1:obj.n_ROIs, obj.bad_ROI_list);
            figure(1034);cla();plot(obj.t, obj.rescaled_traces(:,valid_ROIS));title('rescaled traces');xlabel('time (s)');ylabel('ROIs');set(gcf,'Color','w');
            figure(1033);cla();imagesc(obj.rescaled_traces(:,valid_ROIS)');caxis([0, M]);title('rescaled traces 2D');xlabel('frames');ylabel('ROIs');set(gcf,'Color','w');
        end
        
        function plot_similarity(obj)
            if size(obj.binned_data.median_traces,2) > 1
                f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                f.Tag   = 'Similarity plot'; %for figure saving
                R       = [nanmin(obj.binned_data.median_traces(:)), nanmax(obj.binned_data.median_traces(:))];
                ax1     = subplot(3,1,1); hold on;plot(obj.binned_data.median_traces);title('Cell Signal');ylim([R(1)-(range(R/10)), R(2)+(range(R/10))]);
                ax2     = subplot(3,1,2); hold on;plot(obj.variability.corr_results,'Color',[0.9,0.9,0.9]);
                hold on;plot(nanmean(obj.variability.corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(obj.variability.corr_results, 1)],[-1,-1],'k--');plot([0,size(obj.variability.corr_results, 1)],[1,1],'k--')
                ax3     = subplot(3,1,3); hold on;plot(obj.variability.precision);title('Precision (1/Var)');set(ax3, 'YScale', 'log');plot([0,size(obj.variability.precision, 1)],[1,1],'k--')
                linkaxes([ax1,ax2, ax3],'x')
            else
                f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                annotation('textbox','String','N/A','FitBoxToText','on','FontSize',16,'FontWeight','bold')
            end
        end 
        
%         function plot_events(obj)
%             for region = 1:size(obj.events.fitting)
%                 figure();plot(squeeze(events(:,2,:)))
%             end
%             
%         end
        
        function plot_dimensionality_summary(obj, weights_to_show, weighted_averages)
                if nargin < 2 || isempty(weights_to_show) % number or list of factor to display
                weights_to_show = 1:obj.dimensionality.n_factors;
            end
            if nargin < 3 || isempty(weighted_averages)
                rescaled_traces = obj.rescaled_traces(:, obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx));
                all_weights         = {}; weighted_averages   = [];
                for w = weights_to_show
                    all_weights{w}          = obj.dimensionality.LoadingsPM(:,w)/sum(obj.dimensionality.LoadingsPM(:,w));
                    weighted_averages(w, obj.dimensionality.mask) = nanmean(rescaled_traces(obj.dimensionality.mask,:)'.* all_weights{w}, 1);
                end
            end
            figure(1024);cla();plot(obj.t(obj.dimensionality.mask), obj.dimensionality.F(:,weights_to_show));set(gcf,'Color','w');title(['Components ',strjoin(strsplit(num2str(weights_to_show),' '),'-'),' per ROI']);xlabel('t');
            figure(1017);cla();imagesc(obj.dimensionality.LoadingsPM(:,weights_to_show)); colorbar;set(gcf,'Color','w');title(['Components ',strjoin(strsplit(num2str(weights_to_show),' '),'-'),' per ROI']);xlabel('component');ylabel('ROI');
            figure(1021);clf();
            for w = weights_to_show
                hold on; plot(obj.t(obj.dimensionality.mask), weighted_averages(w, obj.dimensionality.mask));
            end
            legend();set(gcf,'Color','w');
            title('Weighted signal average per component')
        end
        
    end
end

