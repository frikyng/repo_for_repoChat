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
            
        function plot_original_data(obj)
            figure(1055);clf();plot(obj.t, obj.extracted_traces_conc,'Color',[0.9,0.9,0.9]); hold on
            plot(obj.t, nanmean(obj.extracted_traces_conc,2),'r');            
        end

        function plot_median_traces(obj, rescaled) 
            if nargin < 2 || isempty(rescaled)
                rescaled = true;
            end
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
        
       function plot_rescaling_info(obj)
            n_gp = numel(obj.rescaling_info.individual_scaling{1});
            n_rec = numel(obj.rescaling_info.individual_scaling);
            figure(1012);cla();hold on;
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

            figure(1013);cla();hold on;set(gcf,'Color','w');
            xlim([1,n_gp]);hold on;
            title('global scaling factor and offset per ROI');
            xlabel('ROI');
            plot(obj.rescaling_info.scaling,'ko-'); hold on;ylabel('Consensus Scaling Factor');    
            yyaxis right;plot(obj.rescaling_info.offset,'ro-');ylabel('Consensus Baseline Percentile'); 
            ax = gca;
            ax.YAxis(1).Color = 'r';
            ax.YAxis(2).Color = 'k';
        end
        
        function plot_correlation_results(obj, cross_corr)
            if nargin < 2 || isempty(cross_corr)
                cross_corr = obj.crosscorr;
            end
            imAlpha=ones(size(cross_corr));
            imAlpha(isnan(cross_corr))=0;
            figure(1008);clf();imagesc(cross_corr, 'AlphaData',imAlpha); hold on;set(gcf,'Color','w');
            set(gca,'color',0.8*[1 1 1]);
            caxis([0,1]); hold on;xticks(1:size(cross_corr, 1));yticks(1:size(cross_corr, 1))
            colorbar; hold on;
            if contains(obj.cc_mode, 'pop')
                pop_label = num2cell(1:size(obj.extracted_pop_conc,2));
            else
                pop_label = [];
            end
            if contains(obj.cc_mode, 'groups')
                plot([1.5,1.5],[0.5,size(cross_corr, 1)+0.5],'k-');xticklabels(['Soma/Proximal seg',obj.binned_data.bin_legend]);xtickangle(45);
                plot([0.5,size(cross_corr, 1)+0.5],[1.5,1.5],'k-');yticklabels(['Soma/Proximal seg',obj.binned_data.bin_legend]);
                part_1 = ' groups';
                N_reg = numel(obj.binned_data.bin_legend)+1;
            else
                xticklabels(['Soma/Proximal seg', num2cell(obj.ref.indices.valid_swc_rois'),pop_label]);xtickangle(45);
                yticklabels(['Soma/Proximal seg',num2cell(obj.ref.indices.valid_swc_rois'),pop_label]);
                N_reg = numel(obj.ref.indices.valid_swc_rois)+1;
                part_1 = ' ROIs';
            end
            if contains(obj.cc_mode, 'pop')
                part_2 = ' and population';
                ax = gca;
                ax.XTickLabel((N_reg+1):end) = cellfun(@(x) ['\color{red}', x], ax.XTickLabel((N_reg+1):end),'uni',false);
                ax.YTickLabel((N_reg+1):end) = cellfun(@(x) ['\color{red}', x], ax.YTickLabel((N_reg+1):end),'uni',false);
            else
                part_2 = '';
            end
            title(['Correlation between',part_1,part_2]);
            arrangefigures(0);

            %% Project correlation value onto the tree
            obj.plot_corr_tree();
        end
        
        function plot_detected_events(obj)
            raw_med = nanmedian(obj.extracted_traces_conc, 2);raw_med = raw_med - prctile(raw_med,1); %% QQ CHCK WHAT@S THE THING USED IN DETTECT_EVETS
            figure(1035);clf(); subplot(2,1,1);plot(obj.t,raw_med , 'r'); hold on; plot(obj.t, obj.binned_data.global_median); hold on;scatter(obj.t(vertcat(obj.event.peak_time{:})), vertcat(obj.event.peak_value{:}), 'filled');xlabel('time (s)');ylabel('median signal (A.U)');
            subplot(2,1,2); plot(obj.t, obj.event.mean_pairwise_corr); hold on;scatter(obj.t(obj.event.t_corr), obj.event.mean_pairwise_corr(obj.event.t_corr), [], obj.event.globality_index, 'filled');xlabel('time (s)');ylabel('Mean pairwise correlation');set(gcf,'Color','w');
        end
        
        function plot_rescaled_traces(obj)
            valid_ROIS = ~ismember(1:obj.n_ROIs, obj.bad_ROI_list);
            rescaled_traces = obj.rescaled_traces(:,valid_ROIS);            
            figure(1034);cla();plot(obj.t, rescaled_traces);title('rescaled traces');xlabel('time (s)');ylabel('ROIs');set(gcf,'Color','w');
            figure(1033);cla();imagesc(rescaled_traces');caxis([prctile(reshape(rescaled_traces, [], 1),1), prctile(reshape(rescaled_traces, [], 1),99)]);title('rescaled traces 2D');xlabel('frames');ylabel('ROIs');set(gcf,'Color','w');
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
        

        function plot_cluster_tree(obj, tree_handle, map_handle, trace_handle)
            if nargin < 2 || isempty(tree_handle)
                tree_handle = figure(10200);tree_handle = gcf();
            end
            if nargin < 3 || isempty(map_handle)
                map_handle = figure(999);map_handle = gca();
            end
            if nargin < 4 || isempty(trace_handle)
                trace_handle = figure(7777);trace_handle = gca();
            end            
            
            tree_handle
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(obj.dimensionality.cluster_idx, find(obj.dimensionality.valid_trace_idx), obj.default_handle, 'Clusters','',tree_handle, 'regular', 'jet');
            colorbar('Ticks',1:nanmax(obj.dimensionality.cluster_idx));colormap(jet(nanmax(obj.dimensionality.cluster_idx)));  

            %% Display rearranged Loadings and Cluster limits    
            cla(map_handle);
            imagesc(map_handle, obj.dimensionality.LoadingsPM(obj.dimensionality.sorted_idx,:));caxis([0,1]);xlabel('Factors');hold(map_handle, 'on')
            axis(map_handle,'tight')
            for el = 1:numel(unique(obj.dimensionality.cluster_idx))
                start = find(obj.dimensionality.cluster_idx(obj.dimensionality.sorted_idx,:) == el, 1, 'last');
                if ~isempty(start) && obj.dimensionality.N_clust
                    plot(map_handle, [0.5,obj.dimensionality.n_factors+0.5],[start,start],'w--','Linewidth',2);hold(map_handle, 'on')
                end
            end
            title(map_handle, 'Factor/Loadings/Components');ylabel('ROI (sorted)')
            
            %% Display average weightings
            cla(trace_handle);
            cluster_traces = obj.get_cluster_traces();
            for gp = 1:size(cluster_traces, 1)
                plot(trace_handle, obj.t, cluster_traces(gp,:));hold(trace_handle, 'on')
            end
            title(trace_handle, 'Average traces per cluster');xlabel('Time (s)');colororder(trace_handle, jet(obj.dimensionality.N_clust))
        end
        
        function cluster_traces = get_cluster_traces(obj)
            rescaled_traces     = obj.rescaled_traces(:, obj.dimensionality.valid_trace_idx);
            cluster_traces      = [];
            for gp = sort(unique(obj.dimensionality.cluster_idx))'
                cluster_traces(gp,:) = nanmean(rescaled_traces(:,obj.dimensionality.cluster_idx == gp),2);
            end
        end

        function plot_factor_tree(obj, weights_to_show, weighted_averages, tree_handle, map_handle, trace_handle)
            % rplace : plot_dimensionality_summary
            if nargin < 2 || isempty(weights_to_show) % number or list of factor to display
                weights_to_show = 1:obj.dimensionality.n_factors;
            end
            if nargin < 3 || isempty(weighted_averages)
                rescaled_traces = obj.rescaled_traces(:, obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx));
                all_weights         = {}; weighted_averages   = [];
                for w = weights_to_show
                    all_weights{w}          = obj.dimensionality.LoadingsPM(:,w)/sum(obj.dimensionality.LoadingsPM(:,w));
                    weighted_averages(w, :) = nanmean(rescaled_traces'.* all_weights{w}, 1);
                end
            end
            if nargin < 4 || isempty(tree_handle)
                tree_handle = 10200;
            end
            if nargin < 5 || isempty(map_handle)
                map_handle = figure(1017);map_handle = gca();
            end
            if nargin < 6 || isempty(trace_handle)
                trace_handle = figure(7777);trace_handle = gca();
            end  
            
            %% Plot strongest component per ROI
            obj.plot_strongest_comp_tree('', tree_handle);
            
%             [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(obj.dimensionality.cluster_idx, find(obj.dimensionality.valid_trace_idx), obj.default_handle, 'Clusters','',tree_handle, 'regular', 'jet');
%             colorbar('Ticks',1:nanmax(obj.dimensionality.cluster_idx));colormap(jet(nanmax(obj.dimensionality.cluster_idx)));  


            figure(1024);cla();plot(obj.t(obj.dimensionality.mask), obj.dimensionality.F(:,weights_to_show));set(gcf,'Color','w');title(['Components ',strjoin(strsplit(num2str(weights_to_show),' '),'-'),' per ROI']);xlabel('t');
            
            
            %% Display rearranged Loadings and Cluster limits    
            cla(map_handle);
            imagesc(map_handle, obj.dimensionality.LoadingsPM(:,weights_to_show));caxis([0,1]);xlabel('Factors');hold(map_handle, 'on')
            %colorbar;set(gcf,'Color','w');
            axis(map_handle,'tight')
            title(map_handle, ['Components ',strjoin(strsplit(num2str(weights_to_show),' '),'-'),' per ROI']);
            ylabel('ROI');xlabel('component')

            
            %% Display average weightings
            cla(trace_handle);
            for w = weights_to_show
                plot(trace_handle, obj.t, weighted_averages(w, :));hold(trace_handle, 'on');
            end
            title(trace_handle, 'Signal average using Loadings as weight');xlabel('Time (s)');
            %legend();set(gcf,'Color','w');
        end
    end
end

