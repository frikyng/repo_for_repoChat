classdef arboreal_scan_plotting < handle
    %% subclass of arboreal_scan_experiment 
    properties
        rendering = true;
    end
    
    methods
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

        function plot_median_traces(obj, smoothing)
            if nargin < 2 || isempty(smoothing)
                smoothing = [0,0];
            end
            %% Plot the mean trace for each bin  
            figure(1001);cla();
            traces = smoothdata(obj.binned_data.median_traces, 'gaussian',smoothing);
            plot(obj.t, traces); hold on;
            legend(obj.binned_data.bin_legend); hold on;
            title('median scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
            figure(1023);plot(traces - nanmean(traces,2));hold on;title('group traces median - overall median');xlabel('time (s)');set(gcf,'Color','w');
        end
        
        function plot_similarity(obj)
            if size(obj.binned_data.median_traces,2) > 1
                f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                f.Tag   = 'Similarity plot'; %for figure saving
                ax1     = subplot(3,1,1); hold on;plot(obj.binned_data.median_traces);title('Cell Signal');ylim([-50,range(obj.binned_data.median_traces(:))]);
                ax2     = subplot(3,1,2); hold on;plot(obj.variability.corr_results,'Color',[0.9,0.9,0.9]);
                hold on;plot(nanmean(obj.variability.corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(obj.variability.corr_results, 1)],[-1,-1],'k--');plot([0,size(obj.variability.corr_results, 1)],[1,1],'k--')
                ax3     = subplot(3,1,3); hold on;plot(obj.variability.precision);title('Precision (1/Var)');set(ax3, 'YScale', 'log');plot([0,size(obj.variability.precision, 1)],[1,1],'k--')
                linkaxes([ax1,ax2, ax3],'x')
            else
                f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                annotation('textbox','String','N/A','FitBoxToText','on','FontSize',16,'FontWeight','bold')
            end
        end 
        
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

