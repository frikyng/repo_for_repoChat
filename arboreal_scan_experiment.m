classdef arboreal_scan_experiment < handle    
    properties
        source_folder
        need_update
        extracted_data_paths % i.e. data_folders 
        arboreal_scans
        ref % first extracted arboreal scan, for conveniency

        rendering = true;
        demo = 0;
        
        filter_win      = 0;
        filter_type     = 'gaussian';
                general_info
        extracted_traces
        current_segmentation;
        binned_data
        timescale
        event_fitting
        crosscorr
        dimensionality
        external_variables
        
        default_handle
    end
    
    methods
        function obj = arboreal_scan_experiment(source_folder)
%             if nargin < 2 || isempty(export_folder)
%                 export_folder = [export_folder, '/meta/'];
%             end
            
            %% Fix paths
            obj.source_folder       = parse_paths(source_folder);
%             obj.export_folder       = parse_paths(export_folder);

%             %% Create export folder if necessary
%             if ~isdir(obj.export_folder)
%                 mkdir(obj.export_folder);
%             end

            %% Get all analysed data folder in the source folder, for the type specified
            all_recordings          = dir([obj.source_folder,'/**/*-*-*_exp_*_*-*-*']);
            all_recordings          = all_recordings(~[all_recordings(:).isdir]);
            all_recordings          = all_recordings(~(arrayfun(@(x) strcmp(x.name, '.'), all_recordings) | arrayfun(@(x) strcmp(x.name, '..'), all_recordings)));
            
            names                   = [vertcat(all_recordings.folder), repmat('/',numel(all_recordings), 1), vertcat(all_recordings.name)];
            obj.extracted_data_paths= cellstr(names)';
            obj.need_update         = false(1, numel(obj.extracted_data_paths));    
            
            obj.update();
        end 
        
        function update(obj)
            for rec = 1:numel(obj.extracted_data_paths)
                obj.arboreal_scans{rec} = load(obj.extracted_data_paths{rec});
                obj.arboreal_scans{rec} = obj.arboreal_scans{rec}.obj;
                obj.arboreal_scans{rec}.analysis_params.data = []; % clear full data. If you change t you'll need to rextract everything anyway
            end   
        end
        
        function extracted_traces = get.extracted_traces(obj)
            extracted_traces = cellfun(@(x) x.data, obj.arboreal_scans, 'UniformOutput', false);  
            extracted_traces = cellfun(@(x) x - prctile(x, 1), extracted_traces, 'UniformOutput', false);             
        end

        function timescale = get.timescale(obj)
            %% Prepare timescale for each recording and concatenated timescale
            timescale = {};
            timescale.sr  = 1./cellfun(@(x) x.analysis_params.measured_points_per_s, obj.arboreal_scans);
            if any(isnan(timescale.sr))
                warning('some timscale are estimated and not measured')
                estimated = 1./cellfun(@(x) x.analysis_params.estimated_points_per_s, obj.arboreal_scans);
                timescale.sr(isnan(timescale.sr))  = estimated(isnan(timescale.sr));
            end
            timescale.tp  = cellfun(@(x) x.analysis_params.timepoints, obj.arboreal_scans);
            timescale.t_start = [];
            timescale.durations = timescale.sr.*timescale.tp; % same as cellfun(@(x) x.analysis_params.duration, obj.arboreal_scans)
            timescale.rec_timescale     = arrayfun(@(x, y) linspace(0, y*x, x), timescale.tp, timescale.sr, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) diff(x), timescale.rec_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) [x(1), x], timescale.global_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cumsum(horzcat(timescale.global_timescale{:}));
        end
        
        function ref = get.ref(obj)            
            ref = obj.arboreal_scans{1};
        end
        
        
        function prepare_binning(obj, condition)
            if isempty(obj.current_segmentation) || isempty(obj.binned_data) || isempty(obj.current_segmentation.condition) || (~strcmp(obj.current_segmentation.condition{1}, condition{1}) || obj.current_segmentation.condition{2} ~= condition{2})
            	obj.current_segmentation = {};
                obj.current_segmentation.condition = condition;
            elseif strcmp(obj.current_segmentation.condition{1}, condition{1}) && obj.current_segmentation.condition{2} == condition{2}
                return
            end
            
            %% Define current binning rule. See arboreal_scan.get_ROI_groups for more info % CC and peak extractions are based on this binning    
            [obj.binned_data.groups, obj.binned_data.metrics, obj.binned_data.bin_legend] = obj.ref.get_ROI_groups(obj.current_segmentation.condition);
            obj.need_update(:)            = true; 
            obj.binned_data.readmap       = sort(unique([obj.binned_data.groups{:}])); % ROIs_per_subgroup_per_cond values corresponds to real ROIs, but not column numbers, so we need a readout map
        end
        
        function set.current_segmentation(obj, current_segmentation)
            obj.current_segmentation    = current_segmentation;            
            caller                      = dbstack('-completenames'); 
            if isempty(current_segmentation) && ~contains([caller.name],'MatFile') % detect reset
                warning('segmentation conditions were changed - this will reset all extracted fields relying on it')
                %% Clear dependant analyses
                %             obj.crosscorr{expe}                 = [];
                %             obj.event_fitting{expe}             = [];
                %             obj.general_info{expe}.offset       = zeros(1,size(obj.extracted_traces{expe}{1}, 2)); 
                %             obj.general_info{expe}.scaling      = ones(1,size(obj.extracted_traces{expe}{1}, 2)); 
            end
        end
        
        function rescale_traces(obj) 
            [best_scal_f, best_offset]    = scale_every_recordings(obj.extracted_traces, obj.demo);
            obj.general_info.scaling      = best_scal_f;
            obj.general_info.offset       = best_offset;
            if obj.rendering
                figure();plot(obj.timescale.global_timescale,nanmedian(obj.get_rescaled_traces(),2));xlabel('Scan Time (s)')
                arrangefigures([1,2]);
            end
        end
        
        function rescaled_traces = get_rescaled_traces(obj)
            %% Using obj.general_info.scaling & obj.general_info.offset, rescale the traces
            all_traces_per_rec = obj.extracted_traces;
            
            %% Now rescale each trace with a unique value across recordings (which would be spectific of that region of the tree).
            for trace = 1:size(obj.extracted_traces{1}, 2) 
                temp = cellfun(@(x) x(:, trace) ,all_traces_per_rec, 'UniformOutput', false); % 'end' or any group you want to test
                for rec = 1:numel(temp)
                    temp{rec}(1:2) = NaN; % first 2 points are sometimes a bit lower than the rest. We NaN them. This is also useful later on to identify trials
                end
                concat_trace    = vertcat(temp{:});
                off = prctile(concat_trace, obj.general_info.offset(trace));
                temp = cellfun(@(x) x - off, temp, 'UniformOutput', false);
                temp = cellfun(@(x) x / obj.general_info.scaling(trace), temp, 'UniformOutput', false);
                for rec = 1:numel(temp)
                    all_traces_per_rec{rec}(:, trace) = temp{rec};
                end
            end 
            rescaled_traces = vertcat(all_traces_per_rec{:});           
        end
        
        function set_median_traces(obj)            
            rescaled_traces = obj.get_rescaled_traces();            
            %rescaled_traces = rescaled_traces - prctile(nanmedian(wavelet_denoise(rescaled_traces),2),5);
            
            %% Create median trace per bins
            all_traces_per_bin = cell(1, numel(obj.binned_data.groups));
            for gp = 1:numel(obj.binned_data.groups)
                columns                     = ismember(obj.binned_data.readmap, obj.binned_data.groups{gp});
                all_traces_per_bin{gp}      = nanmedian(rescaled_traces(:,columns), 2);
            end    
            
            obj.binned_data.median_traces = cell2mat(all_traces_per_bin);%cellfun(@(x) cell2mat(x) ,all_traces_per_bin, 'UniformOutput', false); 
            
            if obj.rendering
                %% Plot the mean trace for each bin  
                figure(1001);cla();
                plot(obj.timescale.global_timescale, obj.binned_data.median_traces); hold on;
                legend(obj.binned_data.bin_legend); hold on;
                title('mean scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
                arrangefigures([1,2]); 
            end
        end
        
        function similarity_plot(obj)            
            win                         = ceil(1./median(diff(obj.timescale.global_timescale)));
            comb                        = nchoosek(1:size(obj.binned_data.median_traces,2),2);
            corr_results                = {};
            for pair = 1:size(comb,1)
                corr_results{pair}      = movcorr(obj.binned_data.median_traces(:,comb(pair, 1)),obj.binned_data.median_traces(:,comb(pair, 2)),[win, 0]);
            end

            corr_results                = cell2mat(corr_results);
            precision                   = 1./nanvar(corr_results,[], 2);
            out                         = nanmax(precision(:)) / 10;
            precision(precision > out)  = NaN;

            if obj.rendering
                %figure();hist(nanmean(corr_results, 2),-1:0.01:1,'r')
                f       = figure(1006);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                f.Tag   = 'Similarity plot'; %for figure saving
                ax1     = subplot(3,1,1); hold on;plot(obj.binned_data.median_traces);title('Cell Signal');ylim([-50,range(obj.binned_data.median_traces(:))]);
                ax2     = subplot(3,1,2); hold on;plot(corr_results,'Color',[0.9,0.9,0.9]);
                hold on;plot(nanmean(corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(corr_results, 1)],[-1,-1],'k--');plot([0,size(corr_results, 1)],[1,1],'k--')
                ax3     = subplot(3,1,3); hold on;plot(precision);title('Precision (1/Var)');set(ax3, 'YScale', 'log');plot([0,size(precision, 1)],[1,1],'k--')
                linkaxes([ax1,ax2, ax3],'x')
                arrangefigures([1,2]); 
            end
        end
    
        function norm_cumsum = plot_events_distribution(obj)
            
            peaks = obj.event_fitting.post_correction_peaks;
            
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
                %figure();plot(global_timescale,obj.binned_data.median_traces(:,gp));hold on;scatter(obj.event_fitting.peak_times,obj.event_fitting.post_correction_peaks(:,gp), 'kv')
                if ~isempty(peaks)
                    histogram(peaks(:,gp),0:bin_size:max_peaks, 'FaceAlpha', 0.8,'EdgeColor','none', 'FaceColor', cmap(gp, :));hold on;
                end
                title(obj.binned_data.bin_legend(gp));
            end

            %% Plot all individual peaks to see if they covary
            figure(1005);cla();plot(peaks);legend(obj.binned_data.bin_legend);title('peaks values per subgroups'); xlabel('event #'); ylabel('Amplitude');set(gcf,'Color','w');
            
            %% Plot the mean amplitude per subgroup of peaks, per distance          
            %bin_step = 10;
            %norm_cumsum = cumsum(obj.binned_data.median_traces) ./ nanmax(cumsum(obj.binned_data.median_traces));
            norm_cumsum = cumsum(peaks) ./ nanmax(cumsum(peaks));
            figure(1007);cla();plot(norm_cumsum); hold on;set(gcf,'Color','w');xlabel('Event #');ylabel('normalized cumulative amplitude')
            title('cumulative sum of peaks'); ylim([0, 1]);legend(obj.binned_data.bin_legend,'Location','southeast');
        end
        
        function assess_variability(obj)            
            %sr = nanmedian(diff(obj.timescale{expe}.global_timescale));
            vmr = nanvar(obj.event_fitting.post_correction_peaks,[],2)./nanmean(obj.event_fitting.post_correction_peaks, 2);
            %cv  = nanstd(obj.event_fitting.post_correction_peaks,[],2)./nanmean(obj.event_fitting.post_correction_peaks, 2); % (maybe chack snr at one point?  mu / sigma)
            %fano = []; % windowed VMR. usually for spike trains
            [~, idx] = sort(vmr,'descend');
            
            %figure(123); cla();ylim([0, nanmax(obj.event_fitting.post_correction_peaks(:))]); hold on;
            %     for event = idx'
            %         show_event(obj.binned_data.median_traces, round(obj.event_fitting.peak_times/sr), event);
            %         drawnow;%pause(0.1)
            %     end

            figure(1009);cla();plot(obj.event_fitting.post_correction_peaks(idx, :)); title('Events sorted by Index of dispersion'); ylabel('Amplitude'); xlabel('event #');set(gcf,'Color','w')

            R = max(range(obj.binned_data.median_traces));
            %norm_vmr = vmr/range(vmr);
            norm_vmr = vmr/mean(vmr);
            figure(1010);clf();
            ax1 = subplot(2,1,1);plot(obj.timescale.global_timescale, obj.binned_data.median_traces); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
            ax2 = subplot(2,1,2);plot(obj.event_fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
            plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
            plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
            hold on;plot([obj.timescale.global_timescale(1), obj.timescale.global_timescale(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');
            linkaxes([ax1, ax2], 'x');

            %figure();histogram(norm_vmr, round(10*range(norm_vmr)/std(norm_vmr)))

            figure(1011);cla();scatter(nanmedian(obj.event_fitting.post_correction_peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
        end

        
        function plot_cc(obj)            
            if obj.rendering
                figure(1008);cla();imagesc(obj.crosscorr); hold on;set(gcf,'Color','w');
                caxis([0,1]); hold on;
                xticks([1:size(obj.crosscorr, 1)])
                yticks([1:size(obj.crosscorr, 1)])
                colorbar; hold on;
                plot([1.5,1.5],[0.5,size(obj.crosscorr, 1)+0.5],'k-');xticklabels(['whole cell',obj.binned_data.bin_legend]);xtickangle(45);
                plot([0.5,size(obj.crosscorr, 1)+0.5],[1.5,1.5],'k-');yticklabels(['whole cell',obj.binned_data.bin_legend]);
                title('peak amplitude correlation per subgroup');
                arrangefigures([1,2]); 
            end
        end
        
        function f_handle = get.default_handle(obj)
            f_handle                = @(x) load_several_experiments(x, cellfun(@(x) x.data_folder, obj.arboreal_scans, 'UniformOutput', false));
        end
        
        function crosscorr = get.crosscorr(obj)
            if isempty(obj.event_fitting)
                %error_box('Events were not extracted. Returning CC of traces', 1)
                mask        = ~isnan(sum(obj.binned_data.median_traces,2))    ; 
                crosscorr   = corrcoef([nanmean(obj.binned_data.median_traces(mask,:), 2), obj.binned_data.median_traces(mask,:)]);
            else
                crosscorr   = corrcoef([nanmean(obj.event_fitting.post_correction_peaks, 2), obj.event_fitting.post_correction_peaks]);
            end
            
        end
        
        function [tree, soma_location, tree_values, mean_bin_cc] = plot_corr_tree(obj) 
            
            %% Identify valid set of traces
            valid_gp            = find(~all(isnan(obj.binned_data.median_traces))); % You get NaN'ed bins if the soma location is not scanned (eg a big pyramidal cell)
            
            %% Reload CC
            cc                  = obj.crosscorr(2:end,2:end);
            
            %% Build tree values per bin
            mean_bin_cc     = [];
            ROIs_list       = [];
            if ~isempty(cc)
                for gp = 1:numel(obj.binned_data.groups)
                    roi_of_gp           = obj.binned_data.groups{gp};
                    v_of_gp             = cc(valid_gp(1),gp);
                    ROIs_list           = [ROIs_list, roi_of_gp];
                    mean_bin_cc         = [mean_bin_cc, repmat(v_of_gp, 1, numel(roi_of_gp))];
                end
            end

            %% Map CC values on the tree
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(mean_bin_cc, ROIs_list, obj.default_handle, 'Correlation with most proximal segment');
            caxis([0,1]);
            col = colorbar; col.Label.String = 'Correlation coeff of peaks amplitude with soma';
        end
        
        function get_dimensionality(obj, cross_validate)
            if nargin < 2 || isempty(cross_validate)
                cross_validate = false;
            end
            
            rescaled_traces     = obj.get_rescaled_traces();
            all_ROIs            = 1:size(rescaled_traces, 2);
            normal_n_NaN        = median(sum(isnan(rescaled_traces))) * 4;
            valid_trace_idx     = sum(isnan(rescaled_traces)) <= normal_n_NaN; % eclude traces with too many NaNs (eg. traces that got masked completely)
            rescaled_traces     = fillmissing(rescaled_traces(:, valid_trace_idx),'spline'); % removed funny traces
            %valid_ROIs          = all_ROIs(valid_trace_idx);

            %% Get single or multiple factor estimate
            if ~cross_validate
                [LoadingsPM, specVarPM, T, stats, F] = factoran(double(rescaled_traces), 10);
            else
                %% Cross validation (thanks Harsha)
                [~,N]              = size(rescaled_traces);
                nFactors                    = 40;%round(0.66*N);
                nCV                         = 5;
                train                       = NaN(nFactors, nCV);
                test                        = NaN(nFactors, nCV);
                fac_steps                   = 1;

                for jj = 1:nCV
                    %% Partition data into Xtrain and Xtest
                    test_idx = randperm(size(rescaled_traces, 1));
                    Xtrain = rescaled_traces(test_idx(1:2:end), :);   
                    Xtest  = rescaled_traces(test_idx(2:2:end), :);

                    [jj, toc]
                    parfor factor_nb = 1:nFactors  
                        if ~rem(factor_nb-1, fac_steps) % every fac_steps steps, starting at 1
                            [LoadingsPM, specVarPM, ~, stats, F]    = factoran(Xtrain, factor_nb);    
                            train(factor_nb, jj)                 = stats.loglike;
                            test(factor_nb, jj)                  = testloglike_factorAnalysis(Xtest, nanmean(Xtrain,1), LoadingsPM, specVarPM);
                        end
                    end
                end
                toc

                tested_modes = 1:fac_steps:nFactors;
                figure();plot(tested_modes,train(tested_modes,1:nCV),'o');hold on;plot(tested_modes,nanmean(train(tested_modes,1:nCV), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('train');
                figure();plot(tested_modes,test(tested_modes,1:nCV),'o');hold on;plot(tested_modes,nanmean(test(tested_modes,1:nCV), 2),'ko-')
                xlabel('number of modes');ylabel('?');set(gcf,'Color','w');title('test');

                error('Add optimal selection here')
            end

            %% Store results
            obj.dimensionality.LoadingsPM         = LoadingsPM;
            obj.dimensionality.specVarPM          = specVarPM;
            obj.dimensionality.T                  = T;                % Rotation matrix
            obj.dimensionality.stats              = stats;            % Factoran stats
            obj.dimensionality.F                  = F;                % components
            obj.dimensionality.all_ROIs           = all_ROIs;         % first occurences
            obj.dimensionality.valid_trace_idx    = valid_trace_idx;  % additional filter for recordings with many NaNs
        end
        
        function [tree, soma_location, tree_values, values] = plot_dim_tree(obj, comp)
            if nargin < 2 || isempty(comp) || ~comp
                [tree, soma_location, tree_values, values] = obj.plot_strongest_comp_tree();
                return
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
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, titl); 
        end
        
        function plot_weight_map(obj)
            
            %% Recover traces
            rescaled_traces = obj.get_rescaled_traces();
            
            %% Reload loadings
            LoadingsPM = obj.dimensionality.LoadingsPM;
            Valid_ROIs = obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx);
            rescaled_traces = rescaled_traces(:, Valid_ROIs);

            
            figure(1017);cla();imagesc(LoadingsPM(:,1:5)); colorbar;set(gcf,'Color','w');title('Components 1 to 5 per ROI');xlabel('component');ylabel('ROI');
            [~, loc] = max(LoadingsPM(:,1:5)');

            w1 = LoadingsPM(:,1)/sum(LoadingsPM(:,1));
            w2 = LoadingsPM(:,2)/sum(LoadingsPM(:,2));
            w3 = LoadingsPM(:,3)/sum(LoadingsPM(:,3));
            w4 = LoadingsPM(:,4)/sum(LoadingsPM(:,4));
            w5 = LoadingsPM(:,5)/sum(LoadingsPM(:,5));
            figure(1021);clf();
            for w = {w1, w2, w3, w4, w5}
                hold on; plot(nanmean(rescaled_traces'.* w{1}, 1));
            end
            legend();set(gcf,'Color','w');
            title('Weighted signal average per component')
        end
        
        
        function [tree, soma_location, tree_values, values] = plot_seg_length_tree(obj)   
            
            %% Map dimension weights on the tree
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(Pvec_tree(obj.ref.simplified_tree{1}), '', '', 'Segment length'); 
        end 
        
        
        function [tree, soma_location, tree_values, values] = plot_dist_tree(obj, bin)   
            if nargin < 2 || isempty(bin)
                bin = [];
            end
            
            %% Get the ROIs that were not filtered out 
            ROIs_list               = obj.ref.indices.ROIs_list;
            
            %% Get distance of each segment from soma
            fixed_tree              = obj.ref.simplified_tree_filtered;
            values                  = Pvec_tree(fixed_tree{1}); % distance from soma of each ROIs for non-excluded branches    
            
            %% Map dimension weights on the tree
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, ROIs_list, '', 'Distance from soma', fixed_tree); 
        end 
        

        function [tree, soma_location, tree_values, values] = plot_strongest_comp_tree(obj, n_dim)
            if nargin < 2 || isempty(n_dim)
                n_dim = size(obj.dimensionality.LoadingsPM, 2);
            end
            
            [~, loc]    = nanmax(obj.dimensionality.LoadingsPM(:,1:n_dim),[],2);            
            Valid_ROIs  = find(obj.dimensionality.valid_trace_idx);
            values      = NaN(size(Valid_ROIs));
            for ROI = 1:numel(Valid_ROIs)
                values(ROI)     = loc(ROI);
            end       
            values = values(~isnan(values));
            
            %% Map dimension weights on the tree
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, 'Location of strongest component');       
        end
        
        function save_figures(obj)
            error('to_update')
            p = get(groot,'DefaultFigurePosition');
            folder = [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/'];
            if isfolder(folder)
                rmdir(folder,'s');
            end
            mkdir(folder);
            for fig_idx = [1001:1021, 1081, 1082, (10020+1):(10020+size(obj.event_fitting{expe}.events, 1)), (20020+1):(20020+5)]
                f = figure(fig_idx);
                set(f, 'Position', p)
                try
                    saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                    saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                    savefig(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Children(end).Title.String,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
                catch % for figures with subplots
                    try
                        saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.pdf']);
                        saveas(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.png']);
                        savefig(f, [obj.export_folder, '/',data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'/',f.Tag,' ', data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.fig']);
                    end
                end
            end
            name = parse_paths([folder,'/',f.Tag,data_folders_per_exp{original_expe_idx}(1).name(1:end-9),'.mat']);
            batch_params = obj.general_info{expe}.arboreal_scan.batch_params;

            %% Save every time in case of failure.
            save(parse_paths([obj.export_folder, 'summary.mat']), 'obj','-v7.3')
            %save(name, 'obj.binned_data{expe}.median_traces','all_fits', 'all_pks', 'peak_times','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','-v7.3');   
            %save(name, 'all_traces_per_bin','all_traces_per_bin','all_traces_concat','all_fits', 'all_pks', 'peak_times','cc','all_pks','peak_times','trial_timescale','global_timescale','bin_legend','batch_params','ROIs_list','dimensionality','-v7.3');   
        end



        function process(obj, condition, filter_win, rendering)
            if nargin < 2 || isempty(condition)
                condition = {'distance',100};                
            end  
            if nargin < 3 || isempty(filter_win)
                filter_win = [0,0];                
            end            
            if nargin >= 4 && ~isempty(rendering)
                obj.rendering = rendering;
            end
            save_data       = false;            
            obj.filter_win  = filter_win;

            %% Load and concatenate traces for the selected experiment
            %obj.load_extracted_data();   % this also sets the current expe #
% 
%             %% Check for consistent number of ROIs
%             quality_control(1, obj.extracted_traces)
%             %quality_control(2, cell2mat(all_sr))

            %% Prepare binning of ROIs based on specific grouping condition
            obj.prepare_binning(condition);

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces();

            %% Create median trace per bins
            obj.set_median_traces()

            %% Get summary covariance plots
            obj.similarity_plot();            

            %% Correct for decay to avoid overstimating peak amplitude
            if isempty(obj.event_fitting)
                obj.event_fitting = detect_and_fit_events(obj.binned_data.median_traces, obj.timescale.global_timescale, obj.demo, obj.binned_data.bin_legend);arrangefigures([1,2]);
            end

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.plot_events_distribution();

            %% Study bAP heterogeneity
            obj.assess_variability()

            %% Check how correlated are the peaks between different parts of the tree
            if isempty(obj.crosscorr)
                obj.crosscorr = corrcoef([nanmean(obj.event_fitting.post_correction_peaks, 2), obj.event_fitting.post_correction_peaks]);
            end
            obj.plot_cc(); 

            %% Plot a tree correlation metric     
            %ROIs_list  = unique([all_ROI_ids_per_bin{:}]);
            %[fixed_tree, soma_location, ROIs_list, listing, batch_params] = rebuild_tree(data_folders_per_exp{original_expe_idx}(1), setting_file_path, true);

            %% Assign value of group to these ROIs
            obj.plot_corr_tree(); 

            cross_validate = false;
            obj.get_dimensionality(cross_validate);

            %% Plot weight-tree for each component
            for comp = 1:5
                obj.plot_dim_tree(comp);
            end

            %% Plot a map of tree weights by ROI number for each component
            obj.plot_weight_map();

            %% Plot strongest componenet per ROI
            obj.plot_strongest_comp_tree();

            %% Optionally, if external variables need an update
            %obj.update_external_metrics(60)
            
            figure(999);plot(obj.binned_data.median_traces - nanmean(obj.binned_data.median_traces,2));%set(0,'DefaultAxesColorOrder',magma(size(obj.binned_data.median_traces,2)))

            %% Save figure
            if save_data
                obj.save_figures()
            end
            %close all
        end
        
        function save(obj, auto)
            if nargin < 2 || isempty(auto) || ~auto
                save_folder = uigetdir(obj.source_folder);
            else
                save_folder = obj.source_folder;
            end
            if save_folder
                name = strsplit(obj.source_folder, '/');
                name = name{end-1};
                save([obj.source_folder, name],'obj','-v7.3')
            end
        end
    end
end