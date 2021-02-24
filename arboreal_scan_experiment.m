classdef arboreal_scan_experiment < handle    
    properties
        source_folder
        need_update
        extracted_data_paths % i.e. data_folders 
        arboreal_scans
        ref % first extracted arboreal scan, for conveniency
        t % handl to obj.t

        rendering = true;
        demo = 0;
        
        filter_win      = 0;
        filter_type     = 'gaussian';
        general_info
        peak_thr        = 10;
        
        extracted_traces
        rescaled_traces
        current_segmentation;
        binned_data
        timescale
        event_fitting
        variability
        crosscorr
        dimensionality
        external_variables
        spiketrains
        behaviours        
        default_handle
    end
    
    methods
        function obj = arboreal_scan_experiment(source_folder)
            %% Fix paths
            obj.source_folder       = parse_paths(source_folder);

            %% Get a list of extracted arboreal scans in the source folder
            obj.extracted_data_paths= list_sources(obj);
            obj.need_update         = true(1, numel(obj.extracted_data_paths)); % you're building the object, so they all need an update
            
            %% Load arboreal scans
            obj.update();
        end
        
        function extracted_data_paths = list_sources(obj)
            %% List available arboreal scans
            all_recordings          = dir([obj.source_folder,'/**/*-*-*_exp_*_*-*-*']);
            all_recordings          = all_recordings(~[all_recordings(:).isdir]);
            all_recordings          = all_recordings(~(arrayfun(@(x) strcmp(x.name, '.'), all_recordings) | arrayfun(@(x) strcmp(x.name, '..'), all_recordings)));
            
            names                   = [vertcat(all_recordings.folder), repmat('/',numel(all_recordings), 1), vertcat(all_recordings.name)];
            extracted_data_paths    = cellfun(@(x) parse_paths(x), cellstr(names)', 'UniformOutput', false);
        end
        
        function update(obj)
            %% Update arboreal_scan objects. Identify recordings that changed
            obj.extracted_data_paths = cellfun(@(x) parse_paths(x), obj.extracted_data_paths, 'UniformOutput', false); %temporary
            new_list = obj.list_sources();
            old_list_clean = erase(obj.extracted_data_paths, parse_paths(fileparts(fileparts(obj.extracted_data_paths{1}))));
            new_list_clean = erase(new_list, parse_paths(obj.source_folder));
            added          = ~ismember(new_list_clean, old_list_clean);
            deleted        = ~ismember(old_list_clean, new_list_clean);
            
            %% Clear deleted trees
            obj.arboreal_scans(deleted)         = [];
            obj.need_update(deleted)            = [];
            obj.extracted_data_paths(deleted)   = [];
            
            %% Update required trees first (bc idx are from before the detection of new trees)
            for el = find(obj.need_update)
                add_tree(el);
            end
            
            %% Add new trees            
            for el = new_list(added)
                obj.extracted_data_paths{end+1} = el{1};
                add_tree(numel(obj.extracted_data_paths));                
            end
            
            %% Sort trees by name
            [obj.extracted_data_paths, new_order] = sort(obj.extracted_data_paths);
            obj.arboreal_scans = obj.arboreal_scans(new_order);

            function add_tree(pos)
                obj.arboreal_scans{pos}                         = load(obj.extracted_data_paths{pos});
                obj.arboreal_scans{pos}                         = obj.arboreal_scans{pos}.obj;
                obj.arboreal_scans{pos}.analysis_params.data    = []; % clear full data. If you change t you'll need to rextract everything anyway
                obj.need_update(pos)                            = false;
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
            timescale.durations = timescale.sr.*timescale.tp; % same as cellfun(@(x) x.analysis_params.duration, obj.arboreal_scans)
            timescale.rec_timescale     = arrayfun(@(x, y) linspace(0, y*x, x), timescale.tp, timescale.sr, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) diff(x), timescale.rec_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cellfun(@(x) [x(1), x], timescale.global_timescale, 'UniformOutput', false);
            timescale.global_timescale  = cumsum(horzcat(timescale.global_timescale{:}));
            
            timescale.t_start_nogap = cellfun(@(x) x(end), timescale.rec_timescale); 
            timescale.t_start_nogap = cumsum([0, timescale.t_start_nogap(1:end-1) + timescale.sr(1:end-1)]); 
            timescale.t_start_real  = [];  % to do            
        end
        
        function t = get.t(obj)
            t = obj.timescale.global_timescale;
        end
        
        function ref = get.ref(obj)            
            ref = obj.arboreal_scans{1};
        end
        
        function f_handle = get.default_handle(obj)
            f_handle = @(x) load_several_experiments(x, cellfun(@(x) x.data_folder, obj.arboreal_scans, 'UniformOutput', false));
        end
        
        %% ##############################
        
        function behaviours = get.behaviours(obj)
            behaviours = obj.behaviours();
            behaviours.types            = fieldnames(obj.ref.analysis_params.external_var);  
            behaviours.valid_encoder    = cellfun(@(x) ~isempty(x.analysis_params.external_var.encoder.time), obj.arboreal_scans);
            behaviours.valid_mc_log     = cellfun(@(x) ~isempty(x.analysis_params.external_var.MC.time), obj.arboreal_scans);
            behaviours.valid_behaviours = cell2mat(cellfun(@(y) structfun(@(x) ~isempty(x.time), y.analysis_params.external_var), obj.arboreal_scans, 'UniformOutput', false));
        end
        
        function [raw_beh, downsampd_beh, concat_downsamp_beh] = get_behaviours(obj, type, rendering)
            if nargin < 2 || isempty(type)
                type = 'encoder';
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end
            
            raw_beh                     = cellfun(@(x) x.analysis_params.external_var.(type), obj.arboreal_scans, 'UniformOutput', false);
            downsampd_beh               = {};
            concat_downsamp_beh.time    = [];
            concat_downsamp_beh.value   = [];
            for rec = 1:numel(raw_beh)
                downsampd_beh{rec}.time     = interpolate_to(raw_beh{rec}.time, obj.timescale.tp(rec));
                downsampd_beh{rec}.value    = interpolate_to(raw_beh{rec}.speed, obj.timescale.tp(rec));
                if isempty(downsampd_beh{rec}.time)
                    downsampd_beh{rec}.time = linspace(0, obj.timescale.durations(rec), obj.timescale.tp(rec));
                    downsampd_beh{rec}.value = NaN(1,obj.timescale.tp(rec));
                end
                concat_downsamp_beh.time = [concat_downsamp_beh.time, downsampd_beh{rec}.time + obj.timescale.t_start_nogap(rec)];
                concat_downsamp_beh.value = [concat_downsamp_beh.value, downsampd_beh{rec}.value];
            end 
            
            if rendering
                figure(1026);cla();
                for rec = 1:numel(raw_beh)
                    if ~isempty(downsampd_beh{rec}.time)
                        plot(downsampd_beh{rec}.time + obj.timescale.t_start_nogap(rec), downsampd_beh{rec}.value);hold on;
                    end
                end
            end                
        end     
        
        function [bouts,beh_sm] = get_activity_bout(obj, type, rendering)
            if nargin < 2 || isempty(type)
                type = 'encoder';
            end
            if nargin < 3 || isempty(rendering)
                rendering = false;
            end
            [~, ~, beh] = obj.get_behaviours(type, false);
            beh = beh.value;
            beh_sm = smoothdata(beh, 'gaussian', [50, 0]);
            %beh_sm = detrend(fillmissing(beh_sm,'nearest'),'linear',cumsum(obj.timescale.tp));
            %thr = prctile(beh_sm(beh_sm > 0), 20);
            thr = prctile(beh_sm,20) + range(beh_sm)/20; % 5% of max

            dt = nanmedian(diff(obj.t));
            margin_duration = 3; %in sec
            if numel(margin_duration) == 1
                margin_duration = [margin_duration, margin_duration];
            end            
            
            %% Define bouts
            active_tp   = abs(beh_sm) > abs(thr);
            [starts, stops] = get_limits(active_tp, beh_sm);
            
            %% Add some pts before and after each epoch
            for epoch = starts
                active_tp(max(1, epoch-round(margin_duration(1)*(1/dt))):epoch) = 1;
            end
            for epoch = stops
                active_tp(epoch:min(numel(active_tp), epoch+round(margin_duration(2)*(1/dt)))) = 1;
            end
            [starts, stops] = get_limits(active_tp, beh_sm); %update bouts edges now that w extended the range
            
            if rendering
                figure(1027);plot(beh_sm);hold on;            
                for el = 1:numel(starts)
                    x = [starts(el),starts(el),stops(el),stops(el)];
                    y = [0,nanmax(beh_sm),nanmax(beh_sm),0];
                    patch('XData',x,'YData',y,'FaceColor','red','EdgeColor','none','FaceAlpha',.1);hold on
                end
            end
            
            bouts = sort([starts, stops]);
            
            function [starts, stops] = get_limits(active_tp, beh_sm)
                starts      = find(diff(active_tp) == 1);
                stops       = find(diff(active_tp) == -1);
                if any(starts)
                    if starts(1) > stops(1)
                        starts = [1, starts];
                    end
                    if stops(end) < starts(end)
                        stops = [stops, numel(beh_sm)];
                    end
                elseif all(active_tp)
                    starts = 1; stops = numel(beh_sm);
                end
            end            
        end
        
        %% #############################
        
        function [tree, soma_location, tree_values, values] = plot_dist_tree(obj, bin)  
            if nargin < 2 || isempty(bin)
                bin = [];
            end
            [tree, soma_location, tree_values, values] = obj.ref.plot_dist_tree(bin); 
        end
        
        %         function [tree, soma_location, tree_values, values] = plot_seg_length_tree(obj)
        %             %% Map dimension weights on the tree
        %             [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(Pvec_tree(obj.ref.simplified_tree{1}), '', '', 'Segment length'); 
        %         end 
        
        %% #############################

        function prepare_binning(obj, condition)
%             if isempty(obj.current_segmentation) || isempty(obj.binned_data) || isempty(obj.current_segmentation.condition) || (~strcmp(obj.current_segmentation.condition{1}, condition{1}) || obj.current_segmentation.condition{2} ~= condition{2})
            	obj.current_segmentation            = {};
                obj.current_segmentation.condition  = condition;
%             elseif strcmp(obj.current_segmentation.condition{1}, condition{1}) && obj.current_segmentation.condition{2} == condition{2}
%                 return
%             e
            
            %% Define current binning rule. See arboreal_scan.get_ROI_groups for more info % CC and peak extractions are based on this binning    
            [obj.binned_data.groups, obj.binned_data.metrics, obj.binned_data.bin_legend] = obj.ref.get_ROI_groups(obj.current_segmentation.condition, obj.demo);
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
            traces                        = obj.extracted_traces;
            invalid = ~ismember(1:80, obj.ref.indices.valid_swc_rois');
            for idx = 1:numel(traces)
                traces{idx}(:,invalid) = NaN;
            end
            
            [best_scal_f, best_offset]    = scale_every_recordings(traces, obj.demo); % qq consider checking and deleting "scale_across_recordings"
            obj.general_info.scaling      = best_scal_f;
            obj.general_info.offset       = best_offset;
            if obj.rendering
                %figure();plot(obj.t,nanmedian(obj.rescaled_traces,2));xlabel('Scan Time (s)') % FF may not be necessary
                arrangefigures([1,2]);
            end
        end
        
        function rescaled_traces = get.rescaled_traces(obj)
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
            rescaled_traces = obj.rescaled_traces;
            
            %% Create median trace per bins
            all_traces_per_bin = cell(1, numel(obj.binned_data.groups));
            for gp = 1:numel(obj.binned_data.groups)
                columns                     = ismember(obj.binned_data.readmap, obj.binned_data.groups{gp});
                all_traces_per_bin{gp}      = nanmedian(rescaled_traces(:,columns), 2);
            end    
            
            obj.binned_data.median_traces = cell2mat(all_traces_per_bin);
            
            if obj.rendering
                %% Plot the mean trace for each bin  
                figure(1001);cla();
                plot(obj.t, obj.binned_data.median_traces); hold on;
                legend(obj.binned_data.bin_legend); hold on;
                title('median scaled trace per group');xlabel('time (s)');set(gcf,'Color','w');
                figure(1023);plot(obj.binned_data.median_traces - nanmean(obj.binned_data.median_traces,2));hold on;title('group traces median - overall median');xlabel('time (s)');set(gcf,'Color','w');
                arrangefigures([1,2]); 
            end
        end
        
        function [precision] = similarity_plot(obj)            
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

                if obj.rendering
                    %figure();hist(nanmean(corr_results, 2),-1:0.01:1,'r')
                    f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                    f.Tag   = 'Similarity plot'; %for figure saving
                    ax1     = subplot(3,1,1); hold on;plot(obj.binned_data.median_traces);title('Cell Signal');ylim([-50,range(obj.binned_data.median_traces(:))]);
                    ax2     = subplot(3,1,2); hold on;plot(corr_results,'Color',[0.9,0.9,0.9]);
                    hold on;plot(nanmean(corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(corr_results, 1)],[-1,-1],'k--');plot([0,size(corr_results, 1)],[1,1],'k--')
                    ax3     = subplot(3,1,3); hold on;plot(precision);title('Precision (1/Var)');set(ax3, 'YScale', 'log');plot([0,size(precision, 1)],[1,1],'k--')
                    linkaxes([ax1,ax2, ax3],'x')
                    arrangefigures([1,2]); 
                end
            else
                f       = figure(1022);clf();title('Similarity plot');hold on;set(gcf,'Color','w');hold on;
                annotation('textbox','String','N/A','FitBoxToText','on','FontSize',16,'FontWeight','bold')
            end
        end
    
        function norm_cumsum = get_events_statistics(obj)
            
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
            obj.variability.index_of_disp = norm_vmr;
            figure(1010);clf();
            ax1 = subplot(2,1,1);plot(obj.t, obj.binned_data.median_traces); ylabel('Amplitude'); hold on;set(gcf,'Color','w');ylim([-R/20,R + R/20]);title('bin traces'); hold on;
            ax2 = subplot(2,1,2);plot(obj.event_fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
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
            %             plot(obj.event_fitting.peak_times, norm_vmr, 'ko-'); title('Index of dispersion per event'); hold on;
            %             plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr), mean(norm_vmr)],'r');hold on;
            %             plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)+std(norm_vmr), mean(norm_vmr)+std(norm_vmr)],'r--');
            %             hold on;plot([obj.t(1), obj.t(end)] ,[mean(norm_vmr)-std(norm_vmr), mean(norm_vmr)-std(norm_vmr)],'r--');  
            %  
            %             
            %             temp = behaviour(round(obj.event_fitting.peak_pos));
            %             hold on;plot(obj.t(obj.event_fitting.peak_pos), temp,'-vr');

            figure(1011);cla();scatter(nanmedian(obj.event_fitting.post_correction_peaks, 2), vmr, 'filled'); title('Index of dispersion vs Amplitude'); xlabel('Amplitude'); ylabel('VMR'); hold on;set(gcf,'Color','w')
        end

        %% #############################################
        
        function get_correlations(obj) 
            obj.crosscorr = corrcoef([nanmean(obj.event_fitting.post_correction_peaks, 2), obj.event_fitting.post_correction_peaks]);

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
                
                %% Project correlation value onto the tree
                obj.plot_corr_tree(); 
            end            
        end

        function crosscorr = get.crosscorr(obj)
            if isempty(obj.event_fitting)
                return
                %error_box('Events were not extracted. Returning CC of traces', 1)
                mask        = ~isnan(sum(obj.binned_data.median_traces,2))    ;
                trace       = obj.binned_data.median_traces - repmat(nanmean(obj.binned_data.median_traces,2),1,size(obj.binned_data.median_traces, 2));
                crosscorr   = corrcoef(trace, trace);
            else
               crosscorr   = corrcoef([nanmean(obj.event_fitting.post_correction_peaks, 2), obj.event_fitting.post_correction_peaks],'Rows','complete'); % ignore NaNs
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
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(mean_bin_cc, ROIs_list, obj.default_handle, 'Correlation with most proximal segment','',1018);
            %caxis([0,1]);
            col = colorbar; col.Label.String = 'Correlation coeff of peaks amplitude with soma';
        end
        
        %% ###################
        
        function get_spike_trains(obj)
            %% Use ML spike to get a spike train estimate
            if ~exist('spk_autocalibration.m','file') || ~exist('fn_getfile.m','file')
                error('You need to download the "bricks" and "ml_sikes" toolboxes, and add thn to the path')
            end
        
            %% Formatting calcium. It seems that signal need to be normalized            
            calcium = num2cell(fillmissing(normalize(obj.binned_data.median_traces',[], 'norm_method','signal/max','percentile',10)','constant',0), 1);
            dt      = median(diff(obj.t));
            first_valid = find(~cellfun(@(x) all(isnan(x)), calcium),1,'first');
            
            AMAX        = 10; % wtf
            p = obj.event_fitting.global_fit;f = p(2) / (p(2)+p(4));estimated_tau = p(3)*f+p(5)*(1-f); 

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

            %% Autoestimat doesn't work on single runs. We run batches of 1000pts, ignore failed estimates and get a consensus estimate
            batch_size  = 1000;
            tauest      = [];
            aest        = [];
            sigmaest    = [];
            for r = 1:(floor(numel(calcium{first_valid})/batch_size) - 1)
                try
                    [tauest(r), aest(r), sigmaest(r), evt, par] = spk_autocalibration(calcium{first_valid}((r*batch_size):((r+1)*batch_size)),pax);
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
            for bin = first_valid:numel(calcium)
                [obj.spiketrains.spike_estimate{bin},obj.spiketrains.fit{bin},obj.spiketrains.drift{bin}] = spk_est(calcium{bin},par);
                obj.spiketrains.spike_estimate{bin} = obj.spiketrains.spike_estimate{bin}*2; %bc of the signal upsampling
            end
            obj.plot_inferred_spikes();
        end
        
        function plot_inferred_spikes(obj)
            calcium = num2cell(normalize(obj.binned_data.median_traces',[], 'norm_method','signal/max','percentile',10)', 1);
            first_valid = find(~cellfun(@(x) all(isnan(x)), calcium),1,'first');
            
            figure(1025);cla();plot(obj.t,calcium{1},'k');hold on;title('Spike inference per bin');xlabel('time(s)')
            colors = [repmat([NaN, NaN, NaN], first_valid-1, 1) ;lines(numel(calcium) - first_valid + 1)];
            for bin = first_valid:numel(calcium)
                plot(obj.t, obj.spiketrains.fit{bin},'Color',colors(bin, :));hold on
                scatter(obj.spiketrains.spike_estimate{bin}, repmat(bin/10, 1, numel(obj.spiketrains.spike_estimate{bin})),'o','MarkerEdgeColor',colors(bin, :));hold on;
            end
        end
        
        %% ###################
        
        function get_dimensionality(obj, cross_validate)
            if nargin < 2 || isempty(cross_validate)
                cross_validate = false;
            end
            
            rescaled_traces     = obj.rescaled_traces();
            all_ROIs            = 1:size(rescaled_traces, 2);
            normal_n_NaN        = median(sum(isnan(rescaled_traces))) * 4;
            valid_trace_idx     = sum(isnan(rescaled_traces)) <= normal_n_NaN; % eclude traces with too many NaNs (eg. traces that got masked completely)
            rescaled_traces     = fillmissing(rescaled_traces(:, valid_trace_idx),'spline'); % removed funny traces

            %% Get single or multiple factor estimate
            if ~cross_validate
                [LoadingsPM, specVarPM, T, stats, F] = factoran(double(rescaled_traces), 5);
            else
                %% Cross validation (thanks Harsha)
                [~,N]              = size(rescaled_traces);
                nFactors                    = 20;%round(0.66*N);
                nCV                         = 5;
                train                       = NaN(nFactors, nCV);
                test                        = NaN(nFactors, nCV);
                fac_steps                   = 1;
                tic
                for jj = 1:nCV
                    %% Partition data into Xtrain and Xtest
                    test_idx = randperm(size(rescaled_traces, 1));
                    Xtrain = double(rescaled_traces(test_idx(1:2:end), :));   
                    Xtest  = double(rescaled_traces(test_idx(2:2:end), :));

                    [jj, toc]
                    parfor factor_nb = 1:nFactors  
                        if ~rem(factor_nb-1, fac_steps) % every fac_steps steps, starting at 1
                            [LoadingsPM, specVarPM, ~, stats, F] = factoran(Xtrain, factor_nb);    
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
            
            
            %% Plot weight-tree for each component
            for comp = 1:5
                obj.plot_dim_tree(comp);
            end

            %% Plot a map of tree weights by ROI number for each component
            obj.get_weight_map();

            %% Plot strongest componenet per ROI
            obj.plot_strongest_comp_tree();            
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
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, titl, '',  10200 + comp); 
        end
        
        function weighted_averages = get_weight_map(obj, n_weigths)
            if nargin < 2 || isempty(n_weigths)
                n_weigths = size(obj.dimensionality.F, 2);
            end
            
            %% Recover traces
            rescaled_traces = obj.rescaled_traces;
            
            %% Reload loadings
            LoadingsPM = obj.dimensionality.LoadingsPM;
            Valid_ROIs = obj.dimensionality.all_ROIs(obj.dimensionality.valid_trace_idx);
            rescaled_traces = rescaled_traces(:, Valid_ROIs);
            [~, loc] = max(LoadingsPM(:,1:n_weigths)');
            
            all_weights         = {};
            weighted_averages   = [];
            for w = 1:n_weigths
                all_weights{w} = LoadingsPM(:,w)/sum(LoadingsPM(:,w));
                weighted_averages(w, :) = nanmean(rescaled_traces'.* all_weights{w}, 1);
            end
            
            if obj.rendering
                figure(1017);cla();imagesc(LoadingsPM(:,1:n_weigths)); colorbar;set(gcf,'Color','w');title(['Components 1 to ',num2str(n_weigths),' per ROI']);xlabel('component');ylabel('ROI');
                figure(1021);clf();
                for w = 1:n_weigths
                    hold on; plot(weighted_averages(w, :));
                end
                legend();set(gcf,'Color','w');
                title('Weighted signal average per component')
            end
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
            [f, tree_values, tree, soma_location] = obj.ref.plot_value_tree(values, Valid_ROIs, obj.default_handle, 'Location of strongest component','',10200);       
        end
        
        %% #############################################

        function update_all_signals(obj, new_method)
            %% Update raw signals changing the compression procedure
            if ~any(strcmp(new_method, {'max', 'mean', 'median','min'}))
                error('Only max, mean and median method are supported')
            else
                for rec = 1:numel(obj.arboreal_scans)                    
                    obj.arboreal_scans{rec}.extraction_method = new_method;
                    if isempty(obj.arboreal_scans{rec}.analysis_params.data)
                        temp = load(obj.extracted_data_paths{rec});
                        obj.arboreal_scans{rec}.analysis_params.data = temp.obj.analysis_params.data;
                        clear temp;
                    end
                        
                    obj.arboreal_scans{rec}.update_segment_signals();
                    obj.arboreal_scans{rec}.analysis_params.data = [];
                    obj.need_update(rec)                    = true;
                    % qq need to handle need pdate and clear exisitng data
                end
                error_box('COMPRESSION MODES WERE UPDATED BUT YOU NEED TO SAVE THE ARBOREAL SCANS TO KEP THIS CHANGE FOR NEXT RELOADING')
            end
        end

        function process(obj, condition, filter_win, rendering)
            %% Call all processing steps
            if nargin < 2 || isempty(condition)
                condition = {'distance',Inf};                
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
            % obj.load_extracted_data();   % this also sets the current expe #
            % quality_control(1, obj.extracted_traces)
            % quality_control(2, cell2mat(all_sr))

            %% Prepare binning of ROIs based on specific grouping condition
            obj.prepare_binning(condition);

            %% Rescale each trace with a unique value across recordings (which would be specific of that region of the tree).
            obj.rescale_traces();

            %% Create median trace per bins
            obj.set_median_traces()

            %% Get summary covariance plots of the raw traces
            obj.similarity_plot();            

            %% Correct for decay to avoid overstimating peak amplitude
            obj.event_fitting = detect_and_fit_events(obj.binned_data.median_traces, obj.t, obj.demo, obj.binned_data.bin_legend, obj.peak_thr);arrangefigures([1,2]);

            %% Detect and display peak histogram distribution (mean and individual groups)
            obj.get_events_statistics();

            %% Study bAP heterogeneity
            obj.assess_variability()

            %% Check how correlated are the peaks between different parts of the tree
            obj.get_correlations();

            %% Dimensionality reduction
            obj.get_dimensionality(false); % set to true for cross validation
            
            %% Extract behaviours
            obj.get_behaviours();
            
            %% Spike inference
            try
                obj.get_spike_trains();
            end

            %% Optionally, if external variables need an update
            %obj.update_external_metrics(60)

            %% Save figure
            if save_data
                obj.save_figures()
            end            
        end
        
        function save(obj, auto)
            %% Save all data in a mat file
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
        
        function save_figures(obj)
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
            % 1019 � no title  individual events?
            % 1020 : strongest component
                % 10201 � 1020n : individual components
            % 1021 : Weighted signal average per component
            % 1022 : variability assessment            
            % to find un-numbered figures look for "% FF"

            p = get(groot,'DefaultFigurePosition');
            folder = parse_paths([obj.source_folder, '/figures/']);
            if isfolder(folder)
                rmdir(folder,'s');
            end
            mkdir(folder);
            [~, tag] = fileparts(fileparts(obj.source_folder));
            n_dim = size(obj.dimensionality.T, 1);
            n_groups = numel(obj.binned_data.groups);            
            for fig_idx = [1001:1022, 10051:(10050 + n_dim), 10021:(10020+n_groups)]
                f = figure(fig_idx);
                set(f, 'Position', p)
                try
                    saveas(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.pdf']);
                    saveas(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.png']);
                    savefig(f, [folder,'/',f.Children(end).Title.String,' ', tag,'.fig']);
                catch % for figures with subplots
                    try
                        saveas(f, [folder,'/',f.Tag,' ', tag,'.pdf']);
                        saveas(f, [folder, '/',tag,'/',f.Tag,' ', tag,'.png']);
                        savefig(f, [folder, '/',tag,'/',f.Tag,' ', tag,'.fig']);
                    end
                end
            end

            %% Save every time in case of failure.
            save(parse_paths([folder, 'summary.mat']), 'obj','-v7.3')
        end
    end
end