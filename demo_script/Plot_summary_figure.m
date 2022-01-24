%% Script showing 2 ways of displaying traces and behaviour and tree structure, either from the raw data or from an extracted arboreal can 
demo_case = 'from_raw'

switch demo_case
    case 'from_raw'
        %% #### Case 1
        % load data directly from data_folders

        data_folder = c.current_data_folder;%'D:\Curated Data\2019-10-04\experiment_1\11-40-09';

        base_extraction_settings = analysis_params('source',data_folder,'signal_channel',2,'smoothing',[1, 0],'z_score',-5,'hi_res',true,'norm_method','dF/F0');
        base_extraction_settings.branch_ids = [  1   2   3   4   5   6   7   8   9  13  14  15  16  17  18  19  20  21  22  23  24];%  25  26  27  28  29  30  31  32  33  34  35  36 ]      
        %% Load Experiment Image
        data_image = load_experiment(apply_default_posthoc_mc(base_extraction_settings),'rendering_mode','realistic','fill_gaps',true,'patch_expansion_factor',5,'tracing_source','swc');
        projection = correlation_image(squeeze(data_image{1}(:,:,:,:,1)));
        
        %% Load SWC Data
        [data_tree,~,~,info] = load_experiment(base_extraction_settings,'compression','0D','tracing_source','swc');

        %% Load pop Data
        base_extraction_settings.branch_ids = 0;
        [data_pop]  = load_experiment(base_extraction_settings,'compression','0D', 'tracing_source','marker');
    
        
          %% Detect spines signals from data_image
        [spine_signals, spine_locations, events] = get_spine(data_image{1}, 0, 3, [], false);
        advanced_spine_analysis()
      
        
        %% Get behavioural data
        names = fieldnames(info.external_var);
        whisker = info.external_var.(names{find(contains(names,'_whisker'),1,'first')}).value;
        whisker = (whisker' - nanmin(whisker)) / range(whisker);
        encoder = interpolate_to(info.external_var.encoder.value,size(data_pop, 4));
        encoder = (encoder - nanmin(encoder)) / range(encoder);
        
        %% Get tree morpho
        tree = arboreal_scan(data_folder);
        tree.prepare_extraction(data_folder);
        tree.prepare_tree();
        [tree_only,pop_only]  = tree.split_tree_and_pop();
        [simp_tree] = tree.split_tree_and_pop(false, true);
        
        %% Prepare figures and plot traces
        fig = figure(123);clf();fig.Color = 'w';
        
        %% Plot ref image
        subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
        ax = subplot(3,3,[2,5]);cla();imagesc(projection);axis image;axis off;colormap(ax, Greens);
        
        %% Plot whole cell and pop Ca2+ data
        tr = plot_many_traces(data_tree(:,:,:,:,1,:), subplot(3,3,3), 'k');set(gca, 'XTickLabel', [],'XTick',[]);axis off;text(1, max(tr(:)),'\itTree signal')
        tr = plot_many_traces(data_pop(:,:,:,:,1,:), subplot(3,3,6), 'r');set(gca, 'XTickLabel', [],'XTick',[]);axis off;text(1, max(tr(:)),'\itPopulation signal')
        
        %% Plot Spin Ca2+ data
        tr = plot_many_traces(spine_signals_subset', subplot(3,3,8), 'k');set(gca, 'XTickLabel', [],'XTick',[]);axis off;text(1, max(tr(:)),'\itSelected Spines')

        %% Plot encoder and whisking
        var_of_interrest = [encoder; whisker + 1]';                    
        subplot(3,3,9);cla();plot(var_of_interrest, 'b');set(gca, 'XTickLabel', [],'XTick',[]);axis off;text(1, max(tr(:)),'\itBehavior')

        %% Plot tree shape
        subplot(3,3,[1,4]);cla(); tr = plot_curved_tree(simp_tree,'r');axis equal;axis off;hold on; %[simp_tree, pop_only]
        adjust_tree_plot('Tree morphology', tr, 2); hold on ; plot_curved_tree(pop_only{1},'','','','','','k')
    case 'from object'    
        %% #### Case 2
        % plot from an arboreal_Scan
        as = meta.experiments(4).arboreal_scans{5};

        %% Prepare figures and plot traces
        fig = figure(123);clf();
        fig.Color = 'w';
        subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
        %subplot(3,3,[2,5]);imagesc(correlation_image(squeeze(data1{1}(:,:,:,:,1))));axis image;axis off;colormap(hot);
        [traces, fig]   = plot_many_traces(as.simple_data, subplot(3,3,3), 'k');set(gca, 'XTickLabel', [],'XTick',[]);axis off
        traces          = plot_many_traces(as.simple_pop_data, subplot(3,3,6), 'r');set(gca, 'XTickLabel', [],'XTick',[]);axis off;

        %% Plot encoder and whisking
        beh = as.external_var;
        names = fieldnames(beh);
        whisker = beh.(names{find(contains(names,'_whisker'),1,'first')}).value;
        whisker = (whisker' - nanmin(whisker)) / range(whisker);
        encoder = interpolate_to(beh.encoder.value,numel(as.timescale));
        encoder = (encoder - nanmin(encoder)) / range(encoder);

        var_of_interrest = [encoder; whisker + 1]';                    
        subplot(3,3,9);plot(var_of_interrest, 'b');set(gca, 'XTickLabel', [],'XTick',[]);axis off;
        subplot(3,3,[1,4]); plot_tree(as.original_tree);axis equal;axis off;
end