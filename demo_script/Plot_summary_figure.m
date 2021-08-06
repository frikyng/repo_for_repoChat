%% Script showing 2 ways of displaying traces and behaviour and tree structure, either from the raw data or from an extracted arboreal can 

%% #### Case 1
% plot directly from data_folders

data_folder = 'D:\Curated Data\2019-10-04\experiment_1\11-40-09';

%% Load Experiment Image
data1 = load_experiment('source',data_folder,'rendering_mode','realistic','signal_channel',2,'smoothing',10,'fill_gaps',true,'hi_res',true,'z_score',20);

%% Load SWC Data
[data2,~,~,info] = load_experiment('source',data_folder,'signal_channel',2,'smoothing',10,'compression','0D','norm_method','dF/F0', 'tracing_source','swc','z_score',20);

%% Load tree Data
[data3]  = load_experiment('source',data_folder,'signal_channel',2,'smoothing',10,'compression','0D','norm_method','dF/F0', 'tracing_source','marker','z_score',20);

%% Prepare figures and plot traces
fig = figure(123);clf();
fig.Color = 'w';
subplot = @(m,n,p) subtightplot (m, n, p, [0.01 0.05], [0.1 0.01], [0.1 0.01]);
subplot(3,3,[2,5]);imagesc(correlation_image(squeeze(data1{1}(:,:,:,:,1))));axis image;axis off;colormap(hot);
plot_many_traces(data2(:,:,:,:,1,:), subplot(3,3,3), 'k');set(gca, 'XTickLabel', [],'XTick',[]);axis off
plot_many_traces(data3(:,:,:,:,1,:), subplot(3,3,6), 'r');set(gca, 'XTickLabel', [],'XTick',[]);axis off;

%% Plot encoder and whisking
names = fieldnames(info.external_var);
whisker = info.external_var.(names{find(contains(names,'_whisker'),1,'first')}).value;
whisker = (whisker' - nanmin(whisker)) / range(whisker);
encoder = interpolate_to(info.external_var.encoder.value,size(data3, 4));
encoder = (encoder - nanmin(encoder)) / range(encoder);

var_of_interrest = [encoder; whisker + 1]';                    
subplot(3,3,9);plot(var_of_interrest, 'b');set(gca, 'XTickLabel', [],'XTick',[]);axis off;

tree = arboreal_scan(data_folder);
tree.prepare_extraction(data_folder);
tree.prepare_tree();
tree_only = tree.split_tree_and_pop();
subplot(3,3,[1,4]); plot_tree(tree_only);axis equal;axis off;


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