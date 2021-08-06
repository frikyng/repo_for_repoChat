%% Demo script to show how you could use various sources for 3d coordinates

%% a regular LD-HD set
h = load_header('D:\Curated Data\2019-09-17\experiment_1\11-49-45')
h2 = load_header('D:\Curated Data\2019-09-17\experiment_1\12-39-17')

%% 
% h = load_header('D:\Curated Data\2019-07-26\experiment_4\20-30-04')
% h2 = load_header('D:\Curated Data\2019-07-26\experiment_4\22-33-45')

% LD SET
% h = load_header('D:\Curated Data\2019-09-26\experiment_1\14-57-53');
% h2 = load_header('D:\Curated Data\2019-09-26\experiment_1\16-57-54'); 


%% Overlay full scan and simplifed scan
x_sta = cellfun(@(x) x(1,:), h.start_pixels, 'Uni', false);
x_sta = [x_sta{:}];
y_sta = cellfun(@(x) x(2,:), h.start_pixels, 'Uni', false);
y_sta = [y_sta{:}];
z_sta = cellfun(@(x) x(3,:), h.start_pixels, 'Uni', false);
z_sta = [z_sta{:}];

x_sto = cellfun(@(x) x(1,:), h.stop_pixels, 'Uni', false);
x_sto = [x_sto{:}];
y_sto = cellfun(@(x) x(2,:), h.stop_pixels, 'Uni', false);
y_sto = [y_sto{:}];
z_sto = cellfun(@(x) x(3,:), h.stop_pixels, 'Uni', false);
z_sto = [z_sto{:}];

figure();
ax1 = subplot(1,2,1);hold on;
title('LD FULL SCAN');hold on;
plot3([y_sta;y_sto],512-[x_sta;x_sto],[z_sta;z_sto]/h.xy_to_z_scaling,'Color',[0.8,0.8,0.8]); hold on;
plot_tree(h.trees{1},'r'); hold on


ax2 = subplot(1,2,2);title('HD PARTIAL SCAN');hold on;

x_sta = cell2mat(cellfun(@(x) x(1,:), h2.start_pixels, 'Uni', false));
y_sta = cell2mat(cellfun(@(x) x(2,:), h2.start_pixels, 'Uni', false));
z_sta = cell2mat(cellfun(@(x) x(3,:), h2.start_pixels, 'Uni', false));
x_sto = cell2mat(cellfun(@(x) x(1,:), h2.stop_pixels , 'Uni', false));
y_sto = cell2mat(cellfun(@(x) x(2,:), h2.stop_pixels , 'Uni', false));
z_sto = cell2mat(cellfun(@(x) x(3,:), h2.stop_pixels , 'Uni', false));

plot3([y_sta;y_sto],512-[x_sta;x_sto],[z_sta;z_sto]/h2.xy_to_z_scaling,'Color',[0.8,0.8,0.8]); hold on;

%% Add full simplified neuron tree from HD header
for br = 2:numel(h2.trees)
    plot_tree(h2.trees{br},[1,0.8,0.8]); hold on
end
plot_tree(h2.trees{1},'r'); hold on

%% Add full original tree from HD header
plot_tree(h2.original_trees{1},'b'); hold on

%% TO ADD
% 
% %% Add original tree from tree object
% plot_tree(tree.original_tree{1}, 'r'); hold on
% 
% %% Add original tree minus excluded branches
% plot_tree(tree.original_tree_filtered{1}, 'k'); hold on
% 
% %% Add original population
% plot_curved_tree(tree.original_pop{1}); hold on
% 
% %% Add simplified tree from tree object
% plot_tree(tree.simplified_tree{1}, 'm'); hold on
% 
% %% Add simplified tree filtered from tree object
% plot_tree(tree.simplified_tree_filtered{1}, 'b'); hold on

linkprop([ax1,ax2],{'CameraPosition','CameraUpVector'})


%% watchout for         xy_to_z_scaling: 1.4597
%                      pxl_per_um: 1.4629
%                      resampling: 1
%                    z_resampling: 1





