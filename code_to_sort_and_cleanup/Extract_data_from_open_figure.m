% Script to extract X and Y data from an opened figure
plot_number = gcf; % Get the plot number, usually 1
% Get the axes and data objects
axes_objects = get(plot_number, 'Children'); 
data_objects = get(axes_objects, 'Children');
% Get the data from the data objects
xdata = get(data_objects{2}, 'XData');  % for a subplot, use : xdata = get(data_objects{ii}, 'XData');
ydata = get(data_objects{2}, 'YData');  % % for a subplot, use : ydata = get(data_objects{ii}, 'YData');
% If you have a 3D plot, uncomment the line below
% zdata = get(dataObjs, 'ZData');

xdata = cell2mat(xdata);
ydata = cell2mat(ydata);

ydata(1,:) = [];
range_vec_all = ydata';