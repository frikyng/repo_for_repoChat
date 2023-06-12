%% Suggests an optimal epsilon for clustering based on DBSCAN or HDBSCAN
% 	Given a set of low-dimensional data, this function performs DBSCAN or 
%   HDBSCAN clustering with a range of epsilon values, storing the number 
%   of clusters and noise points for each epsilon. It then suggests an 
%   optimal epsilon based on the data.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[suggested, n_clust] = test_epsilon(obj, Low_Dim_Data, current_signal, 
%                                        all_ROIs, valid_ROIs, n_clust)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	obj(object):
%                                   The main object that contains parameters 
%                                   for clustering such as minimum group size 
%                                   for DBSCAN.
%
% 	Low_Dim_Data(Matrix) - Optional:
%                                   Low dimensional data for clustering.
%
% 	current_signal(Matrix) - Optional:
%                                   Current signals associated with each data point.
%
% 	all_ROIs(Array) - Optional:
%                                   The complete set of Regions of Interest (ROIs).
%
% 	valid_ROIs(Array) - Optional:
%                                   The valid set of Regions of Interest (ROIs).
%
% 	n_clust(Int) - Optional:
%                                   The number of clusters. If not provided or
%                                   empty, it is initialized within the function.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	suggested(Double)
%                                   The suggested optimal epsilon for clustering.
%
% 	n_clust(Int)
%                                   The optimal number of clusters found.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * If the `n_clust` is not provided, the function will look for an epsilon 
%   value where the ratio of noise points is less than 5%.
% * If the `n_clust` is provided, it will find the last epsilon value where 
%   the number of clusters is equal to `n_clust`.
% * In case of no suggested epsilon, the first epsilon in the range is 
%   suggested.
% -------------------------------------------------------------------------
%% Examples:
% * With number of clusters
% 	[suggested, n_clust] = test_epsilon(obj, Low_Dim_Data, current_signal, 
%                                        all_ROIs, valid_ROIs, n_clust);
%
% * Without number of clusters
% 	[suggested, ~] = test_epsilon(obj, Low_Dim_Data, current_signal, 
%                                        all_ROIs, valid_ROIs);
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% paste license here
% -------------------------------------------------------------------------
% Revision Date:
% 	09-06-2023
% -------------------------------------------------------------------------
% See also: 
%	dbscan, hdbscan
%
% TODO : Check the behavior of the function with different types of data
% and different parameter settings. Improve the function's speed and
% efficiency.

function [suggested, n_clust] = test_epsilon(obj, Low_Dim_Data, current_signal, all_ROIs, valid_ROIs, n_clust)
    % Check if n_clust is given as argument, if not, initialize it as empty
    if nargin < 6 || isempty(n_clust)
        n_clust = [];
    end

    %% Define the method for clustering and initialize variables
    method      = 'dbscan';   % The method of clustering, could be 'dbscan' or 'hdbscan'
    rendering   = false;      % Flag for whether to render the results
    % Set test range based on the method chosen
    if strcmp(method, 'dbscan')
        test_range       = logspace(-1,2,150);  % If the method is 'dbscan', test a logarithmic range of values for epsilon
    else
        test_range       = 1:20;                % Otherwise, test a range from 1 to 20
    end

    % Initialize arrays to keep track of the number of clusters (n_gp) and noise points (n_noise_pt) for each epsilon
    n_gp        = [];
    n_noise_pt  = [];

    current_ep  = 1;  % Initialize the current epsilon index
    % Set the minimum cluster size based on a parameter in the 'obj' object
    MIN_CLUSTER_SIZE = obj.dbscan_min_gp_size;% default (Ester et al., 1996), although we may want 2x NDim for High dimesnional data  (Sander et al., 1998)  
    gp          = []; % Initialize the array of unique cluster indices
    tolerance   = 3;  % Set a tolerance level for the number of times the number of clusters can decrease before stopping
    verbose = false;  % Set a flag for verbosity

    %% Loop over the range of epsilon values
    while current_ep < numel(test_range)
        try
            test_value  = test_range(current_ep);  % Get the current epsilon value to test

            %% Get cluster
            if strcmp(method, 'dbscan')  % If the method is 'dbscan'
                if verbose
                    fprintf(['testing epsilon = ',num2str(test_value), '\n']);  % If verbose, print the current epsilon being tested
                end
                % Perform DBSCAN clustering with the current epsilon and minimum cluster size
                cluster_idx = dbscan(Low_Dim_Data , test_value, MIN_CLUSTER_SIZE); % Minpts from (Sander et al., 1998)
            elseif strcmp(method, 'hdbscan')  % If the method is 'hdbscan'
                if verbose
                    fprintf(['testing Min cluster size = ',num2str(test_value), '\n'])  % If verbose, print the current minimum cluster size being tested
                end
                clusterer = HDBSCAN(Low_Dim_Data);  % Create an HDBSCAN object with the low-dimensional data
                % Run HDBSCAN with the current parameters
                clusterer.run_hdbscan(2,1,1,test_value,0,false);
                cluster_idx = clusterer.labels;  % Get the cluster labels
            else
                error('method not implemented')  % If the method is neither 'dbscan' nor 'hdbscan', throw an error
            end

            % Get the unique cluster indices, excluding noise points
            quote("function [suggested, n_clust] = test_epsilon(obj,", "if isempty(suggested)\n        suggested = test_range(1);\n    end\n end")
            prev_gp = gp;  % Store the previous set of unique clusters
            [gp ,~, indices] = unique(cluster_idx(cluster_idx > 0));  % Find the unique clusters, excluding noise points
            % If the number of clusters has decreased and all clusters are equal to 1, decrease the tolerance
            if numel(gp) < numel(prev_gp) && all(gp == 1) 
                tolerance = tolerance - 1;
            end
            % If the number of clusters has decreased and all clusters are equal to 1, and tolerance is 0, break the loop
            if numel(gp) < numel(prev_gp) && all(gp == 1) && ~tolerance
                break
            end

            % Get a color map for the clusters
            colors = viridis(numel(gp));  % Create a color map using the viridis colormap
            colors = colors(randperm(size(colors,1), size(colors,1)),:);  % Shuffle the colors
            colors = colors(indices, :);  % Get the colors corresponding to each cluster

            % Store the number of clusters and noise points for the current epsilon
            n_gp(current_ep) = numel(gp);
            n_noise_pt(current_ep) = sum(cluster_idx <= 0);

            %% Visualization (optional)
            if rendering
                figure(1);clf();title(num2str(test_value));  % Create a new figure and set the title to the current epsilon
                subplot(1,2,1);  % Create the first subplot
                % Plot the noise points
                scatter3(Low_Dim_Data(cluster_idx <= 0,1), Low_Dim_Data(cluster_idx <= 0,2), Low_Dim_Data(cluster_idx <= 0,3), 30, 'MarkerFaceColor' , [0.8,0.8,0.8], 'MarkerEdgeColor' , 'none'); hold on;
                % Plot the clustered points with colors corresponding to their clusters
                scatter3(Low_Dim_Data(cluster_idx > 0,1), Low_Dim_Data(cluster_idx > 0,2), Low_Dim_Data(cluster_idx > 0,3), 30, colors, 'filled'); hold on;
            end

            % Create a blank version of all_ROIs and fill it with the cluster indices, leaving NaN for noise points
            blank_v = NaN(size(all_ROIs));
            blank_v(valid_ROIs) = cluster_idx;
            blank_v(blank_v <= 0) = NaN;
            if obj.use_hd_data
                values = split_values_per_voxel(blank_v, obj.ref.header.res_list(1:obj.ref.indices.n_tree_ROIs,1),find(~all(isnan(current_signal),2)));
            else
                values = blank_v;
            end
            % Create a color map for all ROIs, with colors for valid ROIs and gray for noise points
            updated_colors = repmat([0.8, 0.8, 0.8], numel(all_ROIs),1);
            updated_colors(~(isnan(blank_v)), :) = colors;
            if rendering
                s = subplot(1,2,2);  % Create the second subplot
                obj.ref.plot_value_tree(values, '' ,'','','',s,'classic',updated_colors);  % Plot the ROIs with their colors
                title(num2str(test_value));  % Set the title to the current epsilon
                pause(2)  % Pause for 2 seconds
            end
        end
        current_ep = current_ep + 1;  % Move to the next epsilon
    end

    % If n_clust is not set, find the first epsilon where the ratio of noise points is less than 5%
    if isempty(n_clust)
        max_loc = find((n_noise_pt/size(Low_Dim_Data, 1)) < 0.05, 1, 'first');
        n_clust = test_range(max_loc);
    else
        % If n_clust is set, find the last epsilon where the number of clusters is equal to n_clust
        max_loc = find(n_gp == abs(n_clust), 1, 'last');
        if isempty(max_loc)
            % If there is no such epsilon, find the epsilon with the maximum number of clusters
            [n_clust, max_loc] = max(n_gp);
        end
    end
    % The suggested epsilon is the epsilon at the location found
    suggested = test_range(max_loc);
    % If no epsilon is suggested, suggest the first epsilon in the range
    if isempty(suggested)
        suggested = test_range(1);
    end
 end

