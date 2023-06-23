%% Shows summary predictive scores in a figure
% This function takes in predictive scores in a cell array and presents them
% in a figure. The results are shown in two subplots, each subplot has its 
% results normalized differently.
% -------------------------------------------------------------------------
%% Syntax:
% 	fig_handle = show_summary_predictive_scores(all_results, x_ticks, y_ticks, y_label)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	all_results(cell array): 
%                                   This is a cell array of predictive scores 
%                                   to be displayed.
%
% 	x_ticks(cell array) - Optional:
%                                   This is an optional cell array of labels 
%                                   for the x-axis. If not provided, default
%                                   labels are used.
%
% 	y_ticks(cell array) - Optional:
%                                   This is an optional cell array of labels 
%                                   for the y-axis. If not provided, default
%                                   labels are used.
%
% 	y_label(string) -- Optional:
%                                   This is an optional label for the y-axis.
%                                   If not provided, default label 'Phate #' 
%                                   is used.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	fig_handle(Figure Handle) Handle
%                                   This is a handle to the figure that is 
%                                   created by this function.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * The function takes in all_results as a cell array of numbers and displays
%   them in a figure. The results are normalized and displayed in two subplots.
% * This function also sets the color of the figure to white.
% -------------------------------------------------------------------------
%% Examples:
% * Displaying predictive scores with default labels
% 	fig_handle = show_summary_predictive_scores(all_results);
%
% * Displaying predictive scores with custom labels
% 	fig_handle = show_summary_predictive_scores(all_results, x_ticks, y_ticks, y_label);
% -------------------------------------------------------------------------
%% Author(s):
%   Antoine Valera
%
% -------------------------------------------------------------------------
%                               Notice
%
% This function was initially released as part of The SilverLab MatLab
% Imaging Software, an open-source application for controlling an
% Acousto-Optic Lens laser scanning microscope. The software was 
% developed in the laboratory of Prof Robin Angus Silver at University
% College London with funds from the NIH, ERC and Wellcome Trust.
%
% Copyright © 2015-2020 University College London
%
% Licensed under the Apache License, Version 2.0 (the "License");
% you may not use this file except in compliance with the License.
% You may obtain a copy of the License at
% 
%     http://www.apache.org/licenses/LICENSE-2.0
% 
% Unless required by applicable law or agreed to in writing, software
% distributed under the License is distributed on an "AS IS" BASIS,
% WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
% See the License for the specific language governing permissions and
% limitations under the License. 
% -------------------------------------------------------------------------
% Revision Date:
% 	09-06-2023
% -------------------------------------------------------------------------
% See also: 

% TODO : Extend the function to handle additional normalization methods.

function fig_handle = show_summary_predictive_scores(all_results, x_ticks, y_ticks, y_label)
    x_ticks = fix_labels(x_ticks);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    N_phates = size(results, 1);

    if nargin < 3 || isempty(y_ticks)
        y_ticks = cellfun(@(x) num2str(x), num2cell(1:N_phates), 'uni', false);
        tick_spacing = 1:2:N_phates*2;
    else        
        tick_spacing = 1:numel(y_ticks);
    end
    if nargin < 4 || isempty(y_label)
        y_label = 'Phate #';
    end
    
    

    fig_handle = figure(3);clf();
    
    subplot(1,2,1);
    
    imagesc(results ./ max(results));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility per behaviour)'};
    xticks(1:numel(x_ticks));xticklabels(cellfun(@(x) strrep(x,'_',' '), x_ticks, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');hold on;
    yticks(tick_spacing);hold on;yticklabels(y_ticks);
    ylabel(y_label); hold on;

    subplot(1,2,2);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    imagesc(results ./ max(results, [],2));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility in phate)'};
    xticks(1:numel(x_ticks));xticklabels(cellfun(@(x) strrep(x,'_',' '), x_ticks, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');
    yticks(tick_spacing);hold on;yticklabels(y_ticks);
    ylabel(y_label);
    
    set(fig_handle, 'Color', 'w'); hold on;
end