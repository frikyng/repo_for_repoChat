%% Plots lag chart of behaviours and extracted behaviours
% 	This function creates two figures. The first figure shows the prediction score 
%   as a function of lag for different behaviours. The second figure shows the 
%   autocorrelation as a function of lag for different extracted behaviours.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[fig_handle, fig_handle_cross] = plot_lag_chart(results, lag_list, sr, 
%                                                     behaviours, extracted_beh, 
%                                                     lag_smoothing)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	results(CELL ARRAY):
%                                   Results of prediction scores.
% 									Should be a 1-D cell array.
%
% 	lag_list(ARRAY) - Non-optional:
%                                   List of lag values used in the prediction.
% 									Should be a 1-D array.
%
%   sr(SCALAR) - Non-optional:
%                                   The sampling rate. Used to convert lag_list
%                                   into seconds.
%
%   behaviours(CELL ARRAY OF STRINGS) - Non-optional:
%                                   List of behaviours.
%
%   extracted_beh(STRUCTURE) - Non-optional:
%                                   A structure with field 'raw_behaviours'
%                                   containing the raw extracted behaviours and 
%                                   'original_behaviour_list' containing the 
%                                   list of original behaviours.
%
% 	lag_smoothing(SCALAR) -- Optional - Default value is 0:
%                                   If non-zero, apply a Gaussian smoothing
%                                   to the results. The value of lag_smoothing
%                                   determines the window width.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	fig_handle(FIGURE HANDLE):
%                                   Handle to the first figure plotting the 
%                                   prediction score as a function of lag.
%
% 	fig_handle_cross(FIGURE HANDLE):
%                                   Handle to the second figure plotting the 
%                                   autocorrelation of extracted behaviours.
%
% -------------------------------------------------------------------------
%% Examples:
% * Plot lag chart with lag smoothing
% 	[fig_handle, fig_handle_cross] = plot_lag_chart(results, lag_list, sr, 
%                                                      behaviours, extracted_beh, 
%                                                      lag_smoothing);
%
% * Plot lag chart without lag smoothing
% 	[fig_handle, fig_handle_cross] = plot_lag_chart(results, lag_list, sr, 
%                                                      behaviours, extracted_beh);
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
%	figure, plot, xline, tiledlayout, nexttile, xcorr
%
% TODO : Add support for 2D array input for results and behaviours.


function [fig_handle, fig_handle_cross] = plot_lag_chart(results, lag_list, sr, behaviours, extracted_beh, lag_smoothing)
    if nargin < 6 || isempty(lag_smoothing)
        lag_smoothing = 0;
    end
    results = cell2mat(results);
    if any(lag_smoothing)
        results = smoothdata(results,2,'gaussian',lag_smoothing);
    end
    x = fliplr(lag_list/sr);
    fig_handle = figure(1);clf();
    tiledlayout('flow');
    for beh = 1:numel(behaviours)
        axes(beh) = nexttile    ;  
        patch([x,fliplr(x)], [results(beh, :), zeros(size(results(beh, :)))], 'k');
        title(strrep(behaviours{beh},'_','\_'));xlabel('Lag (s)');ylabel('Prediction score');
        xline(0,'r--');
    end
    linkaxes(axes, 'xy')
    
    fig_handle_cross = figure(2);clf();
    tiledlayout('flow');
    for beh = 1:size(extracted_beh.raw_behaviours, 1)
        axes(beh) = nexttile    ;
        title(strrep(extracted_beh.original_behaviour_list{beh},'_','\_'));hold on
        plot(linspace(min(x), max(x), max(lag_list)*2+1), xcorr(extracted_beh.raw_behaviours(beh,:),max(lag_list),'normalized')); hold on
        xlabel('Lag (s)');xline(0,'r--')
    end
end

