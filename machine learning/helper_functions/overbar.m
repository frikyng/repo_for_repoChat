%% Draws an overbar over data and displays text
% This function plots a horizontal line (overbar) at a specific height
% between two x-values and then places a text label on the overbar.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[hl, ht] = overbar(x1, x2, y, txt)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	x1(Double):
%                                   Starting x-coordinate of the overbar
%
% 	x2(Double):
%                                   Ending x-coordinate of the overbar
%
% 	y(Double):
%                                   Height of the overbar
%
% 	txt(String):
%                                   Text to be displayed on the overbar
%
% -------------------------------------------------------------------------
%% Outputs:
% 	hl(Line Object):
%                                   Handle to the overbar line object
%
% 	ht(Text Object):
%                                   Handle to the text object
%
% -------------------------------------------------------------------------
%% Extra Notes:
% * The size of the hook (d) is set to 1. Adjust this according to your
%   y-axis scaling.
% -------------------------------------------------------------------------
%% Examples:
% * Draw an overbar from x=1 to x=5 at y=3 with text "Hello"
% 	[hl, ht] = overbar(1, 5, 3, 'Hello');
%
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
%	line, text

function [hl,ht] = overbar(x1, x2, y, txt)
    sz = get(gca,'FontSize');
    %bg = get(gca,'Color');
    d = 1; % size of hook, change depending on y axis scaling
    hl = line([x1,x1,x2,x2], [y,y+d,y+d,y], 'LineWidth', 2, 'Color',[0.6,0.6,0.6]);
    ht = text((x1+x2)/2, y+d, txt, ...
              'HorizontalAlignment','center', ...
              'VerticalAlignment','middle', ... 
              'FontSize',sz*2, ...
              'BackgroundColor','none');
end
