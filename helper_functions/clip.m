%% One Line summary
% 	Extended exlanation, context, documentation, and use cases
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[out1, out2, out3 ...] = function(in1, in2, in3, varargin)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	in1(DATA TYPE):
%       Explanation of the role of the input and range of valid values.
%   	Explanation follows this level of identation.
%       * Bullet points can be used (with two extra space indentation if
%         more than one line) 
%
% 	in2(OPTIONAL INPUT DATA TYPE) - Optional - Default is xxx:
%       Explanation of the role of the input and of each available option if
%       several options or data types are valid. In this case, use bullet
%       points for each scenario
%
% 	in3(OPTIONAL INPUT DATA TYPE ALLOWING WITH MULTIPLE INPUT CASE) --
%       Optional - any in {'bla', 'bli', 'blu'}:
%       Explanation of the role of the input and of each available option 
%       if several options are possible
% 		* if 'bla', then ...
% 		* if 'bli', then ...
% 		* if 'blu', then ...
%
% 	varargin(accepted object class AND/OR {'Argument',Value} pairs):
%       Any pair of 'argument' and value from the following list. for 
%       details and default values, see ClassName. This is only to be used 
%       when a ****_params.file is build from varargin. In this case, the 
%       relevant options are cited in the subsection below.
% --------------
% 			{'option name 1' (DATA TYPE)}: Default is 'xxx'. any in
%           	{'xxx','yyy','zzz'}
%               	Explanation of what the option does.
%                   *'xxx' does this
%                   *'yyy' does that
%                   *'zzz' etc...
%
% 			{' option name 2' (DATA TYPE)}: Default is false.
%                   Explanation of what the option does.
%
%
% -------------------------------------------------------------------------
%% Outputs:
% 	out1([FORMAT]) DATA TYPE
%       Explanation of the output
%
% 	out2([FORMAT]) DATA TYPE
%       Explanation of the output
%
% 	out3([FORMAT]) DATA TYPE
%       Explanation of the output
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * Longer explanation of the function behaviour or things you should know
% 	before using it
% -------------------------------------------------------------------------
%% Examples:
% * Example 1 title
% 	out1 = function(in1);
%
% * Example 2 title
% 	[out1, out2] = function(in1, in2, in3);
%
% * Example 3 title
% 	[out1, out2] = function('', '', '', 'Option1', value_option_1);
% -------------------------------------------------------------------------
%% Author(s):
%   Main/Initial Programmer, Other Programmer 2, etc...
%	if using code for elsewhere, mention the original authors here.
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
% 	DD-MM-YYYY
% -------------------------------------------------------------------------
% See also: 
%	parent_related_function_1, parent_related_function_2,
% 	children_related_function_1, children_related_function_2 ...

% TODO : Know issue or limit cases, future developments, TODO list
% 	etc...

function inputArray = clip(inputArray, thresholdValues)
    % Check if thresholdValues is provided, otherwise use default values
    if nargin < 2 || isempty(thresholdValues)
        thresholdValues = [0, 1];
    end
    
    % Clip values in inputArray based on thresholdValues
    if numel(thresholdValues) == 1
        inputArray(inputArray > thresholdValues(1)) = thresholdValues(1);
    else
        inputArray(inputArray < thresholdValues(1)) = thresholdValues(1);
        inputArray(inputArray > thresholdValues(2)) = thresholdValues(2);
    end
end