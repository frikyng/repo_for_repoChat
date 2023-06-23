%% Computes the inverse of the Root Mean Squared Error (RMSE) between actual and predicted values.
% 	Assumes the y_true and y_pred are vectors of the same length. Also, it
%   is indifferent to row or column orientation of vectors. The function 
%   treats row vectors as column vectors.
%
% -------------------------------------------------------------------------
%% Syntax:
% 	[score] = rmse_score(y_true, y_pred, w)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	y_true(Vector) Numeric:
%                                   The vector of actual values. It should 
%                                   be numeric. It could be row or column vector.
%
% 	y_pred(Vector) Numeric:
%                                   The vector of predicted values. It 
%                                   should be numeric. It could be row 
%                                   or column vector.
%
% 	w - Unused
% -------------------------------------------------------------------------
%% Outputs:
% 	score(Scalar) Numeric
%                                   The computed inverse of the RMSE. It 
%                                   is a measure of the differences between
%                                   values (sample or population values) 
%                                   predicted by a model and the values 
%                                   observed in reality.
%
% -------------------------------------------------------------------------
%% Examples:
% * Simple usage
% 	y_true = [1, 2, 3, 4, 5];
% 	y_pred = [1.1, 1.9, 3.1, 3.9, 5.1];
% 	score = rmse_score(y_true, y_pred, []);
%
% * Usage with column vectors
% 	y_true = [1; 2; 3; 4; 5];
% 	y_pred = [1.1; 1.9; 3.1; 3.9; 5.1];
% 	score = rmse_score(y_true, y_pred, []);
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
%	rms, sqrt, mean, isrow

% TODO : Currently the 'w' parameter is not being used in this function. It 
% may be worth considering if there is a future need for this parameter, 
% or if it should be removed.

function [score] = rmse_score(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    score = 1/sqrt(mean((y_true - y_pred).^2)); % Root Mean Squared Error
end