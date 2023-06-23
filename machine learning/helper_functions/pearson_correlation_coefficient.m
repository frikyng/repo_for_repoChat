%% Pearson Correlation Coefficient
%   Calculates the Pearson correlation coefficient between two vectors.
%
% -------------------------------------------------------------------------
%% Syntax:
%   [score] = pearson_correlation_coefficient(y_true, y_pred, w)
%
% -------------------------------------------------------------------------
%% Inputs:
%   y_true(Vector) Double:
%                                   A vector of true values. Should have the
%                                   same length as y_pred.
%
%   y_pred(Vector) Double:
%                                   A vector of predicted values. Should have
%                                   the same length as y_true.
%
%   w(NOT USED IN THIS FUNCTION) -- Optional:
%                                   This input argument is not used in this 
%                                   function. It is included for compatibility 
%                                   with other function signatures.
% -------------------------------------------------------------------------
%% Outputs:
%   score(Scalar) Double:
%                                   The Pearson correlation coefficient between 
%                                   y_true and y_pred. The value will range from 
%                                   -1 (perfect negative correlation) to 1 
%                                   (perfect positive correlation).
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * This function does not handle missing values. If y_true or y_pred
%   contains NaN values, you should handle them prior to calling this function.
% -------------------------------------------------------------------------
%% Examples:
% * Calculating correlation between two vectors
%   score = pearson_correlation_coefficient([1,2,3,4,5], [1,2,3,4,5]);
%
% -------------------------------------------------------------------------
%% Author(s):
%   Main/Initial Programmer
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
% Copyright © 2015-2023 University College London
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
%   09-06-2023
% -------------------------------------------------------------------------
% See also: 
%   corrcoef, corr

% TODO : Handle NaN values within the function

function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    score = corr(y_true, y_pred); % pearson r
end