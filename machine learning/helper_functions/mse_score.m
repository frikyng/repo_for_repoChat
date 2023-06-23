%% MSE (Mean Squared Error) Score Calculator
%   This function computes the Mean Squared Error (MSE) score between two sets
%   of data, typically a ground truth set and a prediction set. It also supports
%   optional weighting.
%
% -------------------------------------------------------------------------
%% Syntax:
%   [score] = mse_score(y_true, y_pred, w)
%
% -------------------------------------------------------------------------
%% Inputs:
%   y_true(Vector): 
%                   The true values or ground truth, as a vector. Column 
%                   vector is expected but row vector is also accepted.
%
%   y_pred(Vector): 
%                   The predicted values, as a vector. Should be of the 
%                   same size as y_true. Column vector is expected but row 
%                   vector is also accepted.
%
%   w(Vector) - Optional: 
%                   An optional weight vector of the same size as y_true and 
%                   y_pred. If not provided, equal weighting is assumed.
%
% -------------------------------------------------------------------------
%% Outputs:
%   score(Scalar) Double
%                   The calculated Mean Squared Error (MSE) score.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%   * This function is mainly used for regression model performance evaluation.
%   * The lower the MSE score, the closer the predicted values are to the 
%     actual values.
%
% -------------------------------------------------------------------------
%% Examples:
% * Compute MSE score for two vectors
%   score = mse_score([1,2,3], [1,2,3]);
%
% * Compute MSE score for two column vectors
%   score = mse_score([1;2;3], [1;2;3]);
%
% * Compute MSE score for two vectors with weights
%   score = mse_score([1,2,3], [1,2,3], [0.1, 0.3, 0.6]);
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
%   09-06-2023
% -------------------------------------------------------------------------
% See also: 


function [score] = mse_score(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    score = 1/mean((y_true - y_pred).^2); % Mean Squared Error
end