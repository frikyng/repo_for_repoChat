%% One Line summary
% 	Performs wavelet denoising on data, handling missing (NaN) values.

% -------------------------------------------------------------------------
%% Syntax:
% 	[data, missing] = wavelet_denoise(data, varargin)

% -------------------------------------------------------------------------
%% Inputs:
% 	data(MATRIX):
%       Data matrix to be denoised. NaN values indicate missing data points.
%       This function operates on each column independently.
%   	
% 	varargin(Passed to wdenoise function):
%       Any parameters to be passed on to the MATLAB wdenoise function. 
%       See the wdenoise documentation for details on available parameters.
%       This could be the type of wavelet, thresholding rule, etc.

% -------------------------------------------------------------------------
%% Outputs:
% 	data([Same as input]) MATRIX
%       The denoised data. Columns are independently denoised. 
%       The data is converted to double for processing and back to single 
%       for output. Missing values (NaN) in the input are linearly 
%       interpolated for processing and then set back to NaN in the output.

% 	missing([Same size as input]) LOGICAL ARRAY
%       A logical array indicating the position of missing values in the 
%       input data.

% -------------------------------------------------------------------------
%% Extra Notes:

% * The function operates on each column of the input data independently.
% * Missing data (NaN values) are linearly interpolated before processing
%   and then set back to NaN in the output. 
% * The 'wdenoise' function from MATLAB's Wavelet Toolbox is used for the 
%   denoising process.

% -------------------------------------------------------------------------
%% Examples:
% * Denoising a 2D data matrix using default parameters
% 	[data, missing] = wavelet_denoise(data);
%
% * Denoising with a custom wavelet type and threshold rule
% 	[data, missing] = wavelet_denoise(data, 'Wavelet', 'sym4', 
%                                      'ThresholdRule', 'Soft');
% -------------------------------------------------------------------------
%%                               Notice
%
%% Author(s):
%   Antoine Valera
%                 -----------------------------------------
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
% 	13-06-2023
% -------------------------------------------------------------------------
% See also: 
%	fillmissing, isnan, wdenoise


function [data, missing] = wavelet_denoise(data, varargin)
    missing             = isnan(data);
    data                = fillmissing(data,'linear');
    keep                = ~all(isnan(data)); % When one trace is full of NaN, we can't analyse it
    data(:,keep)        = single(wdenoise(double(data(:,keep)), varargin{:}));
    data(missing)       = NaN;
end

%     figure();hold on;
%     plot(data,'k');hold on;
%     plot(wdenoise(double(data)),'r');hold on
%     plot(wdenoise(double(data),'DenoisingMethod','BlockJS'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','FDR'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','Minimax'));hold on
%     plot(wdenoise(double(data),'DenoisingMethod','Sure'));hold on
%     legend({'original','Bayes','BlockJS','FDR','Minimax','Sure'}) ; % BAYES, FDR and SURE look good