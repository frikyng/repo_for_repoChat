%% Scale all the recordings and return the best scale and offset for each.
%  This function scales every recordings based on the best scaling factor 
%  per bin matching the overall median. The function also optimizes 
%  the offset and scale of each recording to minimize the error 
%  with respect to a reference trace.
%
% -------------------------------------------------------------------------
%% Syntax:
%  [global_scaling, global_offset, best_ind_scal, best_ind_offset, N_pks, 
%  N_pt_bsl] = scale_every_recordings(all_traces_per_rec, demo, pk_locs, 
%  bsl_range, smoothing, weighted)
%
% -------------------------------------------------------------------------
%% Inputs:
% 	all_traces_per_rec(cell array of matrices):
%       Each cell contains a matrix of traces for a single recording. 
%       Each column in the matrix represents an individual trace. 
%
% 	demo(integer) - Optional - Default is 0:
%       Debugging parameter. Controls the verbosity level.
% 
% 	pk_locs(vector) - Optional - Default is autoestimated:
%       A vector containing the locations of the peaks in the data. 
%       If not provided, peaks will be auto-estimated.
% 
% 	bsl_range(vector) - Optional - Default is autoestimated:
%       A vector containing the range of baseline data points. 
%       If not provided, the range will be autoestimated.
% 
% 	smoothing(integer) - Optional - Default is 0:
%       The window size for Gaussian smoothing. If zero or empty, 
%       no smoothing is performed.
% 
% 	weighted(boolean) - Optional - Default is true:
%       If true, a weighted mean is used for calculating global scaling 
%       and offset. Otherwise, a simple median is used.
%
% -------------------------------------------------------------------------
%% Outputs:
% 	global_scaling(vector) - double:
%       A vector containing the best overall scaling factor for each trace.
%
% 	global_offset(vector) - double:
%       A vector containing the best overall offset for each trace.
%
% 	best_ind_scal(cell array of vectors) - double:
%       Each cell contains a vector of the best individual scaling factors
%       for each trace in a given recording.
%
% 	best_ind_offset(cell array of vectors) - double:
%       Each cell contains a vector of the best individual offsets for 
%       each trace in a given recording.
%
% 	N_pks(vector) - integer:
%       Number of peaks in each recording.
%
% 	N_pt_bsl(vector) - integer:
%       Number of baseline points in each recording.
%
% -------------------------------------------------------------------------
%% Extra Notes:
%
% * This function uses parallel processing to speed up the computations. 
%   The number of cores used is determined by the demo parameter and the 
%   number of cores available on the machine.
%
% * The function uses optimization to determine the best scale and offset 
%   parameters for each trace in a recording.
%
% -------------------------------------------------------------------------
%% Examples:
% * Apply scaling to a set of recordings and display intermediate steps
% 	[global_scaling, global_offset] = scale_every_recordings(all_traces_per_rec, 2);
%
% * Apply scaling without smoothing and without weighted mean
% 	[global_scaling, global_offset] = scale_every_recordings(all_traces_per_rec, 0, [], [], 0, false);
%
% * Specify peak locations and baseline range
% 	[global_scaling, global_offset] = scale_every_recordings(all_traces_per_rec, 0, pk_locs, bsl_range);
% -------------------------------------------------------------------------
%%                               Notice
%
%% Author(s):
%   Main/Initial Programmer, Other Programmer 2, etc...
%	if using code for elsewhere, mention the original authors here.
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
% 	DD-MM-YYYY
% -------------------------------------------------------------------------
% See also: 
%	parent_related_function_1, parent_related_function_2,
% 	children_related_function_1, children_related_function_2 ...

%% all_traces_per_rec
% 1 x N trials cell array of t x M ROI matrices to rescale
%% demo
% if 1 ...
% if 2 ...
% pk_locs
%% If provided, rescaling is done using  peak location
% bsl_range
%% If provided, offt is computed on these datapoints only
% smoothing
%% If ~0 , presmooth traces before scaling 
% weighted
% If true, offset and peaks estimate will be weighted based on the number
% of datapoints / events (respectively)


% This function scales all the provided recordings, such that their individual
% peaks and baselines match a common reference derived from all recordings.

function [global_scaling, global_offset, best_ind_scal, best_ind_offset, N_pks, N_pt_bsl] = scale_every_recordings(all_traces_per_rec, demo, pk_locs, bsl_range, smoothing, weighted)
    
    % Check if demo parameter is provided, else set default
    if nargin < 2 || isempty(demo)
        demo = 0;
    end

    % Check if pk_locs is provided, else auto-estimate peak detection level
    if nargin < 3 || isempty(pk_locs)

        % Filter signal and remove noise
        representative_trace    = fillmissing(representative_trace,'linear');
        ref_signal              = wdenoise(double(representative_trace));
        noise                   = representative_trace - ref_signal;

        % Calculate envelope of noise for threshold
        [top_noise, ~]          = envelope(noise,50,'rms'); % used to be 20, 'peak'
        initial_thr             = mean(top_noise)*5;    

        % Detect peaks in the denoised signal
        [pks, pk_locs]             = findpeaks(wavelet_denoise(representative_trace')', 'MinPeakProminence', initial_thr);
    end

    % Check if bsl_range is provided, else throw error for now
    if nargin < 4 || isempty(bsl_range)
        error('to do')    
        bsl_percentile = '?';
    else
        bsl_percentile = 50; % Use the median value since we already use a range for the baseline
    end

    % Apply smoothing if smoothing argument is given
    if nargin < 5 || isempty(smoothing) || ~any(smoothing)
        % no smoothing, pass
    else
        all_traces_per_rec = cellfun(@(x) smoothdata(x,'gaussian',smoothing), all_traces_per_rec, 'UniformOutput', false);
    end

    % Check if weighted parameter is provided, else set default
    if nargin < 6 || isempty(weighted)
        weighted = true;
    end

    % Function to obtain reference trace
    consensus_trace_func    = @nanmedian;

    % Define options for optimization
    options                 = optimset;

    % Prepare for parallel processing. Use number of available cores, unless in demo mode
    ncores                  = (feature('numcores') * double(~(demo >= 2)));
    tp                      = [cellfun(@(x) size(x, 1), all_traces_per_rec), inf];
    N_pks                   = zeros(size(all_traces_per_rec))'; % for weights, if you use them
    N_pt_bsl                = zeros(size(all_traces_per_rec))'; % for weights, if you use them
    
    parfor (rec = 1:numel(all_traces_per_rec), ncores)
    %for rec = 1:numel(all_traces_per_rec)

        all_traces_in_rec       = all_traces_per_rec{rec};
        representative_trace    = consensus_trace_func(all_traces_in_rec, 2);

        tp_range        = [(sum(tp(1:rec-1))+1),sum(tp(1:rec))];
        pk_tp_current   = pk_locs(pk_locs > tp_range(1) & pk_locs < tp_range(2)) - tp_range(1)+1;
        bsl_tp_current  = bsl_range(bsl_range > tp_range(1) & bsl_range < tp_range(2)) - tp_range(1)+1;      
%         figure();plot(representative_trace); hold on; scatter(pk_tp_current, representative_trace(pk_tp_current),'ko', 'filled')
%         figure();plot(representative_trace); hold on; scatter(bsl_tp_current, representative_trace(bsl_tp_current),'kv', 'filled')
        N_pks(rec)      = numel(pk_tp_current);
        N_pt_bsl(rec)   = numel(bsl_tp_current);
        %N_pt_bsl        = 1
        
        
        %% Set representative_trace baseline close to 0
        % Adjust offset so both traces baselines are at 0     
        if isempty(bsl_tp_current)
            bsl_tp_current = prctile(representative_trace, 10);
            bsl_tp_current = find(representative_trace < bsl_tp_current);
        end   
        % Loop over each trace in the recording
        bsl = zeros(size(all_traces_in_rec(:,1)));

        best_offset_for_ref     = fminbnd(@(f) find_traces_offset_func(f, bsl(bsl_tp_current), representative_trace(bsl_tp_current), demo == 2, bsl_percentile), 0.01, 99.9999, options); 
        best_offset_for_ref     = prctile(representative_trace(bsl_tp_current), best_offset_for_ref);
        representative_trace    = representative_trace - best_offset_for_ref;
     
        best_ind_scal{rec}    = [];
        best_ind_offset{rec}  = [];

        for trace_idx =  1:size(all_traces_in_rec, 2) 

            % Skip if trace is all NaNs
            if all(isnan(all_traces_in_rec(:,trace_idx)))
                best_ind_offset{rec}(trace_idx) = NaN;
                best_ind_scal{rec}(trace_idx)   = NaN;
            else
                % Find best percentile to subtstract baseline traces, then
                % convert to baseline offset
                best_ind_offset{rec}(trace_idx)     = fminbnd(@(f) find_traces_offset_func(f, bsl(bsl_tp_current), all_traces_in_rec(bsl_tp_current,trace_idx), demo == 2 && trace_idx == subset_for_demo, bsl_percentile), 0.01, 99.9999, options); 
                best_ind_offset{rec}(trace_idx)     = prctile(all_traces_in_rec(bsl_tp_current,trace_idx), best_ind_offset{rec}(trace_idx));
                all_traces_in_rec(bsl_tp_current,trace_idx) = all_traces_in_rec(bsl_tp_current,trace_idx) - best_ind_offset{rec}(trace_idx);

                if ~isempty(pk_tp_current)
                    try
                        best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, representative_trace, all_traces_in_rec(:,trace_idx), demo >= 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                    catch 
                        best_ind_scal{rec}(trace_idx) = fminbnd(@(f) scale_trace_func(f, representative_trace, all_traces_in_rec(:,trace_idx), demo >= 2 && trace_idx == subset_for_demo, pk_tp_current), 1e-3, 100, options); % scaling factor must be > 0
                    end
                    if demo >= 3
                        uiwait(figure(6663))
                    end
                else                    
                    best_ind_scal{rec}(trace_idx) = NaN;
                end
                
            end
        end
    end
    clear data_f_idx % just to remove parfor warning

    % Compute the best scaling factor per bin matching the overall median
    if ~weighted
        global_scaling   = nanmedian(vertcat(best_ind_scal{:}), 1);
        global_offset    = nanmedian(vertcat(best_ind_offset{:}), 1);
        N_pks(:)         = NaN;
        N_pt_bsl(:)      = NaN;
    else
        global_scaling = w_mean(vertcat(best_ind_scal{:}), N_pks);
        if N_pt_bsl == 0 
            N_pt_bsl = 1
        end
        global_offset   = w_mean(vertcat(best_ind_offset{:}), N_pt_bsl); 
    end    
end

% Function to calculate weighted mean
function out = w_mean(var, weights)
    dim = 1;
    out = nansum(weights.*var,dim)./nansum(weights,dim);
    out(out == 0) = NaN;
end
