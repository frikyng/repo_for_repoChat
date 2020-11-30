%% Check for problems in your meta analysis

function quality_control(varargin)
    if varargin{1} == 1
        %% Check for consistent number of bins
        all_traces_per_bin = varargin{2};
%         if numel(unique(cellfun(@(x) size(x, 2),all_traces_per_bin))) > 1
%             error_box('Some of your recordings have a different number of ROIs... you should check that')
%         end

    elseif varargin{1} == 2
        %% Recordings with different sampling rates can cause some issues
        all_sr = varargin{2};
        if numel(unique(round(1./all_sr, 1))) > 1
            error_box('Very variable sampling rate... you should check that')
        end
    end

end

