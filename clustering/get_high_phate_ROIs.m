function ROI_groups = get_high_phate_ROIs(source, N_phates, cutoff, bidirectional, no_duplicates, valid_trace_idx)
    if isa(source, 'arboreal_scan_experiment')
        data            = source.dimensionality.LoadingsPM;
        valid_trace_idx = source.dimensionality.valid_trace_idx;
    else % ROI x N matrix, with N at least > N_Phates
        data            = source;
        valid_trace_idx = 1:size(source, 1); % note : if you do not pass a set list, ROis are relative indexes, not absolute.
    end
    if nargin < 2 || isempty(N_phates)
        N_phates = size(data, 2);
    end
    if nargin < 3 || isempty(cutoff) % set cutoff to 0 to get just the min or the max
        cutoff = 20;
    else
        if cutoff > 100
            error('cutoff must be a percentage between 0 and 100, or 0 to indicate just the min and the max, or a negative value to force the number of ROIs')
        end
    end
    if nargin < 4 || isempty(bidirectional)
        bidirectional = true;
    end
    if nargin < 5 || isempty(no_duplicates)
        no_duplicates = false;
    end
    if nargin >= 6 && ~isempty(valid_trace_idx)
        valid_trace_idx = valid_trace_idx;  % note : if you do not pass a set list, ROIs are relative indexes, not absolute. 
    end 

    ROI_groups = {};
    counter = 0;
    if cutoff >= 0
        cutoff = 100/cutoff;
    end
    if bidirectional
        direction = [1,2];
    else
        direction = 1;
    end
    for phat_comp = 1:N_phates
        for direc = direction
            counter     = counter + 1;
            loadings    = data(:,phat_comp);
            if direc == 1
                if cutoff >= 0
                    win         = range(data(:,phat_comp))/cutoff;
                    up_thr      = max(data(:,phat_comp)) - win;
                    rois        = data(:,phat_comp) >= up_thr;
                else
                    [~, rois]   = sort(data(:,phat_comp), 'descend');
                    rois        = rois(1:abs(cutoff));
                end
            elseif direc == 2
                if cutoff >= 0
                    win         = range(data(:,phat_comp))/cutoff;
                    low_thr     = min(data(:,phat_comp)) + win;
                    rois        = data(:,phat_comp) <= low_thr;
                else
                    [~, rois]   = sort(data(:,phat_comp));
                    rois        = rois(1:abs(cutoff));
                end
            end
            
            %% Store ROIs numbers
            valid_ROIs_idx          = find(valid_trace_idx); % ROIs that were actually used in the dim reduction
            ROI_groups{counter}     = valid_ROIs_idx(rois);
        end
    end
    
    if no_duplicates
        %% Find groups that are exact duplicates because that would make 2 times the same predictor
        for el = 1:(numel(ROI_groups)-1)
            for el2 = (el+1):numel(ROI_groups)
                if numel(ROI_groups{el}) == numel(ROI_groups{el2}) && all(ROI_groups{el} == ROI_groups{el2})
                    ROI_groups{el2} = [];
                end
            end
        end
    end    
end