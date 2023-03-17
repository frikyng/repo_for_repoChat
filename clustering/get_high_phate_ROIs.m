function ROI_groups = get_high_phate_ROIs(obj, N_phates, cutoff, bidirectional, no_duplicates)
    if nargin < 2 || isempty(N_phates)
        N_phates = size(obj.dimensionality.LoadingsPM, 2);
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
            loadings    = obj.dimensionality.LoadingsPM(:,phat_comp);
            if direc == 1
                if cutoff >= 0
                    win         = range(obj.dimensionality.LoadingsPM(:,phat_comp))/cutoff;
                    up_thr      = max(obj.dimensionality.LoadingsPM(:,phat_comp)) - win;
                    rois        = obj.dimensionality.LoadingsPM(:,phat_comp) >= up_thr;
                else
                    [~, rois]   = sort(obj.dimensionality.LoadingsPM(:,phat_comp), 'descend');
                    rois        = rois(1:abs(cutoff));
                end
            elseif direc == 2
                if cutoff >= 0
                    win         = range(obj.dimensionality.LoadingsPM(:,phat_comp))/cutoff;
                    low_thr     = min(obj.dimensionality.LoadingsPM(:,phat_comp)) + win;
                    rois        = obj.dimensionality.LoadingsPM(:,phat_comp) <= low_thr;
                else
                    [~, rois]   = sort(obj.dimensionality.LoadingsPM(:,phat_comp));
                    rois        = rois(1:abs(cutoff));
                end
            end
            
            %% Store ROIs numbers
            valid_ROIs_idx      = find(obj.dimensionality.valid_trace_idx); % ROIs that were actually used in the dim reduction
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