function groups = get_high_phate_ROIs(obj, N_phates, cutoff, bidirectional)
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

    groups = {};
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
            groups{counter}     = find(rois);
        end
    end
end