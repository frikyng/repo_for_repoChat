function events = extract_events(all_traces, peak_loc, pre_peak_delay, post_peak_delay, exact_t)
    if nargin < 2 || isempty(peak_loc)
        [~, peak_loc, ~, ~] = detect_events(nanmedian(all_traces, 2));
    end
    if nargin < 3 || isempty(pre_peak_delay) 
        pre_peak_delay = 20;
    end
    if nargin < 4 || isempty(post_peak_delay) 
        post_peak_delay = 50;
    end
    if nargin < 5 || isempty(exact_t)
        exact_t = [];
    end
        
    events = [];    
    
    %% In case you are in a case with no event > thr
    if isempty(peak_loc) || all(isnan(peak_loc))
        return
    end
    
    %% For 99% of the other cases
    for peak = 1:size(peak_loc, 2)
        
        pr_pd = pre_peak_delay;
        po_pd = post_peak_delay;
        
        if (peak_loc(peak)-pr_pd) <= 1 
            pr_pd = peak_loc(peak) - 1;
        elseif (peak_loc(peak)+po_pd) > size(all_traces, 1)
            po_pd = size(all_traces, 1)- peak_loc(peak);
        end
        
        if isempty(exact_t)
            % standard case, centred around peak location of the median
            % or mean trace
            new = all_traces((peak_loc(peak)-pr_pd):(peak_loc(peak)+po_pd),:);
            if pr_pd ~= pre_peak_delay
               new = [NaN(pre_peak_delay - pr_pd, size(new, 2)); new];
            elseif po_pd ~= post_peak_delay
               new = [new; NaN(post_peak_delay - po_pd, size(new, 2))];
            end
        else
            % custom case, obtained after fitting individual events,
            % and adjusted for jitter
            new = repmat((peak_loc(peak)-pr_pd):(peak_loc(peak)+po_pd), size(all_traces, 2), 1) + exact_t{peak};                
        end

        events = cat(3, events, new);
    end
end