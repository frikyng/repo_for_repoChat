%% copied from explore_factors.mlapp .. to merge
%% see also load_bAP

function load_specific_trial(obj, pts_of_origin, ROIs_to_load)
    if ~iscell(ROIs_to_load)
        ROIs_to_load = {ROIs_to_load};
    end
    
    %% Get point of interest
    pt = ginput(1);
    
    %% Identify recod where this point belong
    tp_timescale = find(pts_of_origin);
    t = tp_timescale(round(pt(1)));
    tp_starts = cumsum([1, obj.timescale.tp]);
    record = find(tp_starts < t, 1, 'last');
    
    %% Identify the timepoint within this record
    tp_in_record = t - tp_starts(record);
    
    %% Get path for reload
    source = obj.arboreal_scans{record}.data_folder;    
    
    for set = ROIs_to_load
        %% Set analysis_params, for loading
        ap = obj.arboreal_scans{record}.analysis_params;
        p = analysis_params(   'data_folder'    ,ap.data_folder,...
                               'data_type'      , 'concatenated',...
                               'rendering_mode' , 'mosaic',...
                               'tracing_source' , 'swc'         , ...
                               'signal_channel' , 2             ,...
                               'registration_channel', 2,.... % to fix
                               'registration'   , true          ,...
                               'registration_reload','auto_offsets.mat',...
                               'smoothing'      , [10,0]        ,...
                               'mask_method'    , ''   ,...
                               'mask_value'     , ''     ,...
                               'compression'    , '2D'          ,...
                               'rendering'      , true         ,...
                               'ROIs'           ,set{1},...
                               'z_score'        ,-5);
        data = load_ribbon_scan(p); 
    end
end

