%% see also load_specific_trial

function im_out = load_bap(expe, record, pt_in_record, extra_window)
    if nargin < 4 || isempty(extra_window)
        extra_window = [0,0];
    elseif numel(extra_window) == 1
        extra_window = [extra_window, extra_window];
    elseif numel(extra_window) ~= 2
        error('the window input must be an 1x1 INT (for symetric windows) or a 2x1 INT for asymetric windows')
    end

    ap = expe.arboreal_scans{record}.analysis_params;
    p = analysis_params('data_folder',ap.data_folder,...
                        'data_type'      , 'concatenated',...
                           'rendering_mode',  'mosaic',...
                           'tracing_source' , 'swc'         , ...
                           'signal_channel' , 2             ,...
                           'registration_channel', 2,.... % to fix
                           'registration'   , true          ,...
                           'registration_reload','auto_offsets.mat',...
                           'smoothing'      , [10,0]        ,...
                           'mask_method'    , ''   ,...
                           'mask_value'     , ''     ,...
                           'compression'    , '2D'          ,...
                           'rendering'      , false         ,...
                           'z_score'        ,-5,...
                           'ROIs', 0); 

%     if app.ShowexcludedROIsCheckBox.Value
%         p.ROIs = 0;
%     end

%     if record == app.current_record_number
%         data = app.current_record_data;
%     else
        data = load_ribbon_scan(p);   
%         app.current_record_number = record;
%         app.current_record_data = data;
%     end


    
    im_out = {};
    for ev = 1:numel(pt_in_record)
        pt = pt_in_record{ev};
        pt = (pt(1)-extra_window(1)):(pt(end) + extra_window(2));
        pt = pt(pt > 0 & pt <= size(data{1}, 4));            
        im = max_projection(data{1}(:,:,:,pt,:),'','mean');                
        im_out{ev} = im(:,:,1);
        drawnow()
    end

    d           = reshape(data{1}(:,:,:,:,1),[],size(data{1},4));
    d(isinf(d)) = NaN;
    d           = nanmean(d);
    figure();plot(d); hold on;
    for ev = 1:numel(pt_in_record)
        pt = pt_in_record{ev};
        pt = (pt(1)-extra_window(1)):(pt(end) + extra_window(2));
        pt = pt(pt > 0 & pt <= size(data{1}, 4));            
        scatter(pt, d(pt), 'ro')
    end
end