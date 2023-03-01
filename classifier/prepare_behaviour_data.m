%% In ML code, preparing behaviour takes time. This function can be used to list the behaviours once (No TP filter applied) and pass it to 

function [stuct_out, raw_behaviours, beh_thr, formatted_behaviour_list] = prepare_behaviour_data(obj, behaviour_list, detrend_behaviours, smooth_behaviours)
    if nargin < 2 || isempty(behaviour_list)
        behaviour_list   = obj.behaviours.types;
    elseif ~iscell(behaviour_list)
        behaviour_list = {behaviour_list};
    end
    if nargin < 3 || isempty(detrend_behaviours)    
        detrend_behaviours      = true;
    end
    if nargin < 4 || isempty(smooth_behaviours)    
        smooth_behaviours       = 1;
    end
    
    pt_per_s                = 1/nanmedian(obj.timescale.sr);    
    raw_behaviours          = [];
    beh_thr                 = [];
    for type_idx = 1:numel(behaviour_list)
        type = behaviour_list(type_idx);
        
        %% Get original behaviour
        warning('off')
        if contains(type{1}, 'baseline')
            % Never detrend baseline
            [~,~,original_beh]   = obj.get_behaviours(type{1},'',false,true,true);            
            %original_beh.value = original_beh.value-nanmean(original_beh.value);
        else
            [~,~,original_beh]   = obj.get_behaviours(type{1},'',detrend_behaviours,true,true);
        end
        
        %% Adjust name
        if iscell(type) && iscell(type{1})
            type = type{1}{1};
            behaviour_list{type_idx} = type;
        end
        if isempty(original_beh.value)
            warning(['Behaviour ',type{1},' Not found. Check for Typo and see if it is listed in Obj.behaviours']);
            continue
        end
        
        %% Median smoothing on all behaviour but trigger to remove small blips
        if contains(type, 'trigger') && smooth_behaviours
            current_beh         = smoothdata(original_beh.value, 'gaussian', [smooth_behaviours*pt_per_s/nanmedian(diff(obj.t)),0]); 
        elseif contains(type, 'baseline')
            current_beh         = original_beh.value;
        elseif ~smooth_behaviours
            current_beh         = original_beh.value;
            current_beh(isnan(current_beh)) = 0;
            current_beh(current_beh == 0) = randn(sum(current_beh == 0),1) * (rms(current_beh)/100);
        else
            current_beh         = smoothdata(original_beh.value, 'movmedian', smooth_behaviours*pt_per_s);
            
            %% Gaussian smoothing to denoise behaviour
            current_beh         = smoothdata(current_beh, 'gaussian', smooth_behaviours*pt_per_s*5); 
            
            current_beh(current_beh == 0) = randn(sum(current_beh == 0),1) * (rms(current_beh)/100);
        end

        
        %% Threshold to binarize behaviour for the classifier
        if any(contains(type, {'RT3D_MC','BodyCam_Eye','BodyCam_Laser'}))
            beh_thr(type_idx) = nanmedian(current_beh);
        else
            beh_thr(type_idx) = nanmax(current_beh) - (range(current_beh) * 0.9);
        end

        %% Plot behaviour and threshold
        %         above = current_beh > beh_thr(type_idx);
        %         below = ~above;
        %         above_beh = current_beh;above_beh(below) = NaN;
        %         below_beh = current_beh;below_beh(above) = NaN;
        %figure();plot(above_beh, 'g');hold on;plot(below_beh, 'r');hold on;title(type{1});hold on; plot([0,numel(current_beh)],[thr, thr],'k--');

        raw_behaviours      = [raw_behaviours; current_beh];
    end
    formatted_behaviour_list   = strrep(behaviour_list, '_', '\_'); % reformat strings to be usable in titles and legends
    
    %% Build output structure if you want to re-use it across iterations
    stuct_out.original_behaviour_list  = behaviour_list;
    stuct_out.formatted_behaviour_list = formatted_behaviour_list;
    stuct_out.raw_behaviours           = raw_behaviours;
    stuct_out.beh_thr                  = beh_thr;
end

