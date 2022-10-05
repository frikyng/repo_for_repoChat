%% Classify behaviour for a given experiment, for all behaviours
%% use expe 07-26 to fix gaps in behaviour

%% Cell 2019-09-17_exp_1 anticorelated to running

function out = predict_behaviours(obj, use_classifier, method, type_of_trace, behaviour_list, ROI_groups, varargin)
    if nargin < 1 || isempty(obj)
        obj = '' 
    end
    if nargin < 2 || isempty(use_classifier)
        use_classifier = true; % if false, use regression learner
    end
    if nargin < 3 || isempty(method)
        method = 'svm'; % i.e. fitcsvm vs fitclinear etc
    end
    if nargin < 4 || isempty(type_of_trace)
        type_of_trace   = 'subtracted_peaks'; % ['subtracted' OR 'rescaled' OR 'raw'] AND ['peaks' or '']. eg 'subtracted_peaks' , or 'raw'
    end
    if nargin < 5 || isempty(behaviour_list)
        behaviour_list   = obj.behaviours.types;
    elseif ~iscell(behaviour_list)
        behaviour_list = {behaviour_list};
    end
    if nargin < 6 || isempty(ROI_groups)
        ROI_groups   = num2cell(obj.ref.indices.valid_swc_rois);
    elseif ~iscell(ROI_groups)
        ROI_groups   = num2cell(ROI_groups);
    end    
    if nargin < 7 || isempty(varargin)
        parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        parameters = DEFAULT_CLASSIFIER_OPTION(varargin);
    end
    
    use_hd_data     = false;
    time_filter     = 0;

    %% Make sure preprocessing was done correctly
    rendering       = obj.rendering;
    obj.rendering   = false;
    [obj, source_signal, ~, timepoints] = prepare_phate_analysis(obj, use_hd_data, time_filter, type_of_trace);
    obj.rendering   = rendering;

    %% Get the signal for the selected timepoints and ROIs
    %Valid_ROIs      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); % remove excluded branches AND bad_ROIs based on correlation)
    invalid_ROIs     = ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list); % remove excluded branches AND bad_ROIs based on correlation)
    source_signal(:, invalid_ROIs) = NaN;
    data             = NaN(numel(timepoints), numel(ROI_groups));
    for gp_idx = 1:numel(ROI_groups)
        data(:, gp_idx) =  nanmean(source_signal(timepoints, ROI_groups{gp_idx}),2)';        
    end
    
    data(:, all(isnan(data),1)) = [];
    data = data';

    %     %% Correlation of the different ROIs with each other
    %     figure();imagesc(corr(Ca'))
    % 
    %     %% Correlation of the Ca2+ values per event
    %     figure();imagesc(corr(Ca))

    %% Get all behaviours at these tp
    behaviours              = [];
    raw_behaviours          = [];
    processed_behaviours    = [];
    
    for type_idx = 1:numel(behaviour_list)
        type = behaviour_list(type_idx);
        
        %% Get original bhaviour
        [~,~,original_beh]   = obj.get_behaviours(type{1},'',true,true,true);
        if iscell(type) && iscell(type{1})
            type = type{1}{1};
            behaviour_list{type_idx} = type;
        end
        
        %% Median smoothing on all behaviour but trigger to remove small blips
        if ~contains(type, 'trigger')
            current_beh         = smoothdata(original_beh.value, 'movmedian', 10); 
        else
            current_beh         = original_beh.value;
        end

        %% Gaussian smoothin to denoise behaviour
        current_beh         = smoothdata(current_beh, 'gaussian', 50);

        %% Behaviour normalization
        %         if ~contains(type, 'trigger')
        %             try
        %                 beh = normalize(beh,'medianiqr') ; % this fails when there are a lot of similar values
        %             catch
        %                 beh = beh;
        %             end
        %         else
        %             beh = beh;
        %         end

        %% Threshold to binarize behaviour for the classifier
        if any(contains(type, {'RT3D_MC','BodyCam_Eye','BodyCam_Laser'}))
            thr = nanmedian(current_beh);
        else
            thr = nanmax(current_beh) - (range(current_beh) * 0.9);
        end

        %% Plot behaviour and threshold
        above = current_beh > thr;
        below = ~above;
        above_beh = current_beh;above_beh(below) = NaN;
        below_beh = current_beh;below_beh(above) = NaN;
        %figure();plot(above_beh, 'g');hold on;plot(below_beh, 'r');hold on;title(type{1});hold on; plot([0,numel(current_beh)],[thr, thr],'k--');

        raw_behaviours      = [raw_behaviours; current_beh];
        behaviours          = [behaviours; current_beh(timepoints)];
        if use_classifier
            processed_behaviours= logical([processed_behaviours; current_beh(timepoints) > thr]);
        else
            processed_behaviours= [processed_behaviours; current_beh(timepoints)];
        end
    end
    behaviour_list   = strrep(behaviour_list, '_', '\_'); % reformat strings to be usable in titles and legends

%     %% Correlation of the different behaviours with each other
%     figure();imagesc(corr(all_beh'))
% 
%     %% Correlation of the behaviours values per event
%     figure();imagesc(corr(all_beh))

    %Soma_ROIs   = find(ismember(1:size(data, 1), obj.ref.indices.somatic_ROIs)); % somatic ROIs, but ignoring the NaNs
    All_ROIs    = 1:size(data, 1);
    out = train_and_test(data, processed_behaviours, timepoints, All_ROIs, method, behaviour_list, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs),2), parameters);
end



