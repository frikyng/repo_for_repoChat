%% Classify behaviour for a given experiment, for all behaviours
%% use expe 07-26 to fix gaps in behaviour

%% Cell 2019-09-17_exp_1 anticorelated to running

function [out, data, ROI_groups, meanvalue] = predict_behaviours(obj, use_classifier, method, type_of_trace, behaviour_list, ROI_groups, n_iter,rand_ROI_groups, varargin)
    if nargin < 1 || isempty(obj)
        obj = ''; 
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
        ROI_groups   = num2cell(ROI_groups);    % does this force groups to be cell array regardless of input type?
    end
    if iscolumn(ROI_groups)
        ROI_groups = ROI_groups';
    end
    if nargin < 7 || isempty(n_iter)
        n_iter   = 1; 
    end
    if nargin < 8 || isempty(rand_ROI_groups)
        rand_ROI_groups   = false; % false is no randomization, true is randperm (this is to compare info content in a group of (e.g. PHATE discovered) ROIs with randomly selected ROIs (from anywhere from the tree), esuring the same # of ROIs in each)
    end
    if nargin < 9 || isempty(varargin)
        ml_parameters = DEFAULT_CLASSIFIER_OPTION;
    else
        ml_parameters = DEFAULT_CLASSIFIER_OPTION(varargin);
    end
    
    use_hd_data             = false;
    time_filter             = 0;
    detrend_behaviours      = true;
    smooth_behaviours       = 1;
    pt_per_s                = 1/nanmedian(obj.timescale.sr);

    %% Make sure preprocessing was done correctly
    rendering       = obj.rendering;
    obj.rendering   = false;
    
    if ischar(type_of_trace) && contains(type_of_trace, '_lag')
        % Regular expression to extract number after "lag"
        expr = '(?<=lag)[+-]?\d+';
        % Extract number from each string
        lag = str2num(regexp(type_of_trace, expr, 'match', 'once'));
        type_of_trace = erase(type_of_trace,{['_lag',num2str(lag)],['_lag_',num2str(lag)]});
    else
        lag = 0;
    end

    [obj, source_signal, ~, timepoints] = prepare_phate_analysis(obj, use_hd_data, time_filter, type_of_trace);

    beh_timepoints                              = timepoints + lag; % if lag == 0 then points are the same
    timepoints(beh_timepoints < 1)              = [];
    beh_timepoints(beh_timepoints < 1)          = [];
    timepoints(beh_timepoints > obj.tp)         = [];
    beh_timepoints(beh_timepoints > obj.tp)     = [];
    
    % holder = source_signal(:,1); figure;plot(source_signal(:,1));hold on; scatter(timepoints, holder(timepoints),'r','filled')
    % figure();hist(reshape(source_signal(timepoints,:),[],1),100); % TOMMY uncomment to see distribution of predictors
    obj.rendering   = rendering;
    
    %% Get all behaviours at these tp
    behaviours              = [];
    raw_behaviours          = [];
    processed_behaviours    = [];
    
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
        if contains(type, 'trigger')
            current_beh         = smoothdata(original_beh.value, 'gaussian', [smooth_behaviours*pt_per_s/nanmedian(diff(obj.t)),0]); 
        elseif contains(type, 'baseline')
            current_beh         = original_beh.value;
        else
            current_beh         = smoothdata(original_beh.value, 'movmedian', smooth_behaviours*pt_per_s);
            
            %% Gaussian smoothing to denoise behaviour
            current_beh         = smoothdata(current_beh, 'gaussian', smooth_behaviours*pt_per_s*5); 
            
            current_beh(current_beh == 0) = randn(sum(current_beh == 0),1) * (rms(current_beh)/100);
        end

        
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
        behaviours          = [behaviours; current_beh(beh_timepoints)];
        if use_classifier
            processed_behaviours= logical([processed_behaviours; current_beh(beh_timepoints) > thr]);
        else
            processed_behaviours= [processed_behaviours; current_beh(beh_timepoints)];
        end
    end
    behaviour_list   = strrep(behaviour_list, '_', '\_'); % reformat strings to be usable in titles and legends

    
    %% Get the signal for the selected timepoints and ROIs
    %Valid_ROIs      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); % remove excluded branches AND bad_ROIs based on correlation)
    invalid_ROIs_logical    = ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list);
    invalid_ROIs            = obj.ref.indices.valid_swc_rois(invalid_ROIs_logical); % remove excluded branches AND bad_ROIs based on correlation)
    source_signal(:, invalid_ROIs) = NaN;
    data                    = NaN(numel(ROI_groups), numel(timepoints));
    for gp_idx = 1:numel(ROI_groups)
        ROI_groups{gp_idx}(ismember(ROI_groups{gp_idx}, obj.bad_ROI_list)) = [];
        data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups{gp_idx}),2)';        
    end

    fully_invalid_group = cellfun(@isempty, ROI_groups)' | all(isnan(data),2);
    ROI_groups(fully_invalid_group) = [];
    data(fully_invalid_group, :) = [];

    %     %% Correlation of the different ROIs with each other
    %     figure();imagesc(corr(data'))
    % 
    %     %% Correlation of the Ca2+ values per event
    %     figure();imagesc(corr(data))


    %     %% Correlation of the different behaviours with each other
    %     figure();imagesc(corr(processed_behaviours'))
    % 
    %     %% Correlation of the behaviours values per event
    %     figure();imagesc(corr(processed_behaviours))

    %     Soma_ROIs   = find(ismember(1:size(data, 1), obj.ref.indices.somatic_ROIs)); % somatic ROIs, but ignoring the NaNs

    
    %% FYI find(all(isnan(source_signal))) should be empty, or you have observation with only NaN
    
    All_ROIs    = 1:size(data, 1);

    if ~rand_ROI_groups
        for iter = 1:n_iter
            out{iter}         = train_and_test(data, processed_behaviours, timepoints, All_ROIs, method, behaviour_list, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
        end
    else  
        % now, as long as you don't have all the ROIs in ROI_groups to
        % start with, we build alternative randomized groups
        if ~all(cellfun(@numel, ROI_groups) == 1) && numel(unique(vertcat(ROI_groups{:}))) ~= sum(~fully_invalid_group)
            rand_ROI_pool      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, vertcat(ROI_groups{:})) & ~ismember(obj.ref.indices.valid_swc_rois, invalid_ROIs)); %list of valid ROIs excluding ROI groups
        else % when you pass all ROIs (so all group size == 1), randomization for leftover ROIs makes no sense so we just randomize ROIs
            rand_ROI_pool      = cell2mat(ROI_groups);
        end
        for iter = 1:n_iter
            try
                ROI_groups_rdm     = cellfun(@(x) rand_ROI_pool(randperm(numel(rand_ROI_pool),x)), cellfun(@numel, ROI_groups, 'uni', false), 'UniformOutput', false); % size matched groups from ROI_pool
            catch
                warning('Not enough non-used ROIs available. Using all ROIs instead for randomization')
                ok = obj.ref.indices.valid_swc_rois(~invalid_ROIs_logical);
                ROI_groups_rdm     = cellfun(@(x) ok(randperm(numel(ok),x)), cellfun(@numel, ROI_groups, 'uni', false), 'UniformOutput', false); % size matched groups from ROI_pool
            end
            for gp_idx = 1:numel(ROI_groups_rdm)
                data(gp_idx, :) =  nanmean(source_signal(timepoints, ROI_groups_rdm{gp_idx}),2)';        
            end                
            out{iter}         = train_and_test(data, processed_behaviours, timepoints, 1:numel(ROI_groups_rdm), method, behaviour_list, raw_behaviours, nanmedian(obj.rescaled_traces(:,~invalid_ROIs_logical),2), ml_parameters);
        end
    end

    meanvalue = [];
    if ml_parameters.rendering
        meanvalue = bar_chart(out, 'beh_type');
        title(ml_parameters.title)
    end
end
