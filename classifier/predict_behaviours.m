%% Classify behaviour for a given experiment, for all behaviours
%% use expe 07-26 to fix gaps in behaviour

function out = predict_behaviours(obj, use_classifier, method, type_of_trace)
    if nargin < 1 || isempty(obj)
        obj = '' 
    end
    if nargin < 2 || isempty(use_classifier)
        use_classifier = true; % if false, use regression learner
    end
    if nargin < 3 || isempty(method)
        method = 'svm'; 
    end
    if nargin < 4 || isempty(type_of_trace)
        type_of_trace   = 'subtracted_peaks'; % ['subtracted' OR 'rescaled' OR 'raw'] AND ['peaks' or '']. eg 'subtracted_peaks' , or 'raw'
    end

    use_hd_data     = false;
    time_filter     = 0;

    %% Make sure preprocessing was done correctly
    rendering       = obj.rendering;
    obj.rendering   = false;
    [obj, source_signal, ~, timepoints] = prepare_phate_analysis(obj, use_hd_data, time_filter, type_of_trace);
    obj.rendering   = rendering;


    %% Get the signal for the selected timepoints and ROIs
    Valid_ROIs      = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); % remove excluded branches AND bad_ROIs based on correlation)
    Ca              = source_signal(timepoints, Valid_ROIs)';

    %     %% Correlation of the different ROIs with each other
    %     figure();imagesc(corr(Ca'))
    % 
    %     %% Correlation of the Ca2+ values per event
    %     figure();imagesc(corr(Ca))

    %% Get all behaviours at these tp
    behaviours              = [];
    raw_behaviours          = [];
    processed_behaviours    = [];
    beh_types               = obj.behaviours.types;
    for type = beh_types
        %% Get original bhaviour
        [~,~,original_beh]   = obj.get_behaviours(type{1},'',true);
        
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
            processed_behaviours= logical([processed_behaviours; current_beh(timepoints)]);
        end
    end
    beh_types   = strrep(beh_types, '_', '\_'); % reformat strings to be usable in titles and legends

%     %% Correlation of the different behaviours with each other
%     figure();imagesc(corr(all_beh'))
% 
%     %% Correlation of the behaviours values per event
%     figure();imagesc(corr(all_beh))

    %Soma_ROIs   = find(ismember(ROI_to_study, obj.ref.indices.somatic_ROIs)); % somatic ROIs, but ignoring the NaNs
    All_ROIs    = 1:numel(Valid_ROIs);
    out = train_and_test(Ca, processed_behaviours, timepoints, All_ROIs, method, beh_types, raw_behaviours, nanmedian(obj.rescaled_traces(:,All_ROIs),2));
end




% 
% 
% 
% 
% 
% %% Define analysis otpions
% variable    = 'binary';% 'binary'
% ref_roi     = 'all'; % 'soma' or 'all'
% method      = 'svm' ; % 'svm' or 'forest'
% 
% %% Set ROIs to use
% if strcmp(ref_roi, 'soma')
%     roi_subset  = Soma_ROIs; % this is the one to study
% else
%     roi_subset  = All_ROIs; % this is the one to study
% end
% cc          = [];
% for beh_idx = 1:numel(beh_types)  
%     %% Set training and testing dataset
%     train_R     = rand_tp(1:floor(numel(peak_tp)/2));
%     test_R      = rand_tp(ceil(numel(peak_tp)/2):end);
%     x_train     = double(Ca(roi_subset,train_R))';
%     x_test      = double(Ca(roi_subset,test_R))';
%     
%     %% Set variable to predict
%     current_beh = all_beh(beh_idx,:);
%     
%     %% Run training
%     if strcmp(variable, 'binary')
%         thr = nanmin(current_beh) + range(current_beh) / 5; % 20 %% of min        
%         y_train = all_beh(beh_idx,train_R)' > thr;    
%         y_test = all_beh(beh_idx,test_R)' > thr;
%         
%         if strcmp(method, 'svm')
%             %% SVM classifier
%             MODEL = fitcsvm(x_train,y_train);
%         elseif strcmp(method, 'forest')
%             MODEL = fitcensemble(x_train,y_train);
%         end
%     else
%         y_train = all_beh(beh_idx,train_R)';    
%         y_test = all_beh(beh_idx,test_R)';
%         %% multiclass SVM classifier        
%         MODEL=fitcecoc(x_train,y_train); % fitcecoc
%     end
% 
%     y_pred = predict(MODEL, x_test);
%     figure(beh_idx);clf();plot(y_pred); hold on; plot(y_test); legend({'prediction','measurement'});title(beh_types{beh_idx})
%     if all(y_pred)
%         y_pred(1) = 0;
%     elseif ~any(y_pred)
%         y_pred(1) = 1;
%     end
%     cc(beh_idx) = corr(y_pred, y_test)
% end
% figure();bar(categorical(beh_types), cc)

% if strcmp(method, 'forest')
%     impOOB = oobPermutedPredictorImportance(MODEL);
%     obj.ref.plot_value_tree(impOOB, ROI_to_study)
% end

% 
% 
% 
% 
%     x_train = double(Ca)';
%     y_train = all_beh(1,:)';
% 
% 
%     t = templateTree('NumVariablesToSample','all',...
%         'PredictorSelection','interaction-curvature','Surrogate','on');
%     rng(1); % For reproducibility
%     Mdl = fitrensemble(x_train,y_train,'Method','Bag','NumLearningCycles',200, ...
%         'Learners',t);
% 
%     yHat = oobPredict(Mdl);
%     R2 = corr(Mdl.Y,yHat)^2
% 
%     impOOB = oobPermutedPredictorImportance(Mdl);
% 
%     figure
%     bar(impOOB)
%     title('Unbiased Predictor Importance Estimates')
%     xlabel('Predictor variable')
%     ylabel('Importance')
%     h = gca;
%     h.XTickLabel = Mdl.PredictorNames;
%     h.XTickLabelRotation = 45;
%     h.TickLabelInterpreter = 'none';
% 
%     
%     [impGain,predAssociation] = predictorImportance(Mdl);
% 
%     figure
%     plot(1:numel(Mdl.PredictorNames),[impOOB' impGain'])
%     title('Predictor Importance Estimation Comparison')
%     xlabel('Predictor variable')
%     ylabel('Importance')
%     h = gca;
%     h.XTickLabel = Mdl.PredictorNames;
%     h.XTickLabelRotation = 45;
%     h.TickLabelInterpreter = 'none';
%     legend('OOB permuted','MSE improvement')
%     grid on
%     
%     obj.ref.plot_value_tree(impOOB, find(obj.dimensionality.valid_trace_idx))
% 
% 
% 
% forest = fitrensemble(zscore(x_train), ... % all features as input but max firing rate
%     zscore(y_train), 'crossval', 'on', 'kfold', 10); %%% implement crossvalidation during at training of the forest
% 
% predicted_firingrate = kfoldPredict(forest); % perform crossvalidated regression
% 
% figure; 
% plot(zscore(y_train), predicted_firingrate, 'o'); % scatter plot of actual vs predicted maximal firing rates
% axis square
% xlabel 'Actual zscore of max firing rate'
% ylabel 'Regressed zscore of max firing rate'
% title 'Random forest regression'
