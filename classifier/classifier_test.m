%% use expe 07-26 to fix gaps in behaviour
if ~exist('obj','var')
    obj = arboreal_scan_experiment('', false);

    obj.rescaling_method = 'peaks_trials';
    obj.prepare_binning({'distance',50});
    obj.find_events();
    obj.rescale_traces();
    obj.set_median_traces();
    obj.compute_similarity();
end

%% Get tp of Calcium peaks
peak_tp = find(obj.get_tp_for_condition('peaks'));

%% Get Calcium singals, median subtracted at these tp
Ca              = obj.rescaled_traces(peak_tp, :)';
Ca              = Ca - nanmedian(Ca, 1);
ROI_to_study    = obj.ref.indices.valid_swc_rois(~ismember(obj.ref.indices.valid_swc_rois, obj.bad_ROI_list)); % remove excluded branches AND bad_ROIs based on correlation)
Ca              = Ca(ROI_to_study, :);

%% Correlation of the different ROIs with each other
figure();imagesc(corr(Ca'))

%% Correlation of the Ca2+ values per event
figure();imagesc(corr(Ca))

%% Get all behaviours at these tp
all_beh     = [];
full_beh    = [];
beh_types   = obj.behaviours.types;
for type = beh_types
    [~,~,beh]   = obj.get_behaviours(type{1},'',true);
    beh         = smoothdata(beh.value, 'gaussian', 50);
    full_beh    = [full_beh; beh];
    beh         = beh(peak_tp);
    all_beh     = [all_beh; beh];
end
beh_types   = strrep(beh_types, '_', '\_'); % reformat strings to be usable in titles and legends

%% Normalize Behaviours (TBD)
all_beh = [normalize(all_beh(1:end-1,:), 2,'medianiqr') ; normalize(all_beh(end, :), 2,'zscore')];


thr = min(all_beh,[],2) + ((prctile(all_beh,5, 2) - min(all_beh,[],2))*2)
all_beh_bin = all_beh > thr;

%% Correlation of the different behaviours with each other
figure();imagesc(corr(all_beh'))

%% Correlation of the behaviours values per event
figure();imagesc(corr(all_beh))

%% Get random timepoints for traing/testing
rand_tp     = randperm(numel(peak_tp));

%% Define the Ca2_ ROIs we want to use
All_ROIs    = 1:numel(ROI_to_study);
Soma_ROIs   = find(ismember(ROI_to_study, obj.ref.indices.somatic_ROIs)); % somatic ROIs, but ignoring the NaNs

roi_subset = cl_4; %45:5:70 in 2019-09-26_exp_1
for el = 1:numel(beh_types)
    type = strrep(beh_types{el},' ','_');
    type_corrected = strrep(type,'\_','_');
    current_var = all_beh_bin(el,:);
    eval(['Bvar',type_corrected,'=all_beh_bin(',num2str(el),',:);'])
    eval(['var',type_corrected,'=all_beh(',num2str(el),',:);'])
    %[trainedClassifier, validationAccuracy] = trainClassifier(Ca, eval(['var',type]), var_to_use);

    train_R     = rand_tp(1:floor(numel(peak_tp)/2)); 

    if contains(type, 'trigger')
        cost = [0,1;10,0];
    else
        cost = [0,10;1,0];
    end        
    %phate(double(Ca(roi_subset,:)'),'ndim',10)'
    [y_predict, y_test, score, x_test, x_train, y_train] = csvm_prediction(Ca(roi_subset,:), current_var, train_R, cost);
    
%     all_scores = [];
%     for roi = 1:size(Ca, 1)
%         roi
%         temp_Ca = Ca;
%         temp_Ca(roi,:) = [];
%         [~, ~, all_scores(roi), ~, ~, ~] = csvm_prediction(temp_Ca, current_var, train_R, cost);        
%     end

    plot_prediction(nanmedian(obj.rescaled_traces(:,roi_subset),2), Ca(roi_subset,:), current_var, peak_tp, train_R, y_predict, full_beh(el,:), type, score, el);
end

arrangefigures
%figure(1);clf();bar(categorical(beh_types), accuracy)















%% Define analysis otpions
variable    = 'binary';% 'binary'
ref_roi     = 'all'; % 'soma' or 'all'
method      = 'svm' ; % 'svm' or 'forest'

%% Set ROIs to use
if strcmp(ref_roi, 'soma')
    roi_subset  = Soma_ROIs; % this is the one to study
else
    roi_subset  = All_ROIs; % this is the one to study
end
cc          = [];
for beh_idx = 1:numel(beh_types)  
    %% Set training and testing dataset
    train_R     = rand_tp(1:floor(numel(peak_tp)/2));
    test_R      = rand_tp(ceil(numel(peak_tp)/2):end);
    x_train     = double(Ca(roi_subset,train_R))';
    x_test      = double(Ca(roi_subset,test_R))';
    
    %% Set variable to predict
    current_beh = all_beh(beh_idx,:);
    
    %% Run training
    if strcmp(variable, 'binary')
        thr = nanmin(current_beh) + range(current_beh) / 5; % 20 %% of min        
        y_train = all_beh(beh_idx,train_R)' > thr;    
        y_test = all_beh(beh_idx,test_R)' > thr;
        
        if strcmp(method, 'svm')
            %% SVM classifier
            MODEL = fitcsvm(x_train,y_train);
        elseif strcmp(method, 'forest')
            MODEL = fitcensemble(x_train,y_train);
        end
    else
        y_train = all_beh(beh_idx,train_R)';    
        y_test = all_beh(beh_idx,test_R)';
        %% multiclass SVM classifier        
        MODEL=fitcecoc(x_train,y_train); % fitcecoc
    end

    y_pred = predict(MODEL, x_test);
    figure(beh_idx);clf();plot(y_pred); hold on; plot(y_test); legend({'prediction','measurement'});title(beh_types{beh_idx})
    if all(y_pred)
        y_pred(1) = 0;
    elseif ~any(y_pred)
        y_pred(1) = 1;
    end
    cc(beh_idx) = corr(y_pred, y_test)
end
figure();bar(categorical(beh_types), cc)

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








