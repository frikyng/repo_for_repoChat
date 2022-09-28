function plot_prediction(raw_data, YData, timepoints, train_tp, y_predict, raw_behaviour, behaviour_name, score, Fig_nb)
    if nargin < 6 || isempty(raw_behaviour) || all(isnan(raw_behaviour))
        raw_behaviour = zeros(1, numel(raw_data));
    end
    if nargin < 7 || isempty(behaviour_name)
        behaviour_name = 'Unnamed behaviour';
    end
    if nargin < 8 || isempty(score)
        score = NaN;
    end
    if nargin < 9 || isempty(Fig_nb)
        existing = sort(arrayfun(@(x) x.Number, findobj('type','figure')));
        Fig_nb = existing(find(diff(existing) > 1, 1, 'first')) + 1;        
    end
    
    train_tp    = sort(train_tp);
    all_tp      = 1:numel(timepoints);
    test_tp     = all_tp(~ismember(all_tp, train_tp));
    %     x_train     = double(XData(:,train_tp))';   % Predictors for training
    %     x_test      = double(XData(:,test_tp))';    % Predictors for testing
    y_train     = logical(YData(train_tp));              % Variable to predict used for training  
    y_test      = logical(YData(test_tp));               % Variable to predict used for testing (ground truth)     
    
        
    figure(Fig_nb);clf();sgtitle([behaviour_name, ' accuracy is ',num2str(score, 3),'%']);hold on;
    ax1 = subplot(2,1,1); hold on;       
    plot(raw_data);hold on
    scatter(timepoints(train_tp),raw_data(timepoints(train_tp)), 'k', 'filled');hold on;
    scatter(timepoints(test_tp),raw_data(timepoints(test_tp)), 'g', 'filled');

    ax2 = subplot(2,1,2); 
    plot(raw_behaviour);hold on
    scatter(timepoints(train_tp(y_train)),raw_behaviour(timepoints(train_tp(y_train))), 'k^', 'filled');hold on;
    scatter(timepoints(train_tp(~y_train)),raw_behaviour(timepoints(train_tp(~y_train))), 'kv', 'filled');hold on;
    scatter(timepoints(test_tp(y_test)),raw_behaviour(timepoints(test_tp(y_test))), 'g^', 'filled');hold on;
    scatter(timepoints(test_tp(~y_test)),raw_behaviour(timepoints(test_tp(~y_test))), 'gv', 'filled');hold on;
    errors_pos = (y_predict - YData(test_tp)') ~= 0 & YData(test_tp)';
    errors_neg = (y_predict - YData(test_tp)') ~= 0 & ~YData(test_tp)';
    scatter(timepoints(test_tp(errors_pos)),raw_behaviour(timepoints(test_tp(errors_pos))), 'r^', 'filled');hold on;
    scatter(timepoints(test_tp(errors_neg)),raw_behaviour(timepoints(test_tp(errors_neg))), 'rv', 'filled');hold on;
    
    linkaxes([ax1, ax2],'x')
    drawnow
end