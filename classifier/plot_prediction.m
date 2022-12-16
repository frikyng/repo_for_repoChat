function plot_prediction(raw_data, YData, timepoints, partition, y_predict, raw_behaviour, behaviour_name, score, Fig_nb)
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
        Fig_nb = next_fig;
    end
    
    %% Split data between train and test
    if islogical(partition)
        train_tp    = find(partition);
        test_tp     = find(~partition);
    else
        train_tp    = find(training(partition));
        test_tp     = find(test(partition));
    end
    y_train     = YData(train_tp);              % Variable to predict used for training  
    y_test      = YData(test_tp);               % Variable to predict used for testing (ground truth) 
        
    figure(Fig_nb);clf();
    if numel(score) == 1 || ~islogical(YData)
        sgtitle([behaviour_name, ' accuracy is ',num2str(score(1), 3),'%']);hold on;
    else
        sgtitle({[behaviour_name, ' accuracy is ',num2str(score(1), 3),'%'],['specificity = ',num2str(round(score(2))),' %; sensitivity = ',num2str(round(score(3))),' %; MCC = ',num2str(round(score(4))),' %']});hold on;
    end
    ax1 = subplot(3,1,1); hold on;       
    plot(raw_data);hold on
    s1 = scatter(timepoints(train_tp),raw_data(timepoints(train_tp)), 'k', 'filled');hold on;
    s2 = scatter(timepoints(test_tp),raw_data(timepoints(test_tp)), 'b', 'filled');legend([s1, s2] ,{'train','test'})

    ax2 = subplot(3,1,2); 
    if islogical(YData)
        plot(timepoints(train_tp), y_train, 'ko'); hold on;
        plot(timepoints(test_tp), test_tp, 'go'); hold on;
        wrong = (y_predict - y_test') ~= 0;
        plot(timepoints(test_tp(~wrong)), y_predict(~wrong),'go', 'MarkerFaceColor', 'g'); hold on        
        plot(timepoints(test_tp(wrong)), y_predict(wrong),'ro', 'MarkerFaceColor', 'r'); hold on
        ylim([-0.05, 1.05]);
    else
        plot(timepoints(train_tp), normalize(y_train, 'medianiqr')); hold on;
        plot(timepoints(test_tp), normalize(y_predict, 'medianiqr'),'g'); hold on
        %plot(timepoints(test_tp), normalize(smoothdata(y_predict, 'gaussian' , 10), 'medianiqr'),'b--')
    end

    ax3 = subplot(3,1,3); 
    plot(raw_behaviour);hold on
    if islogical(YData)
        scatter(timepoints(train_tp(y_train)),raw_behaviour(timepoints(train_tp(y_train))), 'k^', 'filled');hold on;
        scatter(timepoints(train_tp(~y_train)),raw_behaviour(timepoints(train_tp(~y_train))), 'kv', 'filled');hold on;
        scatter(timepoints(test_tp(y_test)),raw_behaviour(timepoints(test_tp(y_test))), 'g^', 'filled');hold on;
        scatter(timepoints(test_tp(~y_test)),raw_behaviour(timepoints(test_tp(~y_test))), 'gv', 'filled');hold on;
        errors_pos = (y_predict - YData(test_tp)') ~= 0 & YData(test_tp)';
        errors_neg = (y_predict - YData(test_tp)') ~= 0 & ~YData(test_tp)';
        scatter(timepoints(test_tp(errors_pos)),raw_behaviour(timepoints(test_tp(errors_pos))), 'r^', 'filled');hold on;
        scatter(timepoints(test_tp(errors_neg)),raw_behaviour(timepoints(test_tp(errors_neg))), 'rv', 'filled');hold on;
    else
        scatter(timepoints(train_tp),raw_behaviour(timepoints(train_tp)), 'ko', 'filled');hold on;
        plot([timepoints(test_tp);timepoints(test_tp)], [y_test; y_predict'], 'go-')
        scatter(timepoints(test_tp),y_predict, ' ro', 'filled');hold on;
    end
    

    linkaxes([ax1, ax2, ax3],'x')
    drawnow
end