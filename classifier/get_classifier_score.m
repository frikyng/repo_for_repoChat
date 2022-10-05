function [score, TPR, TNR, MCC] = get_classifier_score(y_test, y_predict)
    y_test = y_test';
    
    if islogical(y_test)
        %% See https://en.wikipedia.org/wiki/Sensitivity_and_specificity
        TP = sum(y_test & y_predict) / numel(y_predict);
        FN = sum(y_test & ~y_predict) / numel(y_predict);
        TN = sum(~y_test & ~y_predict) / numel(y_predict);
        FP = sum(~y_test & y_predict) / numel(y_predict);

        TPR = 100* TP / (TP + FN); % aka Recall or sensitivity: # of correct positive predictions out of all positive predictions (TP and FN) that could be made
        TNR = 100 * TN / (TN + FP); % aka Selectivity or specificity: # of correct negative predictions over the total negative prediction made 
        Balanced_Accuracy = (TPR + TNR) / 2; % use case when dealing with imbalanced data

        %     
        %     PP = TP + FP;
        %     PN = TN + FN;
        %     P = TP + FN;
        %     N = FP + TN;
        %     
        %     PPV = TP / PP;
        %     NPV = TN / PN;
        %     FNR = FN / P;
        %     FPR = FP / N;
        %     FOR = FN / PN;
        %     FDR = FP / PP;
        %% https://en.wikipedia.org/wiki/Phi_coefficient
        MCC = 100*(TP*TN - FP*FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
        score = Balanced_Accuracy;
    else
        if numel(unique(y_predict))  == 1
            warning(['Prediction returned ',num2str(unique(y_predict)),' for all prediction. It was replaced by a random set of values to get correlation score '])
            y_predict = rand(size(y_predict)) * range(y_test);
        end
        [score, TPR, TNR, MCC] = deal(corr(y_test,y_predict)*100);
    end
end