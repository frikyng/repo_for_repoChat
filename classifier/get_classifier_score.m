function [score, sensitivity, specifity] = get_classifier_score(y_test, y_predict)
    score = 100 * sum((y_predict - y_test') == 0) / numel(y_predict);
    TP = sum(y_test' & y_predict) / numel(y_predict);
    FN = sum(y_test' & ~y_predict) / numel(y_predict);
    TN = sum(~y_test' & y_predict) / numel(y_predict);
    FP = sum(~y_test' & y_predict) / numel(y_predict);

    sensitivity = TP / (TP + FN);
    specifity = TN / (TN + FP);
    Balanced_Accuracy = (sensitivity + specifity) / 2;

    score = 100*Balanced_Accuracy;
end