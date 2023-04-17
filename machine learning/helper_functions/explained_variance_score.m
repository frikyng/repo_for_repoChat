%% Score function for th ML algorithm
function [score] = explained_variance_score(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    R = corrcoef(y_pred, y_true);
    score = (1 - (1 - R(1,2)^2)) * 100;
end