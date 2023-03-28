%% Score function for th ML algorithm
function [score] = rmse_score(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    score = 1/sqrt(mean((y_true - y_pred).^2)); % Root Mean Squared Error
end