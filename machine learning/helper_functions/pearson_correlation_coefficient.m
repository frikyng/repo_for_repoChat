%% Score function for th ML algorithm
function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
    if isrow(y_pred)
        y_pred = y_pred';
    end
    if isrow(y_true)
        y_true = y_true';
    end
    score = corr(y_true, y_pred); % pearson r
end