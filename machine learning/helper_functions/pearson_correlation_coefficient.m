%% Score function for th ML algorithm
function [score] = pearson_correlation_coefficient(y_true, y_pred, w)
    score = corr(y_true, y_pred)*100;
end