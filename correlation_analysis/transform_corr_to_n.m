function [n_corr, tp_active] = transform_corr_to_n(data_sm, comb, corr_results, tp, idx, thr)
    n_corr = [];
    for key = 1:size(data_sm,2)
        temp = corr_results(tp, comb(:,2) == key | comb(:,1) == key);
        %% Add the key itself
        temp = [temp(1,1:(key-1)), NaN, temp(1,key:end)];
        temp = temp > thr;
        n_corr(key) = sum(temp);
    end
    tp_active = mean(n_corr);
    n_corr = n_corr(idx);
end
