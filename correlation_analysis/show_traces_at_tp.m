function show_traces_at_tp(data_sm, tp, mean_trace, key)
    if nargin < 4 || isempty(key)
        key = 1:size(data_sm, 2);
    end
    
    %saveas(1081,['test_',num2str(tp),'.png'])
    figure(128);cla();plot(mean_trace,'k','LineWidth',2);hold on;plot(data_sm(:,key),'r','LineWidth',2);hold on;
end