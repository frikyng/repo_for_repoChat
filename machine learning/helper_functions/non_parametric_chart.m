function non_parametric_chart(values, input_labels, errors, paired)
    if nargin < 3 || isempty(errors)
        errors = max(values);
    end
    if nargin < 4 || isempty(paired)
        paired = false;
    end

    if ~iscategorical(input_labels)
        label   = categorical(input_labels);
        label   = reordercats(categorical(label),input_labels); % preserve input label order
    else
        label   = input_labels;
    end
    n_samples = size(values, 2);
    figure();
    bar(label,nanmean(values'));hold on;ylim([-10,100]);
    errorbar(label,nanmean(values'),std(values')/sqrt(n_samples),'k.','CapSize',0);
    hold on; scatter(label,values,'ko', 'filled', 'MarkerFaceAlpha',0.3)
    set(gcf, 'Color', 'w');set(gca,'box','off'); ylabel('Explained Variance');
    stats_results = {};
    disp('some conditions are significantly different from others')
    if ~paired
        hold on; scatter(label,values,'ko', 'filled', 'MarkerFaceAlpha',0.3)
        %[stats_results.p,stats_results.table,stats_results.stats]   = kruskalwallis(values', label, 'off');
        [stats_results.p,stats_results.table,stats_results.stats]   = anova1(values', label, 'off');
        stats_results.multi_comp                                    = multcompare(stats_results.stats,  'Display', 'off');
    else
%         hold on; plot(label,values,'-','Color',[0.7,0.7,0.7]);
%         T = array2table(values', 'VariableNames', {'Group1', 'Group2', 'Group3', 'Group4', 'Group5'});
%         within = table([1 2 3 4 5]','VariableNames',{'Group'});
%         rm = fitrm(T, 'Group1-Group5~1', 'WithinDesign', within);
%         [stats_results.p,stats_results.table,stats_results.stats]   = ranova(rm, 'WithinModel', 'Group');
%         stats_results.multi_comp = multcompare(rm, 'Group', 'ComparisonType', 'bonferroni')
%         stats_results.multi_comp = [stats_results.multi_comp.Group_1, stats_results.multi_comp.Group_2,stats_results.multi_comp.Group_2,stats_results.multi_comp.Group_2,stats_results.multi_comp.Group_2,stats_results.multi_comp.pValue]

        [stats_results.p,stats_results.table,stats_results.stats]   = friedman(values', 1, 'off');
        stats_results.multi_comp                                    = multcompare(stats_results.stats,  'Display', 'off');
    end
    stats_results.multi_comp(stats_results.multi_comp(:,6) == 0,6) = Inf
    %

    ns          = stats_results.multi_comp(stats_results.multi_comp(:,6) >= 0.05, [1,2]);
    p0_01       = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.05 & stats_results.multi_comp(:,6) >= 0.01, [1,2,6]);
    p0_001      = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.01 & stats_results.multi_comp(:,6) >= 0.001, [1,2,6]);
    p_strong    = stats_results.multi_comp(stats_results.multi_comp(:,6) < 0.001, [1,2,6]);

    [M,loc] = max(nanmean(values'));
    M = M + errors(loc);
    step = (100 - M)/ (sum(stats_results.multi_comp(:,6) < 0.05)+1);
    for el = 1:size(p0_01, 1)
        M = M+step;
        overbar(p0_01(el,1) ,p0_01(el,2), M,'*'); %  ['p=',num2str(round(p0_01(el,3),3))]
    end
    for el = 1:size(p0_001, 1)
        M = M+step;
        overbar(p0_001(el,1) ,p0_001(el,2), M, '**'); %['p=',num2str(round(p0_001(el,3),3))]
    end
    for el = 1:size(p_strong, 1)
        M = M+step;
        overbar(p_strong(el,1) ,p_strong(el,2), M, '***'); %  ['p=',num2str(round(p_strong(el,3),3))]
    end
end
