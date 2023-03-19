function fig_handle = show_summary_predictive_scores(all_results, labels)
    labels = fix_labels(labels);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    N_phates = size(results, 1);
    
    fig_handle = figure(3);clf();
    
    subplot(1,2,1);
    
    imagesc(results ./ max(results));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility per behaviour)'};
    xticks(1:numel(labels));xticklabels(cellfun(@(x) strrep(x,'_',' '), labels, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');hold on;
    yticks(1:2:N_phates*2);hold on;yticklabels(cellfun(@(x) num2str(x), num2cell(1:N_phates), 'uni', false));
    ylabel('Phate #'); hold on;

    subplot(1,2,2);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    imagesc(results ./ max(results, [],2));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility in phate)'};
    xticks(1:numel(labels));xticklabels(cellfun(@(x) strrep(x,'_',' '), labels, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');
    yticks(1:2:N_phates*2);hold on;yticklabels(cellfun(@(x) num2str(x), num2cell(1:N_phates), 'uni', false));
    ylabel('Phate #');
    
    set(fig_handle, 'Color', 'w'); hold on;
end