function fig_handle = show_summary_predictive_scores(all_results, x_ticks, y_ticks, y_label)
    x_ticks = fix_labels(x_ticks);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    N_phates = size(results, 1);

    if nargin < 3 || isempty(y_ticks)
        y_ticks = cellfun(@(x) num2str(x), num2cell(1:N_phates), 'uni', false);
        tick_spacing = 1:2:N_phates*2;
    else        
        tick_spacing = 1:numel(y_ticks);
    end
    if nargin < 4 || isempty(y_label)
        y_label = 'Phate #';
    end
    
    

    fig_handle = figure(3);clf();
    
    subplot(1,2,1);
    
    imagesc(results ./ max(results));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility per behaviour)'};
    xticks(1:numel(x_ticks));xticklabels(cellfun(@(x) strrep(x,'_',' '), x_ticks, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');hold on;
    yticks(tick_spacing);hold on;yticklabels(y_ticks);
    ylabel(y_label); hold on;

    subplot(1,2,2);
    results = cell2mat(all_results)';%a(a < 20) = 0;
    imagesc(results ./ max(results, [],2));caxis([0,1]);colormap(viridis);
    cb = colorbar;cb.Label.String = {'Normalized predictive score ','(normed to max predictibility in phate)'};
    xticks(1:numel(x_ticks));xticklabels(cellfun(@(x) strrep(x,'_',' '), x_ticks, 'uni', false));xtickangle(45); hold on;
    xlabel('Behaviour');
    yticks(tick_spacing);hold on;yticklabels(y_ticks);
    ylabel(y_label);
    
    set(fig_handle, 'Color', 'w'); hold on;
end