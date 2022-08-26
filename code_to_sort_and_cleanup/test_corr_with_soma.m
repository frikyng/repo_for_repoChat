% sucess_rate = {};
% mean_corr = {};
% all_corr = {};
% t_somatic = {};
% 
% 

for expe_idx = find(cell2mat(depth) > 800)%numel(a.experiments):-1:1%numel(a.experiments):-1:1
    rendering     = false;
    current_expe  = a.experiments(expe_idx);
    idx           = current_expe.ref.indices.ROIs_list; % will duplicate the start of branches

    %% Get signal
    data        = current_expe.rescaled_traces;
    t           = current_expe.event_fitting.peak_pos;
    mean_trace  = nanmedian(data, 2);
    %t          = t(mean_trace(t) > 50); % use to clip/kep baps
    data_sm     = smoothdata(data, 'gaussian', [20, 0]);
    
    
    data_sm(:,ROIs_to_EXCLUDE{expe_idx}) = NaN;
    
    soma = current_expe.ref.indices.somatic_ROIs;
    if numel(soma) == 1
        soma = soma:soma+3;
    end
    data_sm(:, soma) = repmat(nanmedian(data_sm(:, soma), 2),1, numel(soma));

    
    [~, peak_t_short] = detect_events(nanmean(data_sm(:, soma),2),'',30);
    peak_t_short = peak_t_short{1};
% 
%     %% Get pairwise correlations
     [corr_results, comb] = generate_pairwise_correlations(data_sm);
     corr_results = real(corr_results);
     
     if ~ismember(soma, 1)
         1
     end
     corr_results = corr_results(:,any(ismember(comb, soma(1)),2));
     corr_results = [ones(size(corr_results, 1),1),corr_results];


%     
%      for t = peak_t_short'
%         current_expe.ref.plot_value_tree(corr_results(t, :));
%         drawnow
%      end

    current_expe.ref.plot_value_tree(nanmean(corr_results(peak_t_short, :)));
    
     NaN_mask = isnan(corr_results(peak_t_short, :));
     temp = corr_results(peak_t_short, :); % remove all nans
     temp = double(temp > 0.5);
     temp(NaN_mask) = NaN;
     temp = nanmean(temp);
    
     sucess_rate{expe_idx} = temp;
     current_expe.ref.plot_value_tree(temp);
     mean_corr{expe_idx} = nanmean(corr_results(peak_t_short, :));
     t_somatic{expe_idx} = peak_t_short;
     all_corr{expe_idx} = corr_results(peak_t_short, :);


end


all_corr(cellfun(@isempty, all_corr)) = {NaN};
t_somatic(cellfun(@isempty, t_somatic)) = {NaN};
mean_corr(cellfun(@isempty, mean_corr)) = {NaN};
sucess_rate(cellfun(@isempty, sucess_rate)) = {NaN};


vari = cellfun(@(x) nanmean(double(x(:))), sucess_rate)*100;
vari = cellfun(@(x) nanmean(double(x(:))), mean_corr);
vari = cellfun(@(x) nanmean(sum(x < 0.3 & ~isnan(x))./sum(~isnan(x))), all_corr);



figure();scatter(cell2mat(depth)*-1, vari, 'filled', 'k'); xlabel('soma depth (um)');ylabel('% global'); title('Global event spread vs soma depth')

vari = cellfun(@(x) nanmean(sum(x < 0.8 & ~isnan(x))./sum(~isnan(x))), all_corr);
vari = cellfun(@(x) nanmean(sum(x < 0.8 & ~isnan(x))./sum(~isnan(x))), all_corr);
figure();cla();hold on;xlabel('Soma Depth');ylabel('% global'); whitebg('w'); hold on;set(gcf,'Color','w')
step = 300;
for w = -1600:step:0
    idx = find(cell2mat(depth)*-1 > w & cell2mat(depth)*-1 < w+step);
    if numel(idx) > 1
        bar(w+step/2, nanmean(vari(idx)), step-(step/20), 'k'); hold on;
        errorbar(w+step/2, nanmean(vari(idx)), nanstd(vari(idx))/sqrt(numel(idx)), 'k'); hold on;
    end
end

