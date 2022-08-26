% data = a.experiments(expe_idx).rescaled_traces;
% %data = vertcat(a.experiments(expe_idx).rescaled_traces);
% t = a.experiments(expe_idx).event_fitting.peak_pos;
% mean_trace = nanmedian(data, 2);
% %t = t(mean_trace(t) > 50);
% amp = NaN(size(t,2), size(data,2));
% chance_amp = NaN(size(t,2), size(data,2));
% data_sm = smoothdata(data, 'gaussian', 10);



expe_idx      = 90%21%53;%26;%53;

% fnames = {};
% for expe_idx = numel(meta.experiments):-1:1
%     current_expe  = meta.experiments(expe_idx);
%     fnames{expe_idx} = current_expe.ref.data_folder;
% end

mean_corr_bAP       = {};
mean_corr_nobAP     = {};
mean_corr_ALL       = {};
corr_frac_bAP       = {};
corr_frac_nobAP     = {};
corr_frac_ALL       = {};
correlated_nobAP    = {};
depth               = {};
durations           = arrayfun(@(x) max(x.t), meta.experiments);
corr_level          = {};


test_idx = meta.valid_expe(end);

for expe_idx = 1;%test_idx%[68,70,78]%26:-1:1;%numel(meta.experiments):-1:1%numel(meta.experiments):-1:1
    rendering       = false;
    current_expe    = meta.experiments(expe_idx);
    %idx             = current_expe.ref.indices.ROIs_list; % will duplicate the start of branches

    %% Get signal
    data            = current_expe.rescaled_traces;
    t               = current_expe.event.fitting.peak_pos;
    median_traces   = current_expe.binned_data.median_traces;
    
    %t              = t(median_traces(t) > 50); % use to clip/kep baps
    depth{expe_idx} = current_expe.ref.soma_location(3);
    max_d{expe_idx} = current_expe.ref.max_distance;
    
    %% Smooth data
    data_sm         = smoothdata(data, 'gaussian', [20, 0]);  
    
    %% Generate randomized data
    data_rand       = data_sm;
    for row = 1:size(data,2)
        data_rand(:,row) = circshift(data_sm(:,row), randi(size(data_sm, 1),1));
    end    
    
    %% Hide unrelevant ROIs
    data_sm(:,current_expe.bad_ROI_list)    = NaN;
    data_rand(:,current_expe.bad_ROI_list)  = NaN;
    
    %% Get pairwise correlations
    [corr_results, comb]    = generate_pairwise_correlations(data_sm);   
    [corr_results_rand]     = generate_pairwise_correlations(data_rand);

    % condition == 1 --> all TP
    % condition == 2 --> bAPs
    % condition == 3 --> between bAPs    
    %% Get average pairwise correlation across tree
    [corr_results_sub, n_high_corr] = get_average_pairwise_correlation(current_expe, corr_results, data_sm ,rendering, 1);
    mean_corr_ALL{expe_idx} = corr_results_sub;
    corr_frac_ALL{expe_idx} = n_high_corr/sum(~isnan(n_high_corr));

    %% Get average pairwise correlation across tree during bAPs
    [corr_results_sub, n_high_corr] = get_average_pairwise_correlation(current_expe, corr_results, data_sm ,rendering, 2);
    mean_corr_bAP{expe_idx} = corr_results_sub;
    corr_frac_bAP{expe_idx} = n_high_corr/sum(~isnan(n_high_corr));

    %% Get average pairwise correlation across tree between bAPs
    [corr_results_sub, n_high_corr, correlated_nobAP{expe_idx}] = get_average_pairwise_correlation(current_expe, corr_results, data_sm ,rendering, 3, corr_results_rand);
    mean_corr_nobAP{expe_idx} = corr_results_sub;
    corr_frac_nobAP{expe_idx} = n_high_corr;current_expe.ref.plot_value_tree(n_high_corr/numel(correlated_nobAP{expe_idx}),'',current_expe.default_handle);caxis([0,1]);title([num2str(numel(correlated_nobAP{expe_idx})), ' local events detected'])%caxis([0,numel(correlated_nobAP{expe_idx})]);
    
%     f = gcf();
%     [~, tag] = fileparts(fileparts(current_expe.source_folder));
%     saveas(f, [pwd,'/local/',tag,'_tree.png']);
%     dcm_obj = datacursormode(f);
%     set(dcm_obj,'UpdateFcn',[], 'enable', 'off');
%     savefig(f, [pwd,'/local/',tag,'_tree.fig']);
        
        %         h = current_expe.ref.plot_value_tree(nanmean(real(corr_results_sub)));caxis([0,1])
        %         pause(1)
        %         %if ~rendering
        %             n_high_corr(~n_high_corr) = NaN;
        %             current_expe.ref.plot_value_tree(n_high_corr/sum(~isnan(n_high_corr)));caxis([0,1])
        %             %current_expe.ref.plot_value_tree(n_high_corr/max(n_high_corr));caxis([0,1])
        %         %end
        %         
        %drawnow
   

end

corr_frac_ALL(cellfun(@isempty, corr_frac_ALL)) = {NaN};
corr_frac_nobAP(cellfun(@isempty, corr_frac_nobAP)) = {NaN};
corr_frac_bAP(cellfun(@isempty, corr_frac_bAP)) = {NaN};
depth(cellfun(@isempty, depth)) = {NaN};
mean_corr_nobAP(cellfun(@isempty, mean_corr_nobAP)) = {NaN};
mean_corr_ALL(cellfun(@isempty, mean_corr_ALL)) = {NaN};
mean_corr_bAP(cellfun(@isempty, mean_corr_bAP)) = {NaN};
correlated_nobAP(cellfun(@isempty, correlated_nobAP)) = {NaN};



test2_all = cellfun(@(x) nanmean(double(x(x > 0))), corr_frac_ALL );
test2_nobap = cellfun(@(x) nanmean(double(x(x > 0))),  corr_frac_nobAP);
test2_bap = cellfun(@(x) nanmean(double(x(x > 0))), corr_frac_bAP);




figure();scatter(cell2mat(depth)*-1, test2_all, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean responsive fraction vs depth')
figure();scatter(cell2mat(depth)*-1, test2_bap, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean responsive fraction vs depth during global events')
figure();scatter(cell2mat(depth)*-1, test2_nobap, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean responsive fraction vs depth between global events')



test_bap = cellfun(@(x) nanmean(double(x(:))), mean_corr_bAP);
test_all = cellfun(@(x) nanmean(double(x(:))), mean_corr_ALL);
test_nobap  = cellfun(@(x) nanmean(double(x(:))), mean_corr_nobAP);


figure();scatter(cell2mat(depth)*-1, test_all, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean tree correlation vs depth')
figure();scatter(cell2mat(depth)*-1, test_bap, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean tree correlation vs depth during global events')
figure();scatter(cell2mat(depth)*-1, test_nobap, 'filled', 'k'); xlabel('soma depth (um)');ylabel('mean correlation'); title('mean tree correlation vs depth between global events')



% vari = durations./cellfun(@numel, correlated_nobAP);
% vari = test2_nobap;
% vari = test2_nobap ./ (cellfun(@(x) size(x, 2), mean_corr_nobAP));
% vari = (cellfun(@numel, correlated_nobAP)./durations);

figure();cla();hold on;xlabel('Cell type');ylabel('rate'); whitebg('w');
R = fliplr(-1*[0, 100, 444, 644, 1200]);
vari(isinf(vari)) = NaN;
count = 1;
for w = 1:(numel(R)-1)
    w_size = R(w+1)-R(w);
    w = R(w);
    idx = find(cell2mat(depth)*-1 > w & cell2mat(depth)*-1 < w+w_size)
    bar(count, nanmean(vari(idx)), 0.95, 'k'); hold on;
    errorbar(count, nanmean(vari(idx))*20, nanstd(vari(idx))/sqrt(numel(idx)), 'k'); hold on;
    count = count+1;
end
xticklabels({'','L5b','L5a','L2-3',''})



figure();cla();hold on;xlabel('Soma Depth');ylabel('event spatial extent'); whitebg('w'); hold on;set(gcf,'Color','w')
step = 300;
for w = -1600:step:0
    idx = find(cell2mat(depth)*-1 > w & cell2mat(depth)*-1 < w+step);
    bar(w+step/2, nanmean(vari(idx)), step-(step/20), 'k'); hold on;
    errorbar(w+step/2, nanmean(vari(idx)), nanstd(vari(idx))/sqrt(numel(idx)), 'k'); hold on;
end


vari = test2_nobap ./ (cellfun(@(x) size(x, 2), mean_corr_nobAP));

vari = test2_nobap ;

vari = (cellfun(@numel, correlated_nobAP)./durations);
figure();cla();hold on;xlabel('Soma Depth');ylabel('event rate'); whitebg('w'); hold on;set(gcf,'Color','w')
step = 200;
for w = -1600:step:0
    idx = find(cell2mat(depth)*-1 > w & cell2mat(depth)*-1 < w+step);
    if numel(idx) > 1
        bar(w+step/2, nanmean(vari(idx)), step-(step/20), 'k'); hold on;
        errorbar(w+step/2, nanmean(vari(idx)), nanstd(vari(idx))/sqrt(numel(idx)), 'k'); hold on;
    end
end







rendering = true;

key = 6;
corr_results_sub = corr_results(:, comb(:,2) == key | comb(:,1) == key);
%% Add the key itself
corr_results_sub = [corr_results_sub(:,1:(key-1)), ones(size(corr_results_sub,1),1), corr_results_sub(:,key:end)];

if rendering
    h = current_expe.ref.plot_value_tree(nanmean(real(corr_results_sub)));
    t = title('tp = 0')
    ok = ~isnan(h.XData(1,:)); % Nan are added for curved scan at the end of each branch
end

convert_to_n = true;
if convert_to_n
    caxis([0,size(data_sm,2)])
    tp_active = NaN(1, size(corr_results_sub, 1));
end

for tp = 1:size(corr_results_sub, 1)
    %v = corr_results_sub(tp,idx);
    R = tp:nanmin(tp+9, size(corr_results_sub, 1));
    v = nanmean(corr_results_sub(R,idx),1);  
    v_short = nanmean(corr_results_sub(R,:),1); 

    %% Plot correlation or number of correlated sites
    if convert_to_n % Optional
        [n_corr, tp_active(tp)] = transform_corr_to_n(data_sm, comb, corr_results, tp, idx, 0.8);
        if rendering
            h.CData(:, ok)          = repmat(n_corr, 2,1);
            drawnow
        end
    elseif rendering
        h.CData(:, ok) = repmat(v, 2,1);
        t.String = ['tp = ',num2str(tp)];
        drawnow
    end

%     if (nanmean(v) < 0.3 && sum(v_short > 0.8) > 5 && data_sm(tp,key) > 20) || ismember(tp, current_expe.event.fitting.peak_pos)
%         show_traces_at_tp(data_sm, tp, median_traces, key)
%     end
    
end

if convert_to_n
    figure();plot(median_traces/nanmax(median_traces)); hold on; plot(tp_active/max(tp_active));
end

1

% 
% timescale = current_expe.t;
% figure();plot(timescale, nanmean(corr_results, 2),'r');hold on;plot(timescale, median_traces/80)
% % figure; hold on;plot(corr_results,'Color',[0.9,0.9,0.9]);
% % hold on;plot(nanmean(corr_results, 2),'r');title('pairwise temporal correlation');ylim([-1.1,1.1]);plot([0,size(corr_results, 1)],[-1,-1],'k--');plot([0,size(corr_results, 1)],[1,1],'k--')







%h = current_expe.ref.plot_value_tree(nanmean(real(corr_results_sub)));
% t = title('tp = 0')
% ok = ~isnan(h.XData(1,:)); % Nan are added for curved scan at the end of each branch
% idx = current_expe.ref.indices.ROIs_list % will duplicate the start of branches
% values = [];
% previous = 1;
% for tp = 1:size(corr_results_sub, 1)
%     %v = corr_results_sub(tp,idx);
%     v = nanmean(corr_results_sub(tp:tp+9,idx),1);  
%     v_short = nanmean(corr_results_sub(tp:tp+9,:),1); 
%     if sum(v_short > 0.8) > 5 && nanmean(v_short([1,23,41,43,45,39,36,49])) < 0.3
%         
%        if tp - previous > 1 && ~isempty(values) % empty for t = 1
%            %% save last figure
%            h.CData(:, ok) = repmat(nanmean(values,1), 2,1);
%            t.String = ['tp = ',num2str(previous)];
%            drawnow;
%            saveas(1081,['test_',num2str(previous),'.png']);
%            
%            %% Clear values
%            values = v;
%        else
%            values = [values ; v];
%        end
%        previous = tp;
%     end
% end
% 
% 






