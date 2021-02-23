
func = @nanmean;

LoadingsPM = obj.dimensionality.LoadingsPM;
valid_ROIs = find(obj.dimensionality.valid_trace_idx);
individual_traces = obj.rescaled_traces(:,valid_ROIs);
diff_thr = 1


w1 = LoadingsPM(:,1)/sum(LoadingsPM(:,1));
w2 = LoadingsPM(:,2)/sum(LoadingsPM(:,2));
w3 = LoadingsPM(:,3)/sum(LoadingsPM(:,3));
w4 = LoadingsPM(:,4)/sum(LoadingsPM(:,4));
w5 = LoadingsPM(:,5)/sum(LoadingsPM(:,5));
figure(1024);cla();
hold on; plot(func(individual_traces'.* w1, 1));
hold on; plot(func(individual_traces'.* w2, 1));
hold on; plot(func(individual_traces'.* w3, 1));
hold on; plot(func(individual_traces'.* w4, 1));
hold on; plot(func(individual_traces'.* w5, 1));
legend();set(gcf,'Color','w');
title('Weighted signal average per component')

w_tot = func(LoadingsPM, 2)/sum(nanmean(LoadingsPM, 2));
av = func(individual_traces'.* w_tot, 1);
hold on; plot(av, 'k', 'LineWidth', 2)

diff_1 = func(individual_traces'.* w1, 1)-av;
diff_2 = func(individual_traces'.* w2, 1)-av;
diff_3 = func(individual_traces'.* w3, 1)-av;
diff_4 = func(individual_traces'.* w4, 1)-av;
diff_5 = func(individual_traces'.* w5, 1)-av;

diff_1(abs(diff_1) < diff_thr) = 0;
diff_2(abs(diff_2) < diff_thr) = 0;
diff_3(abs(diff_3) < diff_thr) = 0;
diff_4(abs(diff_4) < diff_thr) = 0;
diff_5(abs(diff_5) < diff_thr) = 0;

figure(1025);cla();
hold on; plot(diff_1);
hold on; plot(diff_2);
hold on; plot(diff_3);
hold on; plot(diff_4);
hold on; plot(diff_5);
legend();set(gcf,'Color','w');
title('Weighted signal average per component')


[~, loc] = max(LoadingsPM(:,1:5)');
v = NaN(1, numel(valid_ROIs));
for ROI = 1:numel(valid_ROIs)
    roi = valid_ROIs(ROI);
    v(roi) = loc(ROI);
end
v = v(~isnan(v));

for comp = 5
    obj.plot_dim_tree(comp)
    [pk, pkloc] = findpeaks(eval(['diff_',num2str(comp)]), 'MinPeakProminence', diff_thr);

    %% find best trial here
    strong_patches = [1:3, valid_ROIs(v == comp)];
    if ~isempty(strong_patches)
        load_several_experiments(strong_patches, cellfun(@(x) x.data_folder, obj.arboreal_scans, 'UniformOutput', false), false);
        sig = nanmean(reshape(params.data{1},[], size(params.data{1}, 4)));
        figure(1025); hold on; scatter(pkloc, repmat(diff_thr, 1, numel(pkloc)), 'rv')
        figure(1060); hold on; scatter(params.timescale(pkloc), sig(pkloc), 'rv')
    end
end







