
func = @nanmean;

w1 = LoadingsPM(:,1)/sum(LoadingsPM(:,1));
w2 = LoadingsPM(:,2)/sum(LoadingsPM(:,2));
w3 = LoadingsPM(:,3)/sum(LoadingsPM(:,3));
w4 = LoadingsPM(:,4)/sum(LoadingsPM(:,4));
w5 = LoadingsPM(:,5)/sum(LoadingsPM(:,5));
figure(1021);cla();
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

diff_1(abs(diff_1) < 1) = 0;
diff_2(abs(diff_2) < 1) = 0;
diff_3(abs(diff_3) < 1) = 0;
diff_4(abs(diff_4) < 1) = 0;
diff_5(abs(diff_5) < 1) = 0;

figure(1023);cla();
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

for comp = 4
    [pk, pkloc] = findpeaks(eval(['diff_',num2str(4)]), 'MinPeakProminence', 1);
    strong_patches = valid_ROIs(v == comp);
    load_several_experiments(strong_patches, folders(1), false);
end







