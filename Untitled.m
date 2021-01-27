figure(1014);
ax = gca
Xs = [];
Ys = [];
for el = 1:numel(ax.Children)
    Xs = [Xs, ax.Children(el).XData];
    Ys = [Ys, ax.Children(el).YData];
end
 
Ys = sr./Ys;

Ym = [];
Ysd= [];
YN = [];
Xm = unique(Xs);
for bin = Xm
    Ym(bin) = nanmean(Ys(Xs == bin));  
    Ysd(bin) = nanstd(Ys(Xs == bin));  
    YN(bin) = numel(Ys(Xs == bin));  
end

figure();plot(Xm,Ym, 'ko-'); hold on
errorbar(Xm, Ym, Ysd ./ sqrt(YN), 'ko-') ;
title(['mean tau per group']);set(gcf,'Color','w');
xlabel('bin')
ylabel('tau decay')
xticklabels(bin_legend);xtickangle(45);
