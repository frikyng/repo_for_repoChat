PHATE_std = NaN(128,1);
PHATE_var = NaN(128,1);
PHATE_CV = NaN(128,1);

for ii = 3:5:128
    Y_PHATE_3D = phate(current_signal, 'ndim', ii, 't', []);
    close all;
   PHATE_std(ii) = std(Y_PHATE_3D(:));
   PHATE_var(ii) = var(Y_PHATE_3D(:));
   PHATE_CV(ii) = std(Y_PHATE_3D(:))/mean(Y_PHATE_3D(:)).*100;
end

figure();plot(PHATE_std,'x-'); hold on; plot(PHATE_var,'o-')
title('Variability metrics across PHATE dimensions')
ylabel('Variability')
xlabel({'N Dim (3 : 5 : 128'})
legend('STD','VAR')

% figure();plot(PHATE_CV,'x')
% title('Variability metrics across PHATE dimensions')
% ylabel('CV')
% xlabel({'N Dim (3 : 6 : 33)'})
% legend('CV')

cs_std = PHATE_std;
cs_var = PHATE_var;
cs_std(1:2) = [];
cs_var(1:2) = [];
cs_std = cumsum(cs_std);
cs_var = cumsum(cs_var);

% look through and put cum_sum data back on original x (that included NaNs)
% pre-allocate
full_cs_std = nan(size(PHATE_std));
full_cs_var = nan(size(PHATE_var));

el_list_std = find(~isnan(PHATE_std));
for ii = cs_std(1:numel(cs_std))
    full_cs_std(el_list_std(1:numel(el_list_std))) = ii;
end

el_list_var = find(~isnan(PHATE_var));
for jj = cs_var(1:numel(cs_var))
    full_cs_var(el_list_var(1:numel(el_list_var))) = jj;
end

figure();plot(cs_std,'x-'); hold on; plot(cs_var,'o-')
title('Cumulative sum of variance metrics')
ylabel('Cumul Sum')
xlabel({'N Dim (3 : 128)'})
legend('std','var')

% updated plot with NaNs
figure();plot(full_cs_std,'x-'); hold on; plot(full_cs_var,'o-')
title('Cumulative sum of variance metrics')
ylabel('Cumul Sum')
xlabel({'N Dim (3 : 5: 128)'})
legend('std','var')