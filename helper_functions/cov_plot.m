[M,N] = size(all_concat);
win = ceil(data.params.points_per_s)*0.5;
numberCovarianceMatrices = M - win + 1;
%cov_results = zeros(N, N, numberCovarianceMatrices);
cov_results = zeros(1, numberCovarianceMatrices);

for nc = 1:numberCovarianceMatrices
    v = cov(all_concat(nc:(nc+win-1),:),1);
    cov_results(nc) = median(v(:));
    %cov_results(:,:,nc) = cov(all_concat(nc:(nc+win-1),:));
    %figure(1);cla();plot(all_concat(nc:(nc+win-1),:));
    %figure(2);cla();imagesc(cov(all_concat(nc:(nc+win-1),:)));colorbar
    
end
figure(15);clf();
ax1 = subplot(3,1,1);plot(all_concat);title('DF/F0')
ax2 = subplot(3,1,2);plot(cov_results,'Color',[0.8,0.8,0.8]);%set(ax2, 'YScale', 'log');
ax3 = subplot(3,1,3);plot(cov_results ./ abs(nanmean(all_concat(1:(end-win+1),:),2)'));ylim([-1,20]);%set(ax3, 'YScale', 'log');
linkaxes([ax1,ax2,ax3],'x')


% win = ceil(data.params.points_per_s);
% comb = nchoosek(1:size(all_concat,2),2);
% corr_results = {};
% for pair = 1:size(comb,1)
%     corr_results{pair} = movcorr(all_concat(:,comb(pair, 1)),all_concat(:,comb(pair, 2)),[win, 0]);
% end
% 
% corr_results = cell2mat(corr_results);
% consensus = 1./var(corr_results');
% out = max(consensus(:)) / 10;
% consensus(consensus > out) = NaN;
% 
% 
% figure(15);clf();
% ax1 = subplot(3,1,1);plot(all_concat);title('DF/F0')
% ax2 = subplot(3,1,2);plot(corr_results,'Color',[0.8,0.8,0.8]); hold on; plot(nanmean(corr_results'),'r');title('pairwise temporal correlation')
% ax3 = subplot(3,1,3);plot(consensus);title('Precision (1/Var)');set(ax3, 'YScale', 'log');
% linkaxes([ax1,ax2, ax3],'x')