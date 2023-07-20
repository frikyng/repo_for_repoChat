% Get size of input array and calculate window size
[numSamples, numChannels] = size(inputArray);
windowSize = ceil(data.params.points_per_s) * 0.5;

% Calculate number of covariance matrices to be calculated
numCovarianceMatrices = numSamples - windowSize + 1;

% Initialize array to store covariance matrix results
covarianceMatrixResults = zeros(1, numCovarianceMatrices);

% Calculate covariance matrix for each window of the input array
for i = 1:numCovarianceMatrices
    % Extract window of input array
    window = inputArray(i:(i+windowSize-1), :);
    
    % Calculate covariance matrix for window
    covarianceMatrix = cov(window, 1);
    
    % Store median value of covariance matrix in results array
    covarianceMatrixResults(i) = median(covarianceMatrix(:));
end

% Generate plot of input array, covariance matrix results, and ratio of
% covariance matrix results to mean of input array
figure(15); clf();
subplot(3, 1, 1); plot(inputArray); title('DF/F0');
subplot(3, 1, 2); plot(covarianceMatrixResults, 'Color', [0.8, 0.8, 0.8]);
subplot(3, 1, 3); plot(covarianceMatrixResults ./ abs(nanmean(inputArray(1:(end-windowSize+1), :), 2)'));
ylim([-1, 20]);
linkaxes([subplot(3, 1, 1), subplot(3, 1, 2), subplot(3, 1, 3)], 'x');

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