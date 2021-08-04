function [corr_results, comb] = generate_pairwise_correlations(data_in, window_size)
    if nargin < 2 || isempty(window_size)
        window_size = 10;
    end
    comb                        = nchoosek(1:size(data_in,2),2);
    corr_results                = cell(1,size(comb,1));
    
    parfor pair = 1:size(comb,1)
    %for pair = 1:size(comb,1)
        corr_results{pair}      = movcorr(data_in(:,comb(pair, 1)),data_in(:,comb(pair, 2)),[window_size, 0]);
    end
    corr_results                = real(cell2mat(corr_results)); % QQ NOT SURE WHY WE HAVE COMPLEX NUMBERS SOMETIMES
end

% function v = movcov(x,y,w)
%     w = w(1)
%     data = [x,y];
%     [M,N] = size(data);
%     numberCovarianceMatrices = M - w + 1;
%     covarianceMatrices = zeros(N,N,numberCovarianceMatrices);
%     for nc = 1:numberCovarianceMatrices
%         covarianceMatrices(:,:,nc) = cov(data(nc:(nc+w-1),:));
%     end
%     v = [nanmean(squeeze(covarianceMatrices(1,:,:)))';NaN(w-1,1)];
% end