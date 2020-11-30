function lle = testloglike_factorAnalysis(Xtest, dataMean, lambda, psi )
% 
% Computes the log-likelihood of the data in Xtest
% under a multivariate Gaussian distribution with mean dataMean and 
% covariance estimated by factor analysis   
%   
% INPUTS: 
%
% Xtest     - data matrix (T x n)
% dataMean  - mean (1 x n)
% lambda    - factor loadings matrix (n x k)
% psi       - noise variance (1 x n) or diagonal matrix (n x n)
%
% where, T= no. of timepoints, n= no. of ROI/neurons, k=no. of factors
%
% OUTPUTS:
%
% lle - log likelihood
% 
% HG @ SilverLab, June 2020

    if any(size(psi))==1, psi = diag(psi); end  %diagonal noise variance

    [T,n] = size(Xtest);
    Xtest = zscore(Xtest); % unit variance because matlab's default factoran returns 
    Xtest = (Xtest - repmat(dataMean, T, 1 )); % centred 
    covMat = lambda*lambda' + psi; %total covariance, here diag==1

    % prob(x| mean=0, cov=C) = 1/[(2pi)^n/2 det(C)^1/2]  + exp( -1/2 * xC^(-1)x' )
    lle = -(1/2)*( n*log(2*pi) + logdet(covMat) ... log of scaling term summed for all timepoints
    	+ sum(sum( Xtest .* ( covMat\Xtest' )' ))/T) ; ... log of exponential term

end
