function [sigmaPoints, meanWeights, covarWeights] = gauss2sigma(mu,P,alpha,beta,kappa)
% GAUSS2SIGMA Compute and return the sigma points, as well as the respective 
%   mean and covariance weight vectors.
%
% Parameters
% ----------
% mu: column vector
%   The (nDim x 1) state mean vector.
% P: matrix
%   The (nDim x nDim) state covariance vector.
% alpha : double, optional
%   Spread of the sigma points. Typically `1e-3`.
%   (default is 1)
% beta : double, optional
%   Used to incorporate prior knowledge of the distribution
%   2 is optimal is the state is normally distributed.
%   (default is 2)
% kappa : double, optional
% Secondary spread scaling parameter
%   (default is calculated as `3-Ns`)
%
% Returns
% -------
% sigmaPoints: matrix
%   The (nDim, 2*nDim+1) matrix of sigma points
% meanWeights: vector
%   The (1,2) vector of mean weights, where meanWeights(1) corresponds to
%   the weight for sigmaPoints(:,1) and meanWeights(2) to the remaining
%   2*nDim sigma points.
% covarWeights: vector
%   The (1,2) vector of covar weights, where covarWeights(1) corresponds to
%   the weight for sigmaPoints(:,1) and covarWeights(2) to the remaining
%   2*nDim sigma points.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    nDim = size(mu,1);
    switch(nargin)
        case(2)
            alpha = 1e-3;
            beta = 2;
            kappa = 3-nDim;
        case(3)
            beta = 2;
            kappa = 3-nDim;
        case(4)
            kappa = 3-nDim;
    end
    
    % Calculate unscented transformation parameters
    lambda = alpha^2 * (nDim + kappa) - nDim;
    c             = alpha^2 * (nDim + kappa); 
    meanWeights   = [1 - nDim/c,  1/(2*c)];
    covarWeights  = [meanWeights(1) + (1 - alpha^2 + beta), meanWeights(2)];
            
    % Compute square root of scaled covariance
    Si = sqrtm(c*P);
    
    % Create sigma points
    sigmaPoints = [mu mu+Si mu-Si];
end