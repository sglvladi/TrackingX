function [mu, cov] = sigma2gauss(sigmaPoints, meanWeights, covarWeights, angleidx)
% SIGMA2GAUSS Calculate estimated mean and covariance from a given set of sigma points
% 
% Parameters
% ----------
% sigma_points : matrix 
%   A (nDim, 2*nDim+1) array containing the sigma point locations
% meanWeights: vector
%   The (1,2) vector of mean weights, where meanWeights(1) corresponds to
%   the weight for sigmaPoints(:,1) and meanWeights(2) to the remaining
%   2*nDim sigma points.
% covarWeights: vector
%   The (1,2) vector of covar weights, where covarWeights(1) corresponds to
%   the weight for sigmaPoints(:,1) and covarWeights(2) to the remaining
%   2*nDim sigma points.
% covarNoise : matrix, optional
%   Additive noise covariance matrix of size (nDim, nDim)
%   (default is 0)
% 
% Returns
% -------
% mu: column vector
%   The (nDim x 1) state mean vector.
% cov: matrix
%   The (nDim x xDim) state covariance vector.
%
% Credit: Paul Horridge

    if(nargin<4)
        angleidx=[];
    end
    [nDim, nSigma] = size(sigmaPoints);
    
    % Expand weights to match number of sigma points
    meanWeights = [meanWeights(1) repmat(meanWeights(2),1, nSigma-1)];
    covarWeights = [covarWeights(1) repmat(covarWeights(2),1, nSigma-1)];
    
    % Using circular mean and anglewrap to ensure stability if angles are
    % included in the state
    mu = circularmean(sigmaPoints, meanWeights, angleidx);
    dx = anglewrap(sigmaPoints - mu, angleidx);
    cov  = bsxfun(@times, covarWeights, dx) * dx';
end