function [y,Pyy,Pxy,Y] = unscentedtf(f,X,Wm,Wc,V)
% UNSCENTEDTF Compute the Unscented Tranform for a set of sigma
% points.
%
% Parameters
% ----------
% f: function handle
%   A (non-linear) function of the form y = f(x).
% X: matrix
%   The (xDim, 2*xDim+1) sigma point matrix.
% Wm: row vector
%   The (1, 2*xDim+1) mean weights vector.
% Wc: row vector
%   The (1, 2*xDim+1) covariance weights vector.
% V: matrix, optional
%   A (yDim, 2*xDim+1) noise point matrix.
%   (default = 0)
%
% Returns
% -------
% y: column vector
%   The (yDim,1) transformed mean vector.
% Pyy: matrix
%   The (yDim x yDim) transformed covariance matrix
% Pxy: matrix 
%   The (xDim x yDim) transformed sigma points matrix.
% Y: matrix
%   The (yDim x 2*xDim+1) matrix of transformed sigma-points.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    if(nargin<5)
        V = 0;
    end

    % Transform points through f
    Y = f(X);%,V);

    % Calculate mean and covariance approximation
    [y, Pyy] = sigma2gauss(Y, Wm, Wc);

    % Calculate cross-covariance
    Pxy = (X(:,2:end)-X(:,1))*(Y(:,2:end))'*Wc(2);

end