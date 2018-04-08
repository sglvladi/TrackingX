function [xTrans,PTrans,PCross,XTrans] = unscentedTransform(f,X,Wm,Wc,OOM)
% UNSCENTEDTRANSFORM Compute the Unscented Tranform of a set of sigma
% points.
%
% Parameters
% ----------
% f: function handle
%   The (non-linear) state transition function.
% X: matrix
%   The (xDim x nSigma) sigma point matrix.
% Wm: row vector
%   The (1 x nSigma) mean weights vector.
% Wc: row vector
%   The (1 x nSigma) covariance weights vector.
% OOM: scalar
%   Order of magnitude associated with Wm and Wc
%
% Returns
% -------
% xTrans: column vector
%   The (xTDim x 1) transformed mean vector.
% PTrans: matrix
%   The (xTim x xDim) transformed state covariance matrix
% XTrans: matrix 
%   The (xTim x nSigma) transformed sigma points matrix.
% PCross: matrix
%   The (xDim x xTDim) unweighted cross-variance matrix
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    % Transform the sigma points
    XTrans = f(X);
    
    xDim = size(X,1);
    Wm = [Wm(1) Wm(2,ones(1,2*xDim))];
    Wc = [Wc(1) Wc(2,ones(1,2*xDim))];
    
    % Calculate transformed mean and covariance
    xTrans = sum(Wm.*XTrans,2);
    % Rescale to the correct order of magnitude
    xTrans = xTrans * OOM;
    
    % Calculate transformed covariance
    PCross = XTrans-xTrans;
    PTrans = Wc(1)*(PCross(:,1)*PCross(:,1)') + Wc(2)*(PCross(:,2:end)*PCross(:,2:end)');
    % Rescale to the correct order of magnitude
    PTrans = OOM * PTrans; 

    if nargout>=3
        
        % Calculate the cross covariance
        X = X(:,2:end) - X(:,1);
        PCross = X * XTrans(:,2:end)';
        
        % Rescale to the correct order of magnitude
        PCross = PCross * (Wc(2) * OOM);
    end
end