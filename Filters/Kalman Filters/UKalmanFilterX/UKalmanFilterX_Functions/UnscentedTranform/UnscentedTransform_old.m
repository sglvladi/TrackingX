function [xPred,XPred,PPred,P1] = UnscentedTransform_old(f,X,Wm,Wc,Q)
% UNSCENTEDTRANSFORM Compute the Unscented Tranform of a set of sigma
% points.
%
% INPUTS:    f     - The (non-linear) state transition function.
%            X     - The (xDim x nSigma) sigma point matrix.
%            Wm    - The (1 x nSigma) mean weights vector.
%            Wc    - The (1 x nSigma) covariance weights vector.
%            Q     - The (xDim x xDim) noise covariance matrix.
%
% OUTPUTS:   xPred - The (xDim x 1) state estimate at the current time-step.
%            XPred - The (xDim x nSigma) predicted sigma point matrix.
%            PPred - The (xDim x xDim) predicted state covariance matrix at 
%                    the current time-step.
%            P1    - The (xDim x xDim) unweighted sigma point variance
%                    matrix
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.                  
    % Propagate sigma points
    XPred = f(X) + Q;

    % Transformed mean
    xPred = sum(Wm.*XPred,2); % Weighted average

    % Transformed variance and covariance
    P1 = (XPred-xPred);       % Variance
    PPred = P1*diag(Wc)*P1';   % Weighted covariance
end