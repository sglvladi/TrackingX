function [yPred, S, Pxy] = KalmanFilterX_PredictObs(xPred,PPred,H,R)
% KALMANFILTERX_PREDICTOBS Perform the discrete-time KF observation prediction 
% step, under the assumption of additive process noise.
%
% Parameters
% ----------
% xPred: column vector
%   The (xDim x 1) predicted state estimate at the current
%   time-step.
% PPred: matrix
%   The (xDim x xDim) predicted state covariance matrix at 
%   the current time-step.
% H: matrix
%   An (xDim x yDim) measurement matrix.
% R: matrix
%   The (yDim x yDim) measurement noise covariance matrix.
%
% Returns
% -------
% yPred: column vector
%   The (yDim x 1) predicted measurement estimate.
% Pxy: matrix
%   The (xDim x yDim) cross-covariance matrix.
% S: matrix
%   The (yDim x yDim) innovation covariance matrix.
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Compute predicted measurement mean and covariance
    yPred   = H*xPred;
    Pxy     = PPred*H'; 
    S       = H*PPred*H' + R;
end