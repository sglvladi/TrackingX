function [yPred, S, Pxy, H] = ExtendedKalmanFilterX_PredictObs(xPred,PPred,h,R)
% EKALMANFILTERX_PREDICTOBS Perform the discrete-time EKF observation prediction 
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
% h: function handle
%   A (non-linear) measurement function.
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
% H: matrix
%   The computed (yDim x yDim) Jacobian measurement matrix
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Prediction for measurement vector and covariance
    [yPred,H] = ExtendedKalmanFilterX_computeJac(h, xPred);    %nonlinear measurement and linearization
    S = H*PPred*H' + R;

    % Cross-covariance
    Pxy = PPred*H';  
end