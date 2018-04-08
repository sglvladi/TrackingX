function [x,P,K] = KalmanFilterX_Update(xPred,PPred,y,yPred,S,Pxy)
% KALMANFILTERX_UPDATE Perform the discrete-time KF update step, under the  
% assumption of additive process noisem for a single measurement.
%
% Parameters
% ----------
% xPred: column vector
%   The (xDim x 1) predicted state estimate.
% PPred: matrix
%   The (xDim x xDim) predicted state covariance matrix.
% y: column vector
%   The (yDim x 1) measurement vector.
% yPred: column vector
%   The (yDim x 1) predicted measurement estimate.
% S: matrix
%   The (yDim x yDim) innovation covariance matrix.
% Pxy: matrix
%   The (xDim x yDim) cross-covariance matrix.
%
% Returns
% -------
% x: column vector
%   The (xDim x 1) state estimate at the current time-step.
% P: matrix
%   The (xDim x xDim) state covariance matrix at the current
%   time-step.
% K: matrix
%   The (xDim x yDim) Kalman gain matrix at the current
%   time-step.
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Compute the Kalman gain
    K = Pxy/(S);

    % Compute the filtered estimates
    x = xPred + K * (y - yPred);
    P = PPred - K*S*K';
end