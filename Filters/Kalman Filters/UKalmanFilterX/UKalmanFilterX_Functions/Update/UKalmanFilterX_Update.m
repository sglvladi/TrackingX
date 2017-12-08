function [x,P,K] = UKalmanFilterX_Update(xPred,PPred,y,yPred,S,Pxy)
% UKALMANFILTERX_UPDATE Perform the discrete-time UKF update step, under the  
% assumption of additive process noisem for a single measurement.
% This is essentially identical to a KF Update step.
%
% INPUTS:   xPred - The (xDim x 1) predicted state estimate.
%           PPred - The (xDim x xDim) predicted state covariance matrix.
%           y     - The (yDim x 1) measurement vector.
%           yPred - The (yDim x 1) predicted measurement estimate.
%           S     - The (yDim x yDim) innovation covariance matrix.
%           Pxy   - The (xDim x yDim) cross-covariance matrix.
%
% OUTPUTS:  x - The (xDim x 1) state estimate at the current time-step.
%           P - The (xDim x xDim) state covariance matrix at the current
%                   time-step.
%           K - The (xDim x yDim) Kalman gain matrix at the current
%                   time-step.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    [x,P,K] = KalmanFilterX_Update(xPred,PPred,y,yPred,S,Pxy);
end