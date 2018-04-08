function [x,P,K] = UKalmanFilterX_UpdatePDA(xPred,PPred,Y,W,yPred,S,Pxy)
% UKALMANFILTERX_UPDATEPDA Perform the discrete-time Probabilistic Data 
% Association (PDA) UKF update step, under the assumption of additive process 
% noise, for multiple measurements (as a Gaussian Mixture)
% This is essentially identical to a PDA KF Update step.
%
% Parameters
% ----------
% xPred: column vector
%   The (xDim x 1) predicted state estimate.
% PPred: matrix
%   The (xDim x xDim) predicted state covariance matrix.
% Y: matrix
%   The (yDim x nY) measurement vector.
% W: row vector
%   The (1 x nY+1) measurement association/mixture weights 
%   vector. (dummy measurement assumed at index 1)
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
    
    [x,P,K] = KalmanFilterX_UpdatePDA(xPred,PPred,Y,W,yPred,S,Pxy);
end