function [x,P,K] = EKalmanFilterX_UpdatePDA(xPred,PPred,Y,W,yPred,S,Pxy)
% KALMANFILTERX_UPDATEPDA Perform the discrete-time Probabilistic Data 
% Association (PDA) EKF update step, under the assumption of additive process 
% noise, for multiple measurements (as a Gaussian Mixture)
% This is essentially identical to a PDA KF Update step.
%
% INPUTS:    xPred - The (xDim x 1) predicted state estimate.
%            PPred - The (xDim x xDim) predicted state covariance matrix.
%            Y     - The (yDim x nY) measurement vector.
%            W     - The (1 x nY+1) measurement association/mixture weights 
%                    vector. (dummy measurement assumed at index 1)
%            yPred - The (yDim x 1) predicted measurement estimate.
%            S     - The (yDim x yDim) innovation covariance matrix.
%            Pxy   - The optional (xDim x yDim) cross-covariance matrix.
%
% OUTPUTS:   x - The (xDim x 1) state estimate at the current time-step.
%            P - The (xDim x xDim) state covariance matrix at the current
%                   time-step.
%            K - The (xDim x yDim) Kalman gain matrix at the current
%                   time-step.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    [x,P,K] = KalmanFilterX_UpdatePDA(xPred,PPred,Y,W,yPred,S,Pxy);
end