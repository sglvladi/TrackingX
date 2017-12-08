function [x,P,K] = KalmanFilterX_UpdatePDA(xPred,PPred,Y,W,yPred,S,Pxy)
%KALMANFILTERX_UPDATEPDA Perform the discrete-time Probabilistic Data 
% Association (PDA) KF update step, under the assumption of additive process 
% noise, for multiple measurements (as a Gaussian Mixture)
%
%INPUTS:    xPred   The (xDim x 1) predicted state estimate.
%           PPred   The (xDim x xDim) predicted state covariance matrix.
%           Y       The (yDim x nY) measurement vector.
%           W       The (1 x nY+1) measurement association/mixture weights 
%                   vector. (dummy measurement assumed at index 1)
%           yPred   The (yDim x 1) predicted measurement estimate.
%           S       The (yDim x yDim) innovation covariance matrix.
%           Pxy     The (xDim x yDim) cross-covariance matrix.
%
%OUTPUTS:   x   The (xDim x 1) state estimate at the current time-step.
%           P   The (xDim x xDim) state covariance matrix at the current
%                   time-step.
%           K   The (xDim x yDim) Kalman gain matrix at the current
%                   time-step.
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
       
    % Get size of observation vector
    [yDim,nY] = size(Y);
    
    % Compute Kalman gain
    K = Pxy/S;  

    % Compute innovation mean and (cross) covariance
    innov_err       = Y - yPred(:,ones(1,nY));
    tot_innov_err   = innov_err*W(2:end)';
    Pc              = PPred - K*S*K';
    Pgag            = K*((innov_err.*W(ones(yDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*K';

    % Compute filtered estimates
    x    = xPred + K*tot_innov_err;  
    P    = W(1)*PPred + (1-W(1))*Pc + Pgag;
end