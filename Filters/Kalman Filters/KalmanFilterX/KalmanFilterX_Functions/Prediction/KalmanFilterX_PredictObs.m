function [yPred, S, Pxy] = KalmanFilterX_PredictObs(xPred,PPred,H,R)
%KALMANFILTERX_PREDICTOBS Perform the discrete-time KF observation prediction 
% step, under the assumption of additive process noise.
%
%INPUTS:    xPred   The (xDim x 1) predicted state estimate at the current
%                   time-step.
%           PPred   The (xDim x xDim) predicted state covariance matrix at 
%                   the current time-step.
%           H       An (xDim x yDim) measurement matrix.
%           R       The (yDim x yDim) measurement noise covariance matrix.
%
%OUTPUTS:   yPred   The (yDim x 1) predicted measurement estimate.
%           Pxy     The (xDim x yDim) cross-covariance matrix.
%           S       The (yDim x yDim) innovation covariance matrix.
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Compute predicted measurement mean and covariance
    yPred   = H*xPred;
    Pxy     = PPred*H'; 
    S       = H*Pxy + R;
end