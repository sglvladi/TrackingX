function [xSmooth_km1, PSmooth_km1, C_km1] = KalmanFilterX_SmoothRTS_Single(xFilt_km1, PFilt_km1, xPred_k, PPred_k, xSmooth_k, PSmooth_k, F_k, B_km1, u_km1)
% KALMANFILTERX_SMOOTHRTS_SINGLE Perform a single iteration of a Rauch–Tung–Striebel 
% Kalman Filter Smoother (from time k, to time k-1)
%
% INPUTS:   xFilt_km1   The (xDim x 1) filtered state vector estimate for time k-1
%           PFilt_km1   The (xDim x xDim) filtered state covariance estimate for time k-1
%           xPred_k   The (xDim x 1) predicted state vector estimate for time k
%           PPred_k   The (xDim x xDim) predicted state vector estimate for time k
%           xSmooth_k   The (xDim x 1) smoothed state vector estimate for time k
%           PSmooth_k   The (xDim x xDim) smoothed state covariance estimate for time k
%           F_k         The (xDim x xDim) transition matrix for time k
%           B_km1       The (xDim x uDim) control input gain matrix for time k-1
%           u_km1       The (uDim x 1) control input vector for time k-1
%
% OUTPUTS:  xSmooth_km1   The (xDim x 1) smoothed state vector estimate for time k-1
%           PSmooth_km1   The (xDim x xDim) smoothed state covariance estimate for time k-1
%           C_km1         The (xDim x xDim) smoothed gain estimate for time k-1
%
% December 2017 Lyudmil Vladimirov, University of Liverpool.

    if(nargin<8)
        B_km1 = 0;
        u_km1 = 0;
    end

    % Perform Rauch–Tung–Striebel Backward Iteration
    C_km1           = PFilt_km1 * F_k' / PPred_k;
    xSmooth_km1     = xFilt_km1 + C_km1 * (xSmooth_k - xPred_k - B_km1*u_km1);
    PSmooth_km1     = PFilt_km1 + C_km1 * (PSmooth_k - PPred_k) * C_km1';                            
end