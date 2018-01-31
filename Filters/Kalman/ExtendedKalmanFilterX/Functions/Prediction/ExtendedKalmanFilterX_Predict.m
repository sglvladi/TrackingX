function [xPred, PPred, yPred, S, Pxy, F, H, B] = ExtendedKalmanFilterX_Predict(x,P,f,Q,h,R,u,b,Qu)
% EKALMANFILTERX_PREDICT Perform the discrete-time EKF state and measurement
% prediction steps, under the assumption of additive process noise.
%
% INPUTS:   x   - The (xDim x 1) state estimate at the previous time-step.
%           P   - The (xDim x xDim) state covariance matrix at the previous
%                 time-step.
%           f   - A (non-linear) state transition function.
%           Q   - The (xDim x xDim) process noise covariance matrix.
%           h   - A (non-linear) measurement function.
%           R   - The (yDim x yDim) measurement noise covariance matrix.
%           u   - A (xDim x 1) control input.
%                 (Optional, Default = 0)
%           b   - A (non-linear) control gain function.
%                 (Optional, Default = 1 if u provided, 0 otherwise)
%           Qu  - A (xDim x xDim) control noise covariance
%                 matrix. (Optional, Default = 0)
%
% OUTPUTS:  xPred - The (xDim x 1) predicted state estimate.
%           PPred - The (xDim x xDim) predicted state covariance matrix.
%           yPred - The (yDim x 1) predicted measurement estimate.
%           Pxy   - The (xDim x yDim) cross-covariance matrix.
%           S     - The (yDim x yDim) innovation covariance matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    switch(nargin)
        case(6) 
            u  = 0;
            b  = 0;
            Qu = 0;
        case(7)
            b  = 1;
            Qu = 0;
        case(8)
            Qu = 0;
    end
    
   [xPred,PPred,F,B]  = ExtendedKalmanFilterX_PredictState(x,P,f,Q,u,b,Qu);
   [yPred,S,Pxy,H]    = ExtendedKalmanFilterX_PredictObs(xPred,PPred,h,R);
end