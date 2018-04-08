function [xPred, PPred, yPred, S, Pxy, F, H, B] = ExtendedKalmanFilterX_Predict(x,P,f,Q,h,R,u,b,Qu)
% EKALMANFILTERX_PREDICT Perform the discrete-time EKF state and measurement
% prediction steps, under the assumption of additive process noise.
%
% Parameters
% ----------
% x: column vector
%   The (xDim x 1) state estimate at the previous time-step.
% P: matrix 
%   The (xDim x xDim) state covariance matrix at the previous
%   time-step.
% f: function handle
%   A (non-linear) state transition function.
% Q: matrix
%   The (xDim x xDim) process noise covariance matrix.
% h: function handle
%   A (non-linear) measurement function.
% R: matrix 
%   The (yDim x yDim) measurement noise covariance matrix.
% u: column vector, optional
%   A optional (xDim x 1) control input.
%   If omitted, no control input is used.
% b: function handle, optional
%   A (non-linear) control gain function.
%   If omitted, B is assumed to be 1.
% O: matrix, optional
%   An optional (xDim x xDim) control noise covariance
%   matrix. If omitted, Q is assumed to be 0.
%
% Returns
% -------
% xPred: column vector
%   The (xDim x 1) predicted state estimate.
% PPred: matrix
%   The (xDim x xDim) predicted state covariance matrix.
% yPred: column vector
%   The (yDim x 1) predicted measurement estimate.
% Pxy: matrix
%   The (xDim x yDim) cross-covariance matrix.
% S: matrix
%   The (yDim x yDim) innovation covariance matrix.
% F: matrix
%   The computed Jacobian transition matrix 
% B: matrix, optional
%   The computed Jacobian control gain matrix
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