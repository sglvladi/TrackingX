function [xPred, PPred, yPred, S, Pxy] = KalmanFilterX_Predict(x,P,F,Q,H,R,u,B,O)
% KALMANFILTERX_PREDICT Perform the discrete-time KF state and measurement
% prediction steps, under the assumption of additive process noise.
%
% Parameters
% ----------
% x: column vector
%   The (xDim x 1) state estimate at the previous time-step.
% P: matrix 
%   The (xDim x xDim) state covariance matrix at the previous
%   time-step.
% F: matrix
%   An (xDim x xDim) state transition matrix.
% Q: matrix
%   The (xDim x xDim) process noise covariance matrix.
% H: matrix
%   A (xDim x yDim) measurement matrix.
% R: matrix 
%   The (yDim x yDim) measurement noise covariance matrix.
% u: column vector, optional
%   A optional (xDim x 1) control input.
%   If omitted, no control input is used.
% B: matrix, optional
%   An optional (xDim x xDim) control gain matrix.
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
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    switch(nargin)
        case(6) 
            u = 0;
            B = 0;
            O = 0;
        case(7)
            B = 1;
            O = 0;
        case(8)
            O = 0;
    end
    
   [xPred, PPred] = KalmanFilterX_PredictState(x,P,F,Q,u,B,O);
   [yPred, S, Pxy] = KalmanFilterX_PredictObs(xPred,PPred,H,R);
end