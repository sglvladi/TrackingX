function [xPred, PPred] = KalmanFilterX_PredictState(x,P,F,Q,u,B,Qu)
% KALMANFILTERX_PREDICTSTATE Perform the discrete-time KF state prediction 
% step, under the assumption of additive process noise.
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
% u: column vector, optional
%   An optional (xDim x 1) control input.
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
%
%October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    switch(nargin)
        case(4) 
            u  = 0;
            B  = 0;
            Qu = 0;
        case(5)
            B  = 1;
            Qu = 0;
        case(6)
            Qu = 0;
    end
    
    % Compute predicted state mean and covariance
    xPred = F*x + B*u;
    PPred =F*P*F' + Q + B*Qu*B';
end