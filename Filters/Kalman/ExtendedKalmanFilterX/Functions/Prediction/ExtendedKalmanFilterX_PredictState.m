function [xPred, PPred, F, B] = ExtendedKalmanFilterX_PredictState(x,P,f,Q,u,b,Qu)
% EKALMANFILTERX_PREDICTSTATE Perform the discrete-time KF state prediction 
% step, under the assumption of additive process noise.
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
% u: column vector, optional
%   An optional (xDim x 1) control input.
%   If omitted, no control input is used.
% b: function handle, optional
%   A (non-linear) control gain function.
%   (Optional, Default = 1 if u provided, 0 otherwise)
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
% F: matrix
%   The computed Jacobian transition matrix
% H: matrix
%   The computed (yDim x yDim) Jacobian measurement matrix
% B: matrix, optional
%   The computed Jacobian control gain matrix
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    switch(nargin)
        case(4) 
            u  = 0;
            b  = 0;
            Qu = 0;
        case(5)
            b  = 1;
            Qu = 0;
        case(6)
            Qu = 0;
    end
    
    % Prediction for state vector and covariance:
    [xPred,F] = ExtendedKalmanFilterX_computeJac(f,x);    %nonlinear update and linearization at current state
    PPred = F*P*F' + Q;                 %partial update

    % Compute Control Input (if applicable)
    [controlInputWithGain,B] = ExtendedKalmanFilterX_computeJac(b,u); 
    
    % Add control input
    xPred = xPred + controlInputWithGain;
    PPred = PPred + B*Qu*B'; 
end