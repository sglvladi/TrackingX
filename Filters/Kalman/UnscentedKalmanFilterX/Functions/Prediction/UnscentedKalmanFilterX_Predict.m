function [xPred, PPred, yPred, S, Pxy] = UnscentedKalmanFilterX_Predict(alpha,kappa,beta,x,P,f,Q,h,R,u,b,Qu)
% UKALMANFILTERX_PREDICT Perform the discrete-time UKF state and measurement
% prediction steps, under the assumption of additive process noise.
% 
% NOTES: 
% * Control inputs are NOT currently supported.
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
    
    xDim = size(x,1);  % State dims
    wDim = xDim;       % State noise dims
    yDim = size(R,1);  % Observation dims
    vDim = yDim;       % Observation noise dims
    aDim = xDim + wDim + vDim; 
    
    Ns = xDim; % # of states
            
    % Calculate unscented transformation parameters
    [c, Wmean, Wcov, OOM] = matlabshared.tracking.internal.calcUTParameters(alpha,beta,kappa,Ns);
    
    % Form the sigma points
    X = formSigmaPoints(x, P, c);           

    % Perform Unscented Transform to get predicted State and Covariance     
    [xPred,PPred] = unscentedTransform(f,X,Wmean,Wcov,OOM);
    % Add uncertainty to our prediction due to process noise
    PPred = PPred + Q;
     
    % Form the sigma points again
    X = formSigmaPoints(xPred, PPred, c);
    
    % Perform Unscented Transform to get predicted measurement mean,
    % covariance and cross-covariance
    [yPred,S,Pxy] = unscentedTransform(h,X,Wmean,Wcov,OOM);
    % Add uncertainty to our prediction due to measurement noise
    S = S + R;
    
end