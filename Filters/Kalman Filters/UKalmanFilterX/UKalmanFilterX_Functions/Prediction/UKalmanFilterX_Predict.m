function [xPred, PPred, yPred, S, Pxy] = UKalmanFilterX_Predict(alpha,kappa,beta,x,P,f,Q,h,R,u,b,Qu)
% UKALMANFILTERX_PREDICT Perform the discrete-time UKF state and measurement
% prediction steps, under the assumption of additive process noise.
% 
% NOTES: 
% * Sigma points are selected based on the augmented state mean and
%   covariance:
%       xa = blckdiag(x,w,v);
%       Pa = blckdiag(P,Q,R);
% * Control inputs are NOT currently supported.
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
    
    xDim = size(x,1);  % State dims
    wDim = xDim;       % State noise dims
    yDim = size(R,1);  % Observation dims
    vDim = yDim;       % Observation noise dims
    aDim = xDim + wDim + vDim; 
    
    % Compute Augmented Sigma Points and Weights
    [Xa, Wm, Wc] = FormAugmentedSigmas(alpha,kappa,beta,x,P,Q,R);

    % Unscented Tranform
    [xPred, XPred, PPred, X1] = UnscentedTransform(f, Xa(1:xDim,:), Wm, Wc, Xa(xDim+1:xDim+wDim,:));
    %[ctr, Ctr, ctrP, P1] = UnscentedTransform(b, Xa(1:xDim,:), Wm, Wc, Xa(xDim+1:xDim+wDim,:));
    [yPred, YPred, S, Y1] = UnscentedTransform(h, XPred, Wm, Wc, Xa(xDim+wDim+1:aDim,:));

    % Cross-covariance
    Pxy = X1*diag(Wc)*Y1';
end