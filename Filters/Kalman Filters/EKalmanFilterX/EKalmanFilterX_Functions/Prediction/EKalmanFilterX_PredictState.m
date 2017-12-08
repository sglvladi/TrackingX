function [xPred, PPred, F, B] = EKalmanFilterX_PredictState(x,P,f,Q,u,b,Qu)
% EKALMANFILTERX_PREDICTSTATE Perform the discrete-time KF state prediction 
% step, under the assumption of additive process noise.
%
% INPUTS:   x   - The (xDim x 1) state estimate at the previous time-step.
%           P   - The (xDim x xDim) state covariance matrix at the previous
%                 time-step.
%           f   - A (non-linear) state transition function.
%           Q   - The (xDim x xDim) process noise covariance matrix.
%           u   - A (xDim x 1) control input.
%                 (Optional, Default = 0)
%           b   - A (non-linear) control gain function.
%                 (Optional, Default = 1 if u provided, 0 otherwise)
%           Qu  - A (xDim x xDim) control noise covariance
%                 matrix. (Optional, Default = 0)
%
% OUTPUTS:  xPred - The (xDim x 1) predicted state estimate.
%           PPred - The (xDim x xDim) predicted state covariance matrix.
%           F     - The computed Jacobian transition matrix 
%           B     - The computed Jacobian control gain matrix
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
    [xPred,F] = jaccsd(f,x);    %nonlinear update and linearization at current state
    PPred = F*P*F' + Q;                 %partial update

    % Compute Control Input (if applicable)
    [controlInputWithGain,B] = jaccsd(b,u); 
    
    % Add control input
    xPred = xPred + controlInputWithGain;
    PPred = PPred + B*Qu*B'; 
end