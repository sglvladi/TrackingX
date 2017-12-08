function [Xa, Wm, Wc] = FormAugmentedSigmas(alpha,kappa,beta,x,P,Q,R)
% FORMAUGMENTEDSIGMAS Compute and return the Augmented Sigma Points matrix
% and the respective mean and covariance sigma weight vectors.
%
% INPUTS:    alpha - UKF scaling factor
%            kappa - UKF scaling factor
%            beta  - UKF scaling factor
%            x     - The (xDim x 1) state mean vector.
%            P     - The (xDim x xDim) state covariance vector.
%            Q     - The (xDim x xDim) process noise covariance matrix.
%            R     - The (yDim x yDim) measurement noise covariance matrix.
%
% OUTPUTS:   xPred - The (xDim x 1) state estimate at the current time-step.
%            XPred - The (xDim x nSigma) predicted sigma point matrix.
%            PPred - The (xDim x xDim) predicted state covariance matrix at 
%                    the current time-step.
%            P1    - The (xDim x xDim) unweighted sigma point variance
%                    matrix
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    xDim = size(x,1);  % State dims
    wDim = xDim;       % State noise dims
    yDim = size(R,1);  % Observation dims
    vDim = yDim;       % Observation noise dims

    % Parameters for the UKF
    aDim = xDim + wDim + vDim;              % Augmented state dims
    lambda = (alpha^2)*(aDim + kappa) - aDim;

    % Augment state and covariance
    xa = [x; zeros(wDim,1); zeros(vDim,1)]; % x^a_{k-1} = [x_k-1 E[w] E[v]]
    Pa = blkdiag(P,Q,R);

    % Scaling parameters and sigma-points
    [Si,flag] = chol((aDim + lambda)*Pa, 'lower');
    if flag ~= 0
        SP = nearestSPD((aDim + lambda)*Pa); % Ensure Pa_km1 is positive semi-definite matrix
        Si = chol(SP, 'lower');
    end
    
    % Form Augmented State Sigma points
    Xa = zeros(aDim,2*aDim+1);
    Xa(:,1) = xa;
    Xa(:,2:aDim+1) = xa*ones(1,aDim) + Si(:,1:aDim);
    Xa(:,aDim+2:2*aDim+1) = xa*ones(1,aDim) - Si(:,1:aDim);
    
    % Compute Sigma Weights
    Wm = [lambda/(aDim + lambda) repmat(1/(2*(aDim + lambda)),1,2*aDim)];
    Wc = Wm;
    Wc(1,1) = Wm(1,1) + (1 - alpha^2 + beta);
end