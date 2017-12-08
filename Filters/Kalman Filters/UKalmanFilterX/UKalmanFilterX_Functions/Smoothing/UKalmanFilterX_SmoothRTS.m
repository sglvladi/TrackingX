function [smoothedEstimates] = UKalmanFilterX_SmoothRTS(filteredEstimates, interval)
% UKALMANFILTERX_SMOOTHRTS Perform Rauch–Tung–Striebel UKF Smoothing
%
% INPUTS:    filteredEstimates - A (1 x N) array of structures with fields:
%            .alpha |
%            .kappa |- The UKF scaling parameters for each time-step
%            .beta  |
%            .x      - The state estimate for each time-step
%            .P      - The state covariance for each time-step
%            .f      - The state transition function for each time-step
%            iterval - An optional interval (<= N) over which to perform 
%                      smoothing
%
% OUTPUTS:   smoothedEstimates - A (1 x N) array of structures with fields:
%            .alpha |
%            .kappa |- The UKF scaling parameters for each time-step
%            .beta  |
%            .x      - The state estimate for each time-step
%            .P      - The state covariance for each time-step
%            .f      - The state transition function for each time-step
%            .C      - The smoothing gain for each time-step
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    % Allocate memory
    N                          = length(filteredEstimates);
    smoothedEstimates          = cell(1,N);
    smoothedEstimates{N}       = filteredEstimates{N}; 

    % Perform Rauch–Tung–Striebel Backward Recursion
    for k = N-1:-1:1

        xDim = size(filteredEstimates{k}.x,1);  % State dims
        wDim = xDim;       % State noise dims
        yDim = size(filteredEstimates{k}.R,1);  % Observation dims
        vDim = yDim;       % Observation noise dims
        aDim = xDim + wDim + vDim;              % Augmented state dims

        % Compute Augmented Sigma Points and Weights
        [Xa, Wm, Wc] = FormAugmentedSigmas(filteredEstimates{k}.alpha,filteredEstimates{k}.kappa,filteredEstimates{k}.beta,...
                                           filteredEstimates{k}.x,filteredEstimates{k}.P,filteredEstimates{k}.Q,filteredEstimates{k}.R);
        % Unscented Tranform
        [xaPred, XaPred, PaPred] = UnscentedTransform( filteredEstimates{k}.f, Xa(1:xDim,:), Wm, Wc, Xa(xDim+1:xDim+wDim,:));

        % Compute smoothed estimates for time k
        smoothedEstimates{k}.C     = (Xa(1:xDim,:)-filteredEstimates{k}.x)*diag(Wc)*(XaPred-xaPred)'/PaPred;%(filteredEstimates{k}.P * filteredEstimates{k+1}.Fjac' / filteredEstimates{k+1}.P_pred;
        smoothedEstimates{k}.x     = filteredEstimates{k}.x + smoothedEstimates{k}.C * (smoothedEstimates{k+1}.x - filteredEstimates{k+1}.xPred);
        smoothedEstimates{k}.P     = filteredEstimates{k}.P + smoothedEstimates{k}.C * (smoothedEstimates{k+1}.P - filteredEstimates{k+1}.PPred) * smoothedEstimates{k}.C';                            
    end
end