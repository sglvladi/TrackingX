function [smoothedEstimates] = EKalmanFilterX_SmoothRTS(filteredEstimates,interval)
% EKALMANFILTERX_SMOOTHRTS Perform Rauch–Tung–Striebel EKF Smoothing
% Essentially identical to KF RTS Smoothing.
%
% INPUTS:    filteredEstimates - A (1 x N) array of structures with fields:
%            .x      - The state estimate for each time-step
%            .P      - The state covariance for each time-step
%            .F      - The (Jacobian) state transition matrix for each 
%                      time-step
%            iterval - An optional interval (<= N) over which to
%                      perform smoothing
%
% OUTPUTS:   smoothedEstimates - A (1 x N) array of structures with fields:
%            .x - The state estimate for each time-step
%            .P - The state covariance for each time-step
%            .F - The state transition matrix for
%                 each time-step
%            .C - The smoothing gain for each
%                 time-step
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    if nargin==1
        smoothedEstimates = KalmanFilterX_SmoothRTS(filteredEstimates);
    else
        smoothedEstimates = KalmanFilterX_SmoothRTS(filteredEstimates,interval);
    end
end