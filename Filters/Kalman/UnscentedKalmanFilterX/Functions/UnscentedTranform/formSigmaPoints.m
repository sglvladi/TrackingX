function sigmaPoints = formSigmaPoints(x,P,lambda)
% FORMSIGMAPOINTS Compute and return the Augmented Sigma Points matrix
% and the respective mean and covariance sigma weight vectors.
%
% Parameters
% ----------
% lambda: scalar
%   UKF scaling factor
% x: column vector
%   The (xDim x 1) state mean vector.
% P: matrix
%   The (xDim x xDim) state covariance vector.
%
% Returns
% -------
% X: matrix
%   The (xDim x 2xDim+1) matrix of sigma points
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.

    xDim = size(x,1);  % State dims

    % Scaling parameters and sigma-points
    [Si,flag] = chol(P, 'lower');
    if flag ~= 0
        SP = nearestSPD(P); % Ensure Pa_km1 is positive semi-definite matrix
        Si = chol(SP, 'lower');
    end
    
    % Scale the covariance
    Si = sqrt(lambda) * Si;
    
    % Create sigma points
    P1 = [Si -Si];
    sigmaPoints = [x x(:,ones(1,2*xDim))+P1];
end