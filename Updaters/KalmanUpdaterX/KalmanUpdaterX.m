classdef KalmanUpdaterX < UpdaterX 
% KalmanUpdaterX class
%
% Summary of KalmanUpdaterX:
% This is a class implementation of a standard Kalman Updater.
%
% KalmanUpdaterX Methods:
%   + KalmanUpdaterX  - Constructor method
%   + update - Performs KF update step 
%   + updatePDA - Performs PDAF-KF update step
%
% (+) denotes puplic properties/methods
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes

      
    methods (Static)
         
        function [x,P,K] = update(xPred,PPred,y,yPred,S,Pxy)
        % UPDATE Perform the discrete-time KF update step, under the  
        % assumption of additive process noisem for a single measurement.
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covariance matrix.
        % y: column vector
        %   The (yDim x 1) measurement vector.
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % S: matrix
        %   The (yDim x yDim) innovation covariance matrix.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covariance matrix.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covariance matrix at the current
        %   time-step.
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.
 
            % Compute the Kalman gain
            K = Pxy/(S);
 
            % Compute the filtered estimates
            x = xPred + K * (y - yPred);
            P = PPred - K*S*K';
        end
 
        function [x,P,K] = updatePDA(xPred,PPred,Y,W,yPred,S,Pxy)
        % UPDATEPDA Perform the discrete-time Probabilistic Data 
        % Association (PDA) KF update step, under the assumption of additive process 
        % noise, for multiple measurements (as a Gaussian Mixture)
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covariance matrix.
        % Y: matrix
        %   The (yDim x nY) measurement vector.
        % W: row vector
        %   The (1 x nY+1) measurement association/mixture weights 
        %   vector. (dummy measurement assumed at index 1)
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % S: matrix
        %   The (yDim x yDim) innovation covariance matrix.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covariance matrix.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covariance matrix at the current
        %   time-step.
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.
 
            % Get size of observation vector
            [yDim,nY] = size(Y);
 
            % Compute Kalman gain
            K = Pxy/S;  
 
            % Compute innovation mean and (cross) covariance
            innov_err       = Y - yPred(:,ones(1,nY));
            tot_innov_err   = innov_err*W(2:end)';
            Pc              = PPred - K*S*K';
            Pgag            = K*((innov_err.*W(ones(yDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*K';
 
            % Compute filtered estimates
            x    = xPred + K*tot_innov_err;  
            P    = W(1)*PPred + (1-W(1))*Pc + Pgag;
        end
    end
end