classdef KalmanPredictorX < PredictorX 
% KalmanPredictorX class
%
% Summary of KalmanPredictorX:
% This is a class implementation of a standard Kalman Predictor.
%
% KalmanPredictorX Methods:
%   + KalmanPredictorX  - Constructor method
%   + predict - Performs full KF prediction step (both state and measurement)
%   + predictState - Performs KF state prediction step
%   + predictMeasurement - Preforms KF measurement prediction step
%
% (+) denotes puplic properties/methods
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    methods (Static)
        
        function [xPred, PPred, yPred, S, Pxy] = predict(x,P,F,Q,H,R,u,B,O)
        % predict Perform the discrete-time KF state and measurement
        % prediction steps, under the assumption of additive process noise.
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
        % H: matrix
        %   A (xDim x yDim) measurement matrix.
        % R: matrix 
        %   The (yDim x yDim) measurement noise covariance matrix.
        % u: column vector, optional
        %   A optional (xDim x 1) control input.
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
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covariance matrix.
        % S: matrix
        %   The (yDim x yDim) innovation covariance matrix.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            switch(nargin)
                case(6) 
                    u = 0;
                    B = 0;
                    O = 0;
                case(7)
                    B = 1;
                    O = 0;
                case(8)
                    O = 0;
            end

           [xPred, PPred] = KalmanPredictorX.predictState(x,P,F,Q,u,B,O);
           [yPred, S, Pxy] = KalmanPredictorX.predictMeasurement(xPred,PPred,H,R);
        end
        
        function [xPred, PPred] = predictState(x,P,F,Q,u,B,Qu)
        % predictState Perform the discrete-time KF state prediction 
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

        function [yPred, S, Pxy] = predictMeasurement(xPred,PPred,H,R)
        % predictMeasurement Perform the discrete-time KF observation prediction 
        % step, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate at the current
        %   time-step.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covariance matrix at 
        %   the current time-step.
        % H: matrix
        %   An (xDim x yDim) measurement matrix.
        % R: matrix
        %   The (yDim x yDim) measurement noise covariance matrix.
        %
        % Returns
        % -------
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covariance matrix.
        % S: matrix
        %   The (yDim x yDim) innovation covariance matrix.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Compute predicted measurement mean and covariance
            yPred   = H*xPred;
            Pxy     = PPred*H'; 
            S       = H*PPred*H' + R;
        end
    end
end