classdef ExtendedKalmanFilterX < KalmanFilterX
% ExtendedKalmanFilterX class
%
% Summary of ExtendedKalmanFilterX:
% This is a class implementation of an Extended Kalman Filter.
%
% ExtendedKalmanFilterX Properties: (**)
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (yDim x 1) matrix used to store the received measurement
%   + ControlInput - A (uDim x 1) matrix used to store the last received control input
%   + KalmanGain - A (xDim x yDim) matrix representing the last computed Kalman Gain
%   + Jacobians - A structure containing the last computed jacobians
%   + Model - An object handle to StateSpaceModelX object
%       + Transition (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass     
%
%   (*)  Signifies properties necessary to instantiate a class object
%   (**) xDim, yDim and uDim denote the dimentionality of the state, measurement
%        and control vectors respectively.
%
% ExtendedKalmanFilterX Methods:
%   + ExtendedKalmanFilterX  - Constructor method
%   + predict                - Performs KF prediction step
%   + update                 - Performs KF update step
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties
        Jacobians
    end
    
    methods
        function this = ExtendedKalmanFilterX(varargin)
        % ExtendedKalmanFilterX Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a GaussianStateX instance, then it will be converted in
        %   one using the extracted mean and covariance.
        %
        % Usage
        % -----
        % * ekf = ExtendedKalmanFilterX(___,Name,Value) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments. 
        %
        %  See also predict, update, smooth.   
        
           % Call SuperClass method
            this@KalmanFilterX(varargin{:});
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the Extended KalmanFilter with a certain 
        % set of parameters. 
        %   
        % DESCRIPTION: 
        % * initialise(ekf, ssm) initialises the ExtendedKalmanFilterX object 
        %   ekf with the provided StateSpaceModelX object ssm.
        % * initialise(ekf,ssm,priorStateMean,priorStateCov) initialises 
        %   the ExtendedKalmanFilterX object kf with the provided StateSpaceModelX 
        %   object ssm and the prior information about the state, provided  
        %   in the form of the prorStateMean and priorStateCov variables.
        % * initialise(ekf,___,Name,Value,___) initialises the ExtendedKalmanFilterX 
        %   object kf with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth.   
           
            initialise@KalmanFilterX(this,varargin{:});
        end
        
        function [statePrediction, measurementPrediction] = predict(this, varargin)
        % PREDICT Perform an Extended Kalman Filter prediction step
        %   
        % Usage
        % ------
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covars.
        %
        % More details
        % ------------
        % * ExtendedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Transition property,
        %   which must be a subclass of TrackingX.Abstract.TransitionModel and
        %   provide the following interface functions:
        %   - Model.Transition.feval(): Returns the model transition matrix
        %   - Model.Transition.covar(): Returns the process noise covar
        % * Measurement prediction and innovation covar calculation is
        %   performed using the Model.Measurement class property, which should be
        %   a subclass of TrackingX.Abstract.TransitionModel and provide the
        %   following interface functions:
        %   - Model.Measurement.feval(): Returns the model measurement matrix
        %   - Model.Measurement.covar(): Returns the measurement noise covar
        %
        %  See also update, smooth.
        
            % Predict state and measurement
            statePrediction = this.predictState(varargin{:});
            if nargin>1 && isa(varargin{1},'StateX')
               % Replace a potential prior with the generated prediction
               % before forwarding the arguments to the measurement
               % prediction. Failure to do so will result in errors!!!
               varargin{1} = statePrediction; 
            end
            measurementPrediction = this.predictMeasurement(varargin{:});                                                               
        end
        
        function statePrediction = predictState(this, varargin)
        % PREDICTSTATE Perform an Extended Kalman Filter state prediction step
        %   
         % Usage
        % -----
        % * predict(this) calculates the predicted system state and covar.
        %
        %  See also update, smooth.
            
            timestamp = [];
            timestamp_old = [];
            for i = 1:min([2,nargin-1])
                if isa(varargin{i},'StateX')
                    this.StatePosterior = varargin{i};
                    timestamp_old = this.StatePosterior.Timestamp;
                elseif isdatetime(varargin{i})
                    timestamp = varargin{i};
                end
            end
            
            if isempty(timestamp)
                dt = this.Model.Transition.TimestepDuration;
                timestamp = this.StatePosterior.Timestamp;
            else
                dt = timestamp - timestamp_old;
            end
            
            % Extract model parameters
            f = @(x) this.Model.Transition.feval(x, false, dt);
            Q = this.Model.Transition.covar(dt);
            if(~isempty(this.Model.Control))
                b   = @(x) this.Model.Control.beval(x);
                Qu  = this.Model.Control.covar();
            else
                this.ControlInput   = 0;
                b   = @(x) 0;
                Qu  = 0;
            end
            
            % Perform prediction
            [statePredictionMean, statePredictionCovar, ...
             this.Jacobians.TransitionMatrix,this.Jacobians.ControlGain] = ...
                this.predictState_(this.StatePosterior.Mean, this.StatePosterior.Covar,...
                                   f, Q, this.ControlInput, b, Qu);    
                               
            statePrediction = GaussianStateX(statePredictionMean, statePredictionCovar, timestamp);
            this.StatePrediction = statePrediction;
        end
        
        function measurementPrediction = predictMeasurement(this, varargin)
        % PREDICTOBS Perform an Extended Kalman Filter mesurement prediction step
        %   
         % Usage
        % -----
        % * predict(this) calculates the predicted measurement,
        %   as well as the associated uncertainty covars.
        %
        % More details
        % ------------
        % * ExtendedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Transition property,
        %   which must be a subclass of TrackingX.Abstract.TransitionModel and
        %   provide the following interface functions:
        %   - Model.Transition.feval(): Returns the model transition matrix
        %   - Model.Transition.covar(): Returns the process noise covar
        % * Measurement prediction and innovation covar calculation is
        %   performed using the Model.Measurement class property, which should be
        %   a subclass of TrackingX.Abstract.TransitionModel and provide the
        %   following interface functions:
        %   - Model.Measurement.feval(): Returns the model measurement matrix
        %   - Model.Measurement.covar(): Returns the measurement noise covar
        %
        %  See also update, smooth.
            
            if nargin>1
                this.StatePrediction = varargin{1};
            end
            
            % Extract model parameters
            h = @(x) this.Model.Measurement.feval(x);
            R = this.Model.Measurement.covar();
                        
            % Perform prediction
            [measurementPredictionMean, measurementPredictionCovar, ...
             this.KalmanGain, this.Jacobians.MeasurementMatrix] = ...
                this.predictMeasurement_(this.StatePrediction.Mean, this.StatePrediction.Covar,h, R);
            
            measurementPrediction = GaussianStateX(measurementPredictionMean,...
                                                   measurementPredictionCovar,...
                                                   this.StatePrediction.Timestamp);
            this.MeasurementPrediction = measurementPrediction;
        end
        
        function posterior = update(this, varargin)
        % UPDATE Perform Extended Kalman Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covar.
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
            
            % Call SuperClass method
            posterior = update@KalmanFilterX(this, varargin{:});
        
        end
        
        function updatePDA(this, assocWeights, varargin)
        % UPDATEPDA - Performs EKF-PDAF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % Usage
        % -----
        %  * updatePDA(assocWeights) Performs KF-PDA update step for multiple 
        %    measurements based on the provided (1-by-Nm+1) association weights 
        %    matrix assocWeights.
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also KalmanFilterX, Predict, Iterate, Smooth, resample.
            updatePDA@KalmanFilterX(this, assocWeights, varargin{:});
        end
    end
    
    methods (Static)
        
        function [xPred, PPred, yPred, S, Pxy, F, H, B] = predict_(x,P,f,Q,h,R,u,b,Qu)
        % PREDICT_ Perform the discrete-time EKF state and measurement
        % prediction steps, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % x: column vector
        %   The (xDim x 1) state estimate at the previous time-step.
        % P: matrix 
        %   The (xDim x xDim) state covar matrix at the previous
        %   time-step.
        % f: function handle
        %   A (non-linear) state transition function.
        % Q: matrix
        %   The (xDim x xDim) process noise covar matrix.
        % h: function handle
        %   A (non-linear) measurement function.
        % R: matrix 
        %   The (yDim x yDim) measurement noise covar matrix.
        % u: column vector, optional
        %   A optional (xDim x 1) control input.
        %   If omitted, no control input is used.
        % b: function handle, optional
        %   A (non-linear) control gain function.
        %   If omitted, B is assumed to be 1.
        % O: matrix, optional
        %   An optional (xDim x xDim) control noise covar
        %   matrix. If omitted, Q is assumed to be 0.
        %
        % Returns
        % -------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covar matrix.
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covar matrix.
        % S: matrix
        %   The (yDim x yDim) innovation covar matrix.
        % F: matrix
        %   The computed Jacobian transition matrix 
        % B: matrix, optional
        %   The computed Jacobian control gain matrix
        %
        % October 2017 Lyudmil Vladimirov, University of Liverpool.

            switch(nargin)
                case(7) 
                    u  = 0;
                    b  = 0;
                    Qu = 0;
                case(8)
                    b  = 1;
                    Qu = 0;
                case(9)
                    Qu = 0;
            end

           [xPred,PPred,F,B]  = ExtendedKalmanFilterX.predictState_(x,P,f,Q,u,b,Qu);
           [yPred,S,Pxy,H]    = ExtendedKalmanFilterX.predictMeasurement_(xPred,PPred,h,R);
        end
        
        function [xPred, PPred, F, B] = predictState_(x,P,f,Q,u,b,Qu)
        % PREDICTSTATE_ Perform the discrete-time EKF state prediction 
        % step, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % x: column vector
        %   The (xDim x 1) state estimate at the previous time-step.
        % P: matrix
        %   The (xDim x xDim) state covar matrix at the previous
        %   time-step.
        % f: function handle
        %   A (non-linear) state transition function.
        % Q: matrix
        %   The (xDim x xDim) process noise covar matrix.
        % u: column vector, optional
        %   An optional (xDim x 1) control input.
        %   If omitted, no control input is used.
        % b: function handle, optional
        %   A (non-linear) control gain function.
        %   (Optional, Default = 1 if u provided, 0 otherwise)
        % O: matrix, optional
        %   An optional (xDim x xDim) control noise covar
        %   matrix. If omitted, Q is assumed to be 0.
        %
        % Returns
        % -------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covar matrix.
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

            % Prediction for state vector and covar:
            [xPred,F] = ExtendedKalmanFilterX.computeJac_(f,x);    %nonlinear update and linearization at current state
            PPred = F*P*F' + Q;                 %partial update

            % Compute Control Input (if applicable)
            [controlInputWithGain,B] = ExtendedKalmanFilterX.computeJac_(b,u); 

            % Add control input
            xPred = xPred + controlInputWithGain;
            PPred = PPred + B*Qu*B'; 
        end
        
        function [yPred, S, K, H] = predictMeasurement_(xPred,PPred,h,R)
        % PREDICTOBS_ Perform the discrete-time EKF observation prediction 
        % step, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate at the current
        %   time-step.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covar matrix at 
        %   the current time-step.
        % h: function handle
        %   A (non-linear) measurement function.
        % R: matrix
        %   The (yDim x yDim) measurement noise covar matrix.
        %
        % Returns
        % -------
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covar matrix.
        % S: matrix
        %   The (yDim x yDim) innovation covar matrix.
        % H: matrix
        %   The computed (yDim x yDim) Jacobian measurement matrix
        %
        % October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Prediction for measurement vector and covar
            [yPred,H] = ExtendedKalmanFilterX.computeJac_(h, xPred);    %nonlinear measurement and linearization
            S = H*PPred*H' + R;

            % Cross-covar
            Pxy = PPred*H';  
            
            % Compute the Kalman gain
            K = Pxy/(S);
        end
        
        function [x,P] = update_(xPred,PPred,y,yPred,S,K)
        % KALMANFILTERX_UPDATE Perform the discrete-time KF update step, under the  
        % assumption of additive process noisem for a single measurement.
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covar matrix.
        % y: column vector
        %   The (yDim x 1) measurement vector.
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % S: matrix
        %   The (yDim x yDim) innovation covar matrix.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covar matrix.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covar matrix at the current
        %   time-step.
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Compute the filtered estimates
            x = xPred + K * (y - yPred);
            P = PPred - K*S*K';
        end
        
        function [x,P] = updatePDA_(xPred,PPred,Y,W,yPred,S,K)
        % KALMANFILTERX_UPDATEPDA Perform the discrete-time Probabilistic Data 
        % Association (PDA) KF update step, under the assumption of additive process 
        % noise, for multiple measurements (as a Gaussian Mixture)
        %
        % Parameters
        % ----------
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covar matrix.
        % Y: matrix
        %   The (yDim x nY) measurement vector.
        % W: row vector
        %   The (1 x nY+1) measurement association/mixture weights 
        %   vector. (dummy measurement assumed at index 1)
        % yPred: column vector
        %   The (yDim x 1) predicted measurement estimate.
        % S: matrix
        %   The (yDim x yDim) innovation covar matrix.
        % Pxy: matrix
        %   The (xDim x yDim) cross-covar matrix.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covar matrix at the current
        %   time-step.
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Get size of observation vector
            [x,P] = updatePDA_@KalmanFilterX(xPred,PPred,Y,W,yPred,S,K);
        end
        
        function [xT,Jac] = computeJac_(fun,x)
        % COMPUTEJAC_ Compute Jacobian through complex step 
        % differentiation
        % 
        % Parameters
        % ----------
        % x: column vector
        %   A state vector
        % fun: function handle
        %   A (non-linear) transition function
        %   Must be of the form "fun = @(x)...." 
        %
        % Returns
        % -------
        % xT: column vector
        %   The transformed vector
        % Jac: matrix
        %   The computed Jacobian

            xT   = fun(x);
            xDim    = numel(x);
            h       = xDim*eps;
            Jac = matlabshared.tracking.internal.numericJacobian(fun, {x});
            %Jac     = imag(fun(x(:,ones(1,xDim))+eye(xDim)*h*1i))./h;
        end
    end
end

