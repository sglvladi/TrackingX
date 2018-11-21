classdef UnscentedKalmanFilterX < KalmanFilterX
% UnscentedKalmanFilterX class
%
% Summary of UnscentedKalmanFilterX:
% This is a class implementation of an Unscented Kalman Filter.
%
% UnscentedKalmanFilterX Properties: (*)
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (yDim x 1) matrix used to store the received measurement
%   + ControlInput - A (uDim x 1) matrix used to store the last received control input
%   + KalmanGain - A (xDim x yDim) matrix representing the last computed Kalman Gain 
%   + Alpha  ||
%   + Kappa  || UKF scaling parameters, as described in [1]
%   + Beta   || 
%   + Model - An object handle to StateSpaceModelX object
%       + feval (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass  
%
%   (*) xDim, yDim and uDim denote the dimentionality of the state, measurement
%       and control vectors respectively.
%
% UnscentedKalmanFilterX Methods:
%   + UnscentedKalmanFilterX  - Constructor method
%   + predict        - Performs UKF prediction step
%   + update         - Performs UKF update step
%
% (+) denotes puplic properties/methods
%
% [1] E. A. Wan and R. Van Der Merwe, "The unscented Kalman filter for nonlinear estimation," 
%     Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and 
%     Control Symposium (Cat. No.00EX373), Lake Louise, Alta., 2000, pp. 153-158.
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
  
    properties
        Alpha = 0.5
        Kappa = 0
        Beta  = 2
    end
    
    methods
        function this = UnscentedKalmanFilterX(varargin)
        % KalmanFilterX Constructor method
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a GaussianStateX instance, then it will be converted in
        %   one using the extracted mean and covariance.
        % Alpha, Beta, Kappa: scalar(s), optional
        %   The UKF scaling parameters
        %   (Defaults: Alpha = 0.5, Beta = 2, Kappa = 0)
        %
        % Usage
        % -----
        % * kf = UnscentedKalmanFilterX(___,Name,Value) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth. 
                 
            % Call SuperClass method
            this@KalmanFilterX(varargin{:});
            
            tmpIndex = 0;
            for i = 1:nargin
                if(~ischar(varargin{i}))
                    tmpIndex = tmpIndex + 1;
                end
            end
            
            if(tmpIndex<nargin)
                % Otherwise, fall back to input parser
                parser = inputParser;
                parser.KeepUnmatched = true;
                parser.addParameter('Alpha',NaN);
                parser.addParameter('Kappa',NaN);
                parser.addParameter('Beta',NaN);
                parser.parse(varargin{:});

                if(~isnan(parser.Results.Alpha))
                    this.Alpha = parser.Results.Alpha;
                end

                if(~isnan(parser.Results.Kappa))
                    this.Kappa = parser.Results.Kappa;
                end

                if(~isnan(parser.Results.Beta))
                    this.Beta = parser.Results.Beta;
                end
            end
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
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addParameter('Alpha',NaN);
            parser.addParameter('Kappa',NaN);
            parser.addParameter('Beta',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.Alpha))
                this.Alpha = parser.Results.Alpha;
            end
            
            if(~isnan(parser.Results.Kappa))
                this.Kappa = parser.Results.Kappa;
            end
            
            if(~isnan(parser.Results.Beta))
                this.Beta = parser.Results.Beta;
            end
        end
        
        function [statePrediction, measurementPrediction] = predict(this, varargin)
        % Predict Perform an Unscented Kalman Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % More details
        % ------------
        % * UnscentedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Transition property,
        %   which must be a subclass of TrackingX.Abstract.TransitionModel and
        %   provide the following interface functions:
        %   - Model.Transition.feval(): Returns the model transition matrix
        %   - Model.Transition.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed using the Model.Measurement class property, which should be
        %   a subclass of TrackingX.Abstract.TransitionModel and provide the
        %   following interface functions:
        %   - Model.Measurement.feval(): Returns the model measurement matrix
        %   - Model.Measurement.covariance(): Returns the measurement noise covariance
        %
        %  See also update, smooth.
        
            % Predict state and measurement
            statePrediction = this.predictState(varargin);
            measurementPrediction = this.predictMeasurement(varargin);
        end
        
        function statePrediction = predictState(this,varargin)
        % predictState Perform an Unscented Kalman Filter prediction step
        %   
        % Usage
        % -----
        % * predictState(this) calculates the predicted system state and covariance.
        %
        %  See also update, smooth.
        
            % Extract model parameters
            f = @(x) this.Model.Transition.feval(x);
            Q = this.Model.Transition.covar();
            if(~isempty(this.Model.Control))
                b   = @(x) this.Model.Control.feval(x);
                Qu  = this.Model.Ctr.covar();
            else
                this.ControlInput   = 0;
                b   = @(x) 0;
                Qu  = 0;
            end
            
            % Perform prediction
            [statePredictionMean, statePredictionCovar] = ...
                this.predictState_(this.Alpha, this.Kappa, this.Beta,...
                              this.StatePosterior.Mean, this.StatePosterior.Covar,...
                              f, Q, this.ControlInput, b, Qu);  
                          
            statePrediction = GaussianStateX(statePredictionMean, statePredictionCovar);
            this.StatePrediction = statePrediction;
        end
        
        function measurementPrediction = predictMeasurement(this, varargin)
        % predictMeasurements Perform an Unscented Kalman Filter measurement 
        %   prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted measurement,
        %   as well as the associated uncertainty covars.
        %
        % More details
        % ------------
        % * UnscentedKalmanFilterX() uses the Model class property, which should be an
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
        
             % Extract model parameters
            h = @(x) this.Model.Measurement.feval(x);
            R = this.Model.Measurement.covar();
            
            % Perform prediction
            [measurementPredictionMean, measurementPredictionCovar, this.KalmanGain] = ...
                this.predictMeasurement_(this.Alpha, this.Kappa, this.Beta,...
                              this.StatePrediction.Mean, this.StatePrediction.Covar, h, R); 
            measurementPrediction = GaussianStateX(measurementPredictionMean, measurementPredictionCovar);
            this.MeasurementPrediction = measurementPrediction;
        end
        
        function posterior = update(this, varargin)
        % ppdate Perform Unscented Kalman Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covar.
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
            
            if(isempty(this.MeasurementPrediction.Mean)||isempty(this.MeasurementPrediction.Covar))
                [measurementPredictionMean, measurementPredictionCovar, this.KalmanGain] = ...
                    this.predictMeasurement_(this.Alpha, this.Kappa, this.Beta,...
                                     this.PredStateMean,this.PredStateCovar,...
                                     this.Model.Measurement.feval(),this.Model.Measurement.covar());
                measurementPrediction = GaussianStateX(measurementPredictionMean, measurementPredictionCovar);
                this.MeasurementPrediction = measurementPrediction;
            end
            
            % Call SuperClass method
            posterior = update@KalmanFilterX(this, varargin);
        end
        
        function posterior = updatePDA(this, assocWeights, varargin)
        % UPDATEPDA - Performs UKF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % DESCRIPTION:
        %  * updatePDA(assocWeights) Performs UKF-PDA update step for multiple 
        %    measurements based on the provided (1-by-Nm+1) association weights 
        %    matrix assocWeights.
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also KalmanFilterX, Predict, Iterate, Smooth, resample.
        
            % Call SuperClass method
            posterior = updatePDA@KalmanFilterX(this, assocWeights);
        end
        
    end
    
     methods (Static)
        
        function [xPred, PPred, yPred, S, K] = predict_(alpha,kappa,beta,x,P,f,Q,h,R,u,b,Qu)
        % PREDICT_ Perform the discrete-time UKF state and measurement
        % prediction steps, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % alpha, kappa, betta: scalar
        %   The UKF scaling parameters
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
        %   The (yDim x yDim) innovation covariance matrix.
        % F: matrix
        %   The computed Jacobian transition matrix 
        % B: matrix, optional
        %   The computed Jacobian control gain matrix
        %
        % October 2017 Lyudmil Vladimirov, University of Liverpool.

            switch(nargin)
                case(10) 
                    u  = 0;
                    b  = 0;
                    Qu = 0;
                case(11)
                    b  = 1;
                    Qu = 0;
                case(12)
                    Qu = 0;
            end

           [xPred,PPred]  = UnscentedKalmanFilterX.predictState_(alpha,kappa,beta,x,P,f,Q,u,b,Qu);
           [yPred,S,K]  = UnscentedKalmanFilterX.predictMeasurement_(alpha,kappa,beta,xPred,PPred,h,R);
        end
        
        function [xPred, PPred] = predictState_(alpha,kappa,beta,x,P,f,Q,u,b,Qu)
        % PREDICTSTATE_ Perform the discrete-time UKF state prediction 
        % step, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % alpha, kappa, betta: scalar
        %   The UKF scaling parameters
        % x: column vector
        %   The (xDim x 1) state estimate at the previous time-step.
        % P: matrix
        %   The (xDim x xDim) state covariance matrix at the previous
        %   time-step.
        % f: function handle
        %   A (non-linear) state transition function.
        % Q: matrix
        %   The (xDim x xDim) process noise covariance matrix.
        % u: column vector, optional
        %   An optional (xDim x 1) control input.
        %   If omitted, no control input is used.
        % b: function handle, optional
        %   A (non-linear) control gain function.
        %   (Optional, Default = 1 if u provided, 0 otherwise)
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
            
            numStateDims = size(x,1);
            
            % Calculate unscented transformation parameters
            [c, Wmean, Wcov, OOM] = ...
                matlabshared.tracking.internal.calcUTParameters(alpha,beta,kappa,numStateDims);

            % Form the sigma points
            X = formSigmaPoints(x, P, c);           

            % Perform Unscented Transform to get predicted State and Covariance     
            [xPred,PPred] = unscentedTransform(f,X,Wmean,Wcov,OOM);
            % Add uncertainty to our prediction due to process noise
            PPred = PPred + Q;
        end
        
        function [yPred, S, K] = predictMeasurement_(alpha,kappa,beta,xPred,PPred,h,R)
        % PREDICTOBS_ Perform the discrete-time UKF observation prediction 
        % step, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % alpha, kappa, betta: scalar
        %   The UKF scaling parameters
        % xPred: column vector
        %   The (xDim x 1) predicted state estimate at the current
        %   time-step.
        % PPred: matrix
        %   The (xDim x xDim) predicted state covariance matrix at 
        %   the current time-step.
        % h: function handle
        %   A (non-linear) measurement function.
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
        % H: matrix
        %   The computed (yDim x yDim) Jacobian measurement matrix
        %
        % October 2017 Lyudmil Vladimirov, University of Liverpool.
            
            numStateDims = size(xPred,1);
            
            % Calculate unscented transformation parameters
            [c, Wmean, Wcov, OOM] = matlabshared.tracking.internal.calcUTParameters(alpha,beta,kappa,numStateDims);
            
            % Form the sigma points
            X = formSigmaPoints(xPred, PPred, c);

            % Perform Unscented Transform to get predicted measurement mean,
            % covariance and cross-covariance
            [yPred,S,Pxy] = unscentedTransform(h,X,Wmean,Wcov,OOM);
            
            % Add uncertainty to our prediction due to measurement noise
            S = S + R; 
            
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

            % Compute the filtered estimates
            x = xPred + K * (y - yPred);
            P = PPred - K*S*K';
        end
        
        function [x,P,K] = updatePDA_(xPred,PPred,Y,W,yPred,S,K)
        % KALMANFILTERX_UPDATEPDA Perform the discrete-time Probabilistic Data 
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