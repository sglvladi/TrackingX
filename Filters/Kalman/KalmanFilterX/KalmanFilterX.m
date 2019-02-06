classdef KalmanFilterX < FilterX 
% KalmanFilterX class
%
% Summary of KalmanFilterX:
% This is a class implementation of a standard Kalman Filter.
%
% KalmanFilterX Properties: (**)
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (yDim x 1) matrix used to store the received measurement
%   + ControlInput - A (uDim x 1) matrix used to store the last received control input
%   + KalmanGain - A (xDim x yDim) matrix representing the last computed Kalman Gain
%   + Model - An object handle to StateSpaceModelX object
%       + Transition (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass     
%
%   (*)  Signifies properties necessary to instantiate a class object
%   (**) xDim, yDim and uDim denote the dimentionality of the state, measurement
%        and control vectors respectively.
%
% KalmanFilterX Methods:
%   + KalmanFilterX  - Constructor method
%   + predict        - Performs KF prediction step
%   + update         - Performs KF update step
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties
        StatePrior
        StatePrediction
        MeasurementPrediction
        StatePosterior
        KalmanGain
        ControlInput
    end
    
    properties (Dependent)
        MeasurementLikelihoods
    end
    
    properties (Access=protected)
        MeasurementLikelihoods_ = [];
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            initialise_@FilterX(this,config);
            if (isfield(config,'StatePrior'))
                this.StatePrior = config.StatePrior;
                this.StatePosterior = this.StatePrior;
            end
        end
        % MeasurementList
        function measurementList = setMeasurementList(this, newMeasurementList)
            measurementList = newMeasurementList;
            this.MeasurementLikelihoods_ = [];
        end
        function StatePrior = setStatePrior(this,newStatePrior)
            if(isa(newStatePrior,'GaussianStateX'))
                StatePrior = newStatePrior;
            else
                StatePrior = GaussianStateX(newStatePrior.Mean, newStatePrior.Covar);
            end
        end
        function StatePrediction = setStatePrediction(this,newStatePrediction)
            if(isa(newStatePrediction,'GaussianStateX'))
                StatePrediction = newStatePrediction;
            else
                StatePrediction = GaussianStateX(newStatePrediction.Mean, newStatePrediction.Covar);
            end
        end
        function MeasurementPrediction = setMeasurementPrediction(this,newMeasurementPrediction)
            if(isa(newMeasurementPrediction,'GaussianStateX'))
                MeasurementPrediction = newMeasurementPrediction;
            else
                MeasurementPrediction = GaussianStateX(newMeasurementPrediction.Mean, newMeasurementPrediction.Covar);
            end
        end
        function MeasurementLikelihoods = getMeasurementLikelihoods(this)
            if(isempty(this.MeasurementLikelihoods_))
                this.MeasurementLikelihoods_ =  mvnpdf(this.MeasurementList.Vectors',this.MeasurementPrediction.Mean',this.MeasurementPrediction.Covar)';
            end
            MeasurementLikelihoods = this.MeasurementLikelihoods_;
        end
        function StatePosterior = setStatePosterior(this,newStatePosterior)
            if(isa(newStatePosterior,'GaussianStateX'))
                StatePosterior = newStatePosterior;
            else
                StatePosterior = GaussianStateX(newStatePosterior.Mean, newStatePosterior.Covar);
            end
        end
    end
    
    methods
        function this = KalmanFilterX(varargin)
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
        %
        % Usage
        % -----
        % * kf = KalmanFilterX(___,Name,Value) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth.   
           
            % Call SuperClass method
            %this@FilterX(varargin{:});
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function initialise(this,varargin)
        % initialise Initialise the KalmanFilter with a certain set of
        %   parameters. 
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: StateX, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a GaussianStateX instance, then it will be converted in
        %   one using the extracted mean and covariance.
        % 
        % Usage
        % -----
        % * initialise(kf,___,Name,Value) initialises the KalmanFilterX 
        %   object kf with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth.   
           
            if(nargin==0)
                error("Not enough input arguments.");
            end
            
            initialise@FilterX(this);
            
            % First check to see if a structure was received
            if(nargin==2)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function [statePrediction, measurementPrediction] = predict(this, varargin)
        % Predict Perform Kalman Filter prediction step
        % 
        % Parameters
        % ----------
        % prior: GaussianStateX, optional
        %   The prior state estimate.
        % timestamp: datetime, optional
        %   A timestamp indicating the time at which prediction is
        %   performed.
        %
        % Returns
        % -------
        % GaussianStateX
        %   The generated state prediction
        % GaussianStateX, optional
        %   The generated measurement prediction
        %
        %  See also update, smooth.
            
            % Predict state and measurement
            statePrediction = this.predictState(varargin{:});
            measurementPrediction = this.predictMeasurement(varargin{:});
        end
        
        function statePrediction = predictState(this,varargin)
        % predictState Perform Kalman Filter state prediction step
        %   
        % Usage
        % -----
        % * predictState(this) calculates the predicted system state and covariance.
        %
        % See also update, smooth.
            
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
            F = this.Model.Transition.matrix(dt);
            Q = this.Model.Transition.covar(dt);
            if(~isempty(this.Model.Control))
                B   = this.Model.Control.feval();
                Qu  = this.Model.Control.covar();
            else
                this.ControlInput   = 0;
                B   = 0;
                Qu  = 0;
            end
            
            % Perform state prediction
            [statePredictionMean, statePredictionCovar] = ...
                this.predictState_(this.StatePosterior.Mean, this.StatePosterior.Covar, F, Q, this.ControlInput, B, Qu); 
            
            statePrediction = GaussianStateX(statePredictionMean, statePredictionCovar, timestamp);
            this.StatePrediction = statePrediction;
        end
        
        function measurementPrediction = predictMeasurement(this, varargin)
        % PREDICTOBS Perform Kalman Filter measurement prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted measurement,
        %   as well as the associated uncertainty covariances.
        %
        % More details
        % ------------
        % * KalmanFilterX uses the Model class property, which should be an
        %   instance of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Transition property,
        %   which must be a subclass of TrackingX.Abstract.TransitionModel and
        %   provide the following interface functions:
        %   - Model.Transition.feval(): Returns the model transition matrix
        %   - Model.Transition.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed usinf the Model.Measurement class property, which should be
        %   a subclass of TrackingX.Abstract.TransitionModel and provide the
        %   following interface functions:
        %   - Model.Measurement.heval(): Returns the model measurement matrix
        %   - Model.Measurement.covariance(): Returns the measurement noise covariance
        %
        % See also update, smooth.
        
            if nargin>1
                this.StatePrediction = varargin{1};
            end
            
            % Extract model parameters
            H = this.Model.Measurement.feval();
            R = this.Model.Measurement.covar();
                        
            % Perform prediction
            [measurementPredictionMean, measurementPredictionCovar, this.KalmanGain] = ...
                this.predictMeasurement_(this.StatePrediction.Mean, this.StatePrediction.Covar, H, R);
            
            measurementPrediction = GaussianStateX(measurementPredictionMean,... 
                                                   measurementPredictionCovar,...
                                                   this.StatePrediction.Timestamp);
            this.MeasurementPrediction = measurementPrediction;
        end
        
        function posterior = update(this, varargin)
        % UPDATE Perform Kalman Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        % See also KalmanFilterX, predict, iterate, smooth.
            
            if nargin>1
                if isa(varargin{1},'MeasurementX')
                    this.MeasurementList = varargin{1};
                elseif isa(varargin{1}, 'StateX')
                    this.StatePrediction = varargin{1};
                    this.MeasurementList = varargin{2};
                end
            end
            if(this.MeasurementList.NumMeasurements)
                timestamp = this.MeasurementList.Timestamp;
            else
                timestamp = this.StatePrediction.Timestamp;
            end
            
            if(isempty(this.MeasurementPrediction.Mean) || isempty(this.MeasurementPrediction.Covar))
                [measurementPredictionMean, measurementPredictionCovar, this.KalmanGain] = ...
                    this.predictMeasurement_(this.StatePrediction.Mean, this.StatePrediction.Covar, H, R);
                measurementPrediction = GaussianStateX(measurementPredictionMean,...
                                                       measurementPredictionCovar,...
                                                       this.StatePrediction.Timestamp);
                this.MeasurementPrediction = measurementPrediction;
            end     
        
            % Perform single measurement update
            [posteriorMean, posteriorCovar] = this.update_(this.StatePrediction.Mean,this.StatePrediction.Covar,...
                                                           this.MeasurementList.Vectors,this.MeasurementPrediction.Mean,...
                                                           this.MeasurementPrediction.Covar,this.KalmanGain);
            
            posterior = GaussianStateX(posteriorMean, posteriorCovar, timestamp);
            this.StatePosterior = posterior;
        end
        
        function posterior = updatePDA(this, assocWeights, varargin)
        % UPDATEPDA Performs KF update step, for multiple measurements
        %           Update is performed according to the generic (J)PDAF equations [1] 
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
        
            if(this.MeasurementList.NumMeasurements)
                timestamp = this.MeasurementList.Timestamp;
            else
                timestamp = this.StatePrediction.Timestamp;
            end
            
            [posteriorMean, posteriorCovar] = ...
                this.updatePDA_(this.StatePrediction.Mean,this.StatePrediction.Covar,this.MeasurementList.Vectors,...
                                assocWeights,this.MeasurementPrediction.Mean,this.MeasurementPrediction.Covar,this.KalmanGain);
            
            posterior = GaussianStateX(posteriorMean,posteriorCovar,timestamp);
            this.StatePosterior = posterior;
        end
        
        function measurementLikelihoods = get.MeasurementLikelihoods(this)
            measurementLikelihoods = getMeasurementLikelihoods(this);
        end
        
        function statePrior = get.StatePrior(this)
            statePrior = this.StatePrior;
        end
        
        function set.StatePrior(this, newStatePrior)
            this.StatePrior = setStatePrior(this, newStatePrior);
        end
        
        function statePrediction = get.StatePrediction(this)
            statePrediction = this.StatePrediction;
        end
        
        function set.StatePrediction(this, newStatePrediction)
            this.StatePrediction = setStatePrediction(this, newStatePrediction);
        end
        
        function measurementPrediction = get.MeasurementPrediction(this)
            measurementPrediction = this.MeasurementPrediction;
        end
        
        function set.MeasurementPrediction(this, newMeasurementPrediction)
            this.MeasurementPrediction = setMeasurementPrediction(this, newMeasurementPrediction);
        end
        
        function statePosterior = get.StatePosterior(this)
            statePosterior = this.StatePosterior;
        end
        
        function set.StatePosterior(this, newStatePosterior)
            this.StatePosterior = setStatePosterior(this, newStatePosterior);
        end
       
    end
    
    methods (Static)
        
        function [xPred, PPred, yPred, S, K] = predict_(x,P,F,Q,H,R,u,B,O)
        % PREDICT_ Perform the discrete-time KF state and measurement
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

           [xPred, PPred] = KalmanFilterX.predictState_(x,P,F,Q,u,B,O);
           [yPred, S, K] = KalmanFilterX.predictMeasurement_(xPred,PPred,H,R);
        end
        
        function [xPred, PPred] = predictState_(x,P,F,Q,u,B,Qu)
        % PREDICTSTATE_ Perform the discrete-time KF state prediction 
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

        function [yPred, S, K] = predictMeasurement_(xPred,PPred,H,R)
        % PREDICTOBS_ Perform the discrete-time KF observation prediction 
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
            K       = Pxy/(S);
        end

        function [x,P] = update_(xPred,PPred,y,yPred,S,K)
        % UPDATE_ Perform the discrete-time KF update step, under the  
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
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covariance matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Compute the filtered estimates
            x = xPred + K * (y - yPred);
            P = PPred - K*S*K';
        end

        function [x,P] = updatePDA_(xPred,PPred,Y,W,yPred,S,K)
        % UPDATEPDA_ Perform the discrete-time Probabilistic Data 
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
        % K: matrix
        %   The (xDim x yDim) Kalman gain matrix at the current
        %   time-step.
        %
        % Returns
        % -------
        % x: column vector
        %   The (xDim x 1) state estimate at the current time-step.
        % P: matrix
        %   The (xDim x xDim) state covariance matrix at the current
        %   time-step.
        %
        %October 2017 Lyudmil Vladimirov, University of Liverpool.

            % Get size of observation vector
            nY = size(Y,2);
            if(nY==0)
                x = xPred;
                P = PPred;
                return;
            end
            
            innov_err = Y - yPred;
            xupd = [xPred, xPred + K*innov_err];
            Pplus = PPred - K*S*K';
            
            try
                x = xupd*W';
            catch
                asd=2;
            end
            v_x = x - xupd;
            P = W(1)*(PPred + v_x(:,1)*v_x(:,1)');
            for j = 2:nY+1
                P = P + W(j)*(Pplus + v_x(:,j)*v_x(:,j)');
            end
            
%             % Compute innovation mean and (cross) covariance
%             innov_err       = Y - yPred(:,ones(1,nY));
%             tot_innov_err   = innov_err*W(2:end)';
%             Pc              = PPred - K*S*K';
%             Pgag            = K*((innov_err.*W(ones(yDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*K';
% 
%             % Compute filtered estimates
%             x    = xPred + K*tot_innov_err;  
%             P    = W(1)*PPred + (1-W(1))*Pc + Pgag;
             P    = (P+P')/2;
        end
        
        function config = getInitConfig()
            
        end
    end
end