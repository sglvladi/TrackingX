classdef KalmanFilterX < FilterX 
% KalmanFilterX class
%
% Summary of KalmanFilterX:
% This is a class implementation of a standard Kalman Filter.
%
% KalmanFilterX Properties: (**)
%   + StateMean - A (xDim x 1) vector used to store the last computed/set filtered state mean  
%   + StateCovar - A (xDim x xDim) matrix used to store the last computed/set filtered state covariance
%   + PredStateMean - A (xDim x 1) vector used to store the last computed prediicted state mean  
%   + PredStateCovar - A (xDim x xDim) matrix used to store the last computed/set predicted state covariance
%   + PredMeasMean - A (yDim x 1) vector used to store the last computed predicted measurement mean
%   + InnovErrCovar - A (yDim x yDim) matrix used to store the last computed innovation error covariance
%   + CrossCovar - A (xDim x yDim) matrix used to store the last computed cross-covariance Cov(X,Y)
%   + KalmanGain - A (xDim x yDim) matrix used to store the last computed Kalman gain%   
%   + Measurement - A (yDim x 1) matrix used to store the received measurement
%   + ControlInput - A (uDim x 1) matrix used to store the last received control input
%   + Model - An object handle to StateSpaceModelX object
%       + Dyn (*)  = Object handle to DynamicModelX SubClass      
%       + Obs (*)  = Object handle to ObservationModelX SubClass 
%       + Ctr (*)  = Object handle to ControlModelX SubClass     
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
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    properties
        StateMean % A (xDim x 1) vector used to store the last computed/set filtered state mean
        StateCovar % A (xDim x xDim) matrix used to store the last computed/set filtered state covariance    
        PredStateMean % A (xDim x 1) vector used to store the last computed prediicted state mean       
        PredStateCovar % A (xDim x xDim) matrix used to store the last computed/set predicted state covariance
        PredMeasMean % A (yDim x 1) vector used to store the last computed predicted measurement mean       
        InnovErrCovar % A (yDim x yDim) matrix used to store the last computed innovation error covariance  
        CrossCovar % A (xDim x yDim) matrix used to store the last computed cross-covariance Cov(X,Y)
        KalmanGain % A (xDim x yDim) matrix used to store the last computed Kalman gain     
        ControlInput % A (uDim x 1) matrix used to store the last received control input   
    end
    
    methods
        function this = KalmanFilterX(varargin)
        % KALMANFILTER Constructor method
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % PriorStateMean: column vector, optional
        %   A (NumStateDims x 1) column vector, representing the prior
        %   state mean, which is copied over to StateMean.
        % PriorStateCovar: matrix, optional  
        %   A (NumStateDims x NumStateDims) matrix, representing the prior
        %   state covariance, which is copied over to StateCovar.
        %
        % Usage
        % -----
        % * kf = KalmanFilterX() returns an unconfigured object handle. Note
        %   that the object will need to be configured at a later instance
        %   before any call is made to it's methods.
        % * kf = KalmanFilterX(ssm) returns an object handle, preconfigured
        %   with the provided StateSpaceModelX object handle ssm.
        % * kf = KalmanFilterX(ssm,priorStateMean,priorStateCov) returns an 
        %   object handle, preconfigured with the provided StateSpaceModel 
        %   object handle ssm and the prior information about the state,  
        %   provided in the form of the prorStateMean and priorStateCov 
        %   variables.
        % * kf = KalmanFilterX(___,Name,Value) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth.   
           
            
            % Call SuperClass method
            this@FilterX(varargin{:});
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    if (isfield(varargin{1},'PriorStateMean'))
                        this.StateMean = varargin{1}.PriorStateMean;
                    end
                    if (isfield(varargin{1},'PriorStateCovar'))
                        this.StateCovar  = varargin{1}.PriorStateCovar;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addRequired('Model');
            parser.addParameter('PriorStateMean',NaN);
            parser.addParameter('PriorStateCovar',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.PriorStateMean))
                this.StateMean = parser.Results.PriorStateMean;
            end
            
            if(~isnan(parser.Results.PriorStateCovar))
                this.StateCovar  = parser.Results.PriorStateCovar;
            end
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the KalmanFilter with a certain set of
        % parameters. 
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % PriorStateMean: column vector, optional
        %   A (NumStateDims x 1) column vector, representing the prior
        %   state mean, which is copied over to StateMean.
        % PriorStateCovar: matrix, optional  
        %   A (NumStateDims x NumStateDims) matrix, representing the prior
        %   state covariance, which is copied over to StateCovar.
        % 
        % Usage
        % -----
        % * initialise(kf, ssm) initialises the KalmanFilterX object kf
        %   with the provided StateSpaceModelX object ssm.
        % * initialise(kf,ssm,priorStateMean,priorStateCov) initialises 
        %   the KalmanFilterX object kf with the provided StateSpaceModelX 
        %   object ssm and the prior information about the state, provided  
        %   in the form of the prorStateMean and priorStateCov variables.
        % * initialise(kf,___,Name,Value,___) initialises the KalmanFilterX 
        %   object kf with the options specified by one or more Name,Value 
        %   pair arguments. 
        %
        %  See also predict, update, smooth.   
           
            if(nargin==0)
                error("Not enough input arguments.");
            end
            
            initialise@FilterX(this);
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    if (isfield(varargin{1},'Model'))
                        this.Model = varargin{1}.Model;
                    end
                    if (isfield(varargin{1},'PriorStateMean'))
                        this.StateMean = varargin{1}.priorStateMean;
                        %this.filtStateMean  = this.priorStateMean;
                    end
                    if (isfield(varargin{1},'PriorStateCovar'))
                        this.StateCovar  = varargin{1}.priorStateCovar;
                        %this.filtStateCov   = this.priorStateCovar;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addParameter('Model',NaN);
            parser.addParameter('PriorStateMean',NaN);
            parser.addParameter('PriorStateCovar',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.Model))
                this.Model = parser.Results.Model;
            end
            
            if(~isnan(parser.Results.PriorStateMean))
                this.StateMean = parser.Results.PriorStateMean;
            end
            
            if(~isnan(parser.Results.PriorStateCovar))
                this.StateCovar  = parser.Results.PriorStateCovar;
            end
        end
        
        function predict(this)
        % PREDICT Perform Kalman Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % More details
        % ------------
        % * KalmanFilterX uses the Model class property, which should be an
        %   instance of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed usinf the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        % See also update, smooth.
            
            % Predict state and measurement
            this.predictState();
            this.predictObs();
        end
        
        function predictState(this)
        % PREDICTSTATE Perform Kalman Filter state prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted system state and covariance.
        %
        % More details
        % ------------
        % * KalmanFilterX uses the Model class property, which should be an
        %   instance of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed usinf the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        % See also update, smooth.
            
            % Extract model parameters
            F = this.Model.Dyn.feval();
            Q = this.Model.Dyn.covariance();
            if(~isempty(this.Model.Ctr))
                B   = this.Model.Ctr.beval();
                Qu  = this.Model.Ctr.covariance();
            else
                this.ControlInput   = 0;
                B   = 0;
                Qu  = 0;
            end
            % Perform state prediction
            [this.PredStateMean, this.PredStateCovar] = ...
                this.predictState_(this.StateMean, this.StateCovar, F, Q, this.ControlInput, B, Qu); 
        end
        
        function predictObs(this)
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
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed usinf the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        % See also update, smooth.
            
            % Extract model parameters
            H = this.Model.Obs.heval();
            R = this.Model.Obs.covariance();
            % Perform prediction
            [this.PredMeasMean, this.InnovErrCovar, this.KalmanGain] = ...
                this.predictObs_(this.PredStateMean, this.PredStateCovar, H, R); 
        end
        
        function update(this)
        % UPDATE Perform Kalman Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        % See also KalmanFilterX, predict, iterate, smooth.
        
%             if(size(this.Measurement,2)>1)
%                 error('[KF] More than one measurement have been provided for update. Use KalmanFilterX.UpdateMulti() function instead!');
%             elseif size(this.Measurement,2)==0
%                 warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
%                 this.StateMean = this.PredStateMean;
%                 this.StateCovar = this.PredStateCovar;
%                 return;
%             end
            
            if(isempty(this.PredMeasMean)||isempty(this.InnovErrCovar)||isempty(this.CrossCovar))
                [this.PredMeasMean, this.InnovErrCovar, this.KalmanGain] = ...
                    this.predictObs_(this.PredStateMean,this.PredStateCovar,this.Model.Obs.heval(),this.Model.Obs.covariance());
            end     
        
            % Perform single measurement update
            [this.StateMean, this.StateCovar] = ...
                this.update_(this.PredStateMean,this.PredStateCovar,...
                                     this.Measurement,this.PredMeasMean,this.InnovErrCovar,this.KalmanGain);
                                 
            update@FilterX(this);
        end
        
        function updatePDA(this, assocWeights)
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
        
            NumData = size(this.Measurement,2);  
            
            if(~NumData)
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                this.StateMean = this.PredStateMean;
                this.StateCovar = this.PredStateCovar;
                return;
            end
            
            if(~exist('assocWeights','var'))
                warning('[KF] No association weights have been supplied to update track! Applying default "assocWeights = [0, ones(1,nData)/nData];"...');
                assocWeights = [0, ones(1,NumData)/NumData]; % (1 x Nm+1)
            end
            
            [this.StateMean,this.StateCovar] = ...
                this.updatePDA_(this.PredStateMean,this.PredStateCovar,this.Measurement,...
                                        assocWeights,this.PredMeasMean,this.InnovErrCovar,this.KalmanGain);
        end
        
        function smoothedEstimates = smooth(this, filteredEstimates, interval)
        % Smooth - Performs KF smoothing on a provided set of estimates
        %   
        %   Inputs:
        %       filteredEstimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %                            
        %   (NOTE: The filtered_estimates array can be computed by running "filtered_estimates{k} = kf.Params" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       kf.Smooth(filteredEstimates);
        %
        %   See also KalmanFilterX, Predict, Update, Iterate.
        
            if(nargin==2)
                smoothedEstimates = KalmanFilterX_SmoothRTS(filteredEstimates);
            else
                smoothedEstimates = KalmanFilterX_SmoothRTS(filteredEstimates,interval);
            end     
        end 
        
        function resetStateEstimates(this)
        % RESETSTATEESTIMATES Reset all the state related class properties
        %   The following properties are reset upon execution:
        %       this.StateMean
        %       this.StateCovar
        %       this.PredStateMean
        %       this.PredStateCovar
        %       this.PredMeasMean
        %       this.InnovErrCovar
        %       this.CrossCovar
        %       this.KalmanGain
        %
        % Usage
        % -----
        % * kf.resetStateEstimates() resets all state related properties
            
            this.StateMean = [];
            this.StateCovar = [];
            this.PredStateMean = [];
            this.PredStateCovar = [];
            this.PredMeasMean = [];
            this.InnovErrCovar = [];
            this.CrossCovar = [];
            this.KalmanGain = [];
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
           [yPred, S, K] = KalmanFilterX.predictObs_(xPred,PPred,H,R);
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

        function [yPred, S, K] = predictObs_(xPred,PPred,H,R)
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
            nY = size(Y,2);
            
            innov_err = Y - yPred;
            xupd = [xPred, xPred + K*innov_err];
            Pplus = PPred - K*S*K';
            
            x = xupd*W';
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
%             P    = (P+P')/2;
        end
    end
end