classdef ExtendedKalmanFilterX<KalmanFilterX
% ExtendedKalmanFilterX class
%
% Summary of ExtendedKalmanFilterX:
% This is a class implementation of an Extended Kalman Filter.
%
% ExtendedKalmanFilterX Properties: (**)
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
%   + Jacobians - A structure containing the last computed jacobians for
%                 the transition, measurement and control gain matrices
%   + Model - An object handle to StateSpaceModelX object
%       - Dyn (*)  = Object handle to DynamicModelX SubClass     | (TO DO: LinearGaussDynModelX) 
%       - Obs (*)  = Object handle to ObservationModelX SubClass | (TO DO: LinearGaussObsModelX)
%       - Ctr (*)  = Object handle to ControlModelX SubClass     | (TO DO: LinearCtrModelX)
%
%   (*)  Signifies properties necessary to instantiate a class object
%   (**) xDim, yDim and uDim denote the dimentionality of the state, measurement
%        and control vectors respectively.
%
% ExtendedKalmanFilterX Methods:
%   + ExtendedKalmanFilterX  - Constructor method
%   + predict        - Performs KF prediction step
%   + update         - Performs KF update step
%   + iterate        - Performs a complete KF iteration (Predict & Update)
%   + smooth         - Performs KF smoothing on a provided set of estimates
%
% (+) denotes puplic properties/methods
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    properties
        Jacobians
    end
    
    methods
        function this = ExtendedKalmanFilterX(varargin)
        % EXTENDEDKALMANFILTER Constructor method
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
        % * ekf = ExtendedKalmanFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * ekf = ExtendedKalmanFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * ekf = ExtendedKalmanFilterX(ssm,priorStateMean,priorStateCov) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the prorStateMean 
        %   and priorStateCov variables.
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
        
        function predict(this)
        % PREDICT Perform an Extended Kalman Filter prediction step
        %   
         % Usage
        % -----
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % More details
        % ------------
        % * ExtendedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed using the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        %  See also update, smooth.
        
            % Predict state and measurement
            this.predictState();
            this.predictObs();                                                                
        end
        
        function predictState(this)
        % PREDICTSTATE Perform an Extended Kalman Filter state prediction step
        %   
         % Usage
        % -----
        % * predict(this) calculates the predicted system state and covariance.
        %
        % More details
        % ------------
        % * ExtendedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed using the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        %  See also update, smooth.
        
            % Extract model parameters
            f = @(x) this.Model.Dyn.feval(x);
            Q = this.Model.Dyn.covariance();
            if(~isempty(this.Model.Ctr))
                b   = @(x) this.Model.Ctr.beval(x);
                Qu  = this.Model.Ctr.covariance();
            else
                this.ControlInput   = 0;
                b   = @(x) 0;
                Qu  = 0;
            end
            
            % Perform prediction
            [this.PredStateMean, this.PredStateCovar, ...
             this.Jacobians.TransitionMatrix,this.Jacobians.ControlGain] = ...
                this.predictState_(this.StateMean, this.StateCovar,...
                                   f, Q, this.ControlInput, b, Qu);                                                                  
        end
        
        function predictObs(this)
        % PREDICTOBS Perform an Extended Kalman Filter mesurement prediction step
        %   
         % Usage
        % -----
        % * predict(this) calculates the predicted measurement,
        %   as well as the associated uncertainty covariances.
        %
        % More details
        % ------------
        % * ExtendedKalmanFilterX() uses the Model class property, which should be an
        %   instance/sublclass of the TrackingX.Models.StateSpaceModel class, in order
        %   to extract information regarding the underlying state-space model.
        % * State prediction is performed using the Model.Dyn property,
        %   which must be a subclass of TrackingX.Abstract.DynamicModel and
        %   provide the following interface functions:
        %   - Model.Dyn.feval(): Returns the model transition matrix
        %   - Model.Dyn.covariance(): Returns the process noise covariance
        % * Measurement prediction and innovation covariance calculation is
        %   performed using the Model.Obs class property, which should be
        %   a subclass of TrackingX.Abstract.DynamicModel and provide the
        %   following interface functions:
        %   - Model.Obs.heval(): Returns the model measurement matrix
        %   - Model.Obs.covariance(): Returns the measurement noise covariance
        %
        %  See also update, smooth.
        
            % Extract model parameters
            h = @(x) this.Model.Obs.heval(x);
            R = this.Model.Obs.covariance();
                        
            % Perform prediction
            [this.PredMeasMean, this.InnovErrCovar, ...
             this.CrossCovar, this.Jacobians.MeasurementMatrix] = ...
                this.predictObs_(this.PredStateMean, this.PredStateCovar,h, R);                                                                  
        end
        
        function update(this)
        % UPDATE Perform Extended Kalman Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
            
            if(isempty(this.PredMeasMean)||isempty(this.InnovErrCovar)||isempty(this.CrossCovar))
                [this.PredMeasMean, this.InnovErrCovar, this.CrossCovar] = ...
                    this.predictObs_(this.PredStateMean,this.PredStateCovar,this.Model.Obs.heval(),this.Model.Obs.covariance());
            end
            
            % Call SuperClass method
            update@KalmanFilterX(this);
        
        end
        
        function updatePDA(this, assocWeights)
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
        
%             ObsNum = size(this.Params.y,2);  
%             ObsDim = size(this.Params.y,1); 
%             
%             if(~ObsNum)
%                 warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
%                 this.Params.x = this.Params.x_pred;
%                 this.Params.P = this.Params.P_pred;
%                 return;
%             end
%             
%             if(~exist('assocWeights','var'))
%                 assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
%             end
%             
%             Compute Kalman gain
%             innov_err      = this.Params.y - this.Params.y_pred(:,ones(1,ObsNum)); % error (innovation) for each sample
%             this.Params.K   = this.Params.P_pred*this.ObsModel.Params.h(this.Params.k)'/this.Params.S;  
% 
%             update
%             Pc              = (eye(size(this.DynModel.Params.f(this.Params.k),1)) - this.Params.K*this.ObsModel.Params.h(this.Params.k)*this.Params.P_pred);
%             Pc              = this.Params.P_pred - this.Params.K*this.Params.S*this.Params.K';
%             tot_innov_err   = innov_err*assocWeights(2:end)';
%             Pgag            = this.Params.K*((innov_err.*assocWeights(ones(ObsDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*this.Params.K';
%             
%             this.Params.x    = this.Params.x_pred + this.Params.K*tot_innov_err;  
%             this.Params.P    = assocWeights(1)*this.Params.P_pred + (1-assocWeights(1))*Pc + Pgag;
            % Call SuperClass method
            updatePDA@KalmanFilterX(this, assocWeights);
        end
        
        function smoothedEstimates = Smooth(this, filteredEstimates, interval)
        % Smooth - Performs EKF smoothing on a provided set of estimates
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %                            
        %   (NOTE: The filtered_estimates array can be computed by running "filtered_estimates{k} = ekf.Params" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       ekf.Smooth(filtered_estimates);
        %
        %   See also ExtendedKalmanFilterX, Predict, Update, Iterate.
        
            % Call SuperClass method
            if nargin==2
                smoothedEstimates = Smooth@KalmanFilterX(this, filteredEstimates);
            else
                smoothedEstimates = Smooth@KalmanFilterX(this, filteredEstimates, interval); 
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
        % * ekf.resetStateEstimates() resets all state related properties
            
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
        
        function [xPred, PPred, yPred, S, Pxy, F, H, B] = predict_(x,P,f,Q,h,R,u,b,Qu)
        % PREDICT_ Perform the discrete-time EKF state and measurement
        % prediction steps, under the assumption of additive process noise.
        %
        % Parameters
        % ----------
        % x: column vector
        %   The (xDim x 1) state estimate at the previous time-step.
        % P: matrix 
        %   The (xDim x xDim) state covariance matrix at the previous
        %   time-step.
        % f: function handle
        %   A (non-linear) state transition function.
        % Q: matrix
        %   The (xDim x xDim) process noise covariance matrix.
        % h: function handle
        %   A (non-linear) measurement function.
        % R: matrix 
        %   The (yDim x yDim) measurement noise covariance matrix.
        % u: column vector, optional
        %   A optional (xDim x 1) control input.
        %   If omitted, no control input is used.
        % b: function handle, optional
        %   A (non-linear) control gain function.
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
           [yPred,S,Pxy,H]    = ExtendedKalmanFilterX.predictObs_(xPred,PPred,h,R);
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

            % Prediction for state vector and covariance:
            [xPred,F] = ExtendedKalmanFilterX.computeJac_(f,x);    %nonlinear update and linearization at current state
            PPred = F*P*F' + Q;                 %partial update

            % Compute Control Input (if applicable)
            [controlInputWithGain,B] = ExtendedKalmanFilterX.computeJac_(b,u); 

            % Add control input
            xPred = xPred + controlInputWithGain;
            PPred = PPred + B*Qu*B'; 
        end
        
        function [yPred, S, Pxy, H] = predictObs_(xPred,PPred,h,R)
        % PREDICTOBS_ Perform the discrete-time EKF observation prediction 
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

            % Prediction for measurement vector and covariance
            [yPred,H] = ExtendedKalmanFilterX.computeJac_(h, xPred);    %nonlinear measurement and linearization
            S = H*PPred*H' + R;

            % Cross-covariance
            Pxy = PPred*H';  
        end
        
        function [x,P,K] = update_(xPred,PPred,y,yPred,S,Pxy)
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

            % Compute the Kalman gain
            K = Pxy/(S);

            % Compute the filtered estimates
            x = xPred + K * (y - yPred);
            P = PPred - K*S*K';
        end
        
        function [x,P,K] = updatePDA_(xPred,PPred,Y,W,yPred,S,Pxy)
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

