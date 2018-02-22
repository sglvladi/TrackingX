classdef ExtendedKalmanFilterX<KalmanFilterX
% ExtendedKalmanFilterX class
%
% Summary of ExtendedKalmanFilterX:
% This is a class implementation of an Extended Kalman Filter.
%
% ExtendedKalmanFilterX Properties: (**)
%   - StateMean         A (xDim x 1) vector used to store the last computed/set filtered state mean  
%   - StateCovar        A (xDim x xDim) matrix used to store the last computed/set filtered state covariance
%   - PredStateMean     A (xDim x 1) vector used to store the last computed prediicted state mean  
%   - PredStateCovar    A (xDim x xDim) matrix used to store the last computed/set predicted state covariance
%   - PredMeasMean      A (yDim x 1) vector used to store the last computed predicted measurement mean
%   - InnovErrCovar     A (yDim x yDim) matrix used to store the last computed innovation error covariance
%   - CrossCovar        A (xDim x yDim) matrix used to store the last computed cross-covariance Cov(X,Y)
%   - KalmanGain        A (xDim x yDim) matrix used to store the last computed Kalman gain%   
%   - Measurement       A (yDim x 1) matrix used to store the received measurement
%   - ControlInput      A (uDim x 1) matrix used to store the last received control input
%   - Jacobians         A structure containing the last computed jacobians for
%                       the transition, measurement and control gain matrices
%   - Model             An object handle to StateSpaceModelX object
%       - Dyn (*)  = Object handle to DynamicModelX SubClass     | (TO DO: LinearGaussDynModelX) 
%       - Obs (*)  = Object handle to ObservationModelX SubClass | (TO DO: LinearGaussObsModelX)
%       - Ctr (*)  = Object handle to ControlModelX SubClass     | (TO DO: LinearCtrModelX)
%
%   (*)  Signifies properties necessary to instantiate a class object
%   (**) xDim, yDim and uDim denote the dimentionality of the state, measurement
%        and control vectors respectively.
%
% ExtendedKalmanFilterX Methods:
%    ExtendedKalmanFilterX  - Constructor method
%    predict        - Performs KF prediction step
%    update         - Performs KF update step
%    iterate        - Performs a complete KF iteration (Predict & Update)
%    smooth         - Performs KF smoothing on a provided set of estimates
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    properties
        Jacobians
    end
    
    methods
        function this = ExtendedKalmanFilterX(varargin)
        % EXTENDEDKALMANFILTER Constructor method
        %   
        % DESCRIPTION: 
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
        % DESCRIPTION: 
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % MORE DETAILS:
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
            h = @(x) this.Model.Obs.heval(x);
            R = this.Model.Obs.covariance();
            if(~isempty(this.Model.Ctr))
                b   = @(x) this.Model.Ctr.beval(x);
                Qu  = this.Model.Ctr.covariance();
            else
                this.ControlInput   = 0;
                b   = @(x) 0;
                Qu  = 0;
            end
            
            % Perform prediction
            [this.PredStateMean, this.PredStateCovar, this.PredMeasMean,...
             this.InnovErrCovar, this.CrossCovar, this.Jacobians.TransitionMatrix,... 
             this.Jacobians.MeasurementMatrix, this.Jacobians.ControlGain] = ...
                ExtendedKalmanFilterX_Predict(this.StateMean, this.StateCovar,...
                                              f, Q, h, R, this.ControlInput, b, Qu);                                                                  
        end   
        
        function update(this)
        % UPDATE Perform Extended Kalman Filter update step
        %   
        % DESCRIPTION: 
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
        
            % Call SuperClass method
            update@KalmanFilterX(this);
        
        end
        
        function updatePDA(this, assocWeights)
        % UPDATEPDA - Performs EKF-PDAF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % DESCRIPTION:
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
    end
end

