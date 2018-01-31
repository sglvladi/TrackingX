classdef KalmanFilterX < FilterX 
% KalmanFilterX class
%
% Summary of KalmanFilterX:
% This is a class implementation of a standard Kalman Filter.
%
% KalmanFilterX Properties: (**)
%   - StateMean            A (xDim x 1) vector used to store the last computed/set filtered state mean  
%   - StateCovar      A (xDim x xDim) matrix used to store the last computed/set filtered state covariance
%   - PredStateMean        A (xDim x 1) vector used to store the last computed prediicted state mean  
%   - PredStateCovar  A (xDim x xDim) matrix used to store the last computed/set predicted state covariance
%   - PredMeasMean         A (yDim x 1) vector used to store the last computed predicted measurement mean
%   - InnovErrCovar   A (yDim x yDim) matrix used to store the last computed innovation error covariance
%   - CrossCovar      A (xDim x yDim) matrix used to store the last computed cross-covariance Cov(X,Y)
%   - KalmanGain           A (xDim x yDim) matrix used to store the last computed Kalman gain%   
%   - Measurement          A (yDim x 1) matrix used to store the received measurement
%   - ControlInput         A (uDim x 1) matrix used to store the last received control input
%   - Model                An object handle to StateSpaceModelX object
%       - Dyn (*)  = Object handle to DynamicModelX SubClass     | (TO DO: LinearGaussDynModelX) 
%       - Obs (*)  = Object handle to ObservationModelX SubClass | (TO DO: LinearGaussObsModelX)
%       - Ctr (*)  = Object handle to ControlModelX SubClass     | (TO DO: LinearCtrModelX)
%
%   (*)  Signifies properties necessary to instantiate a class object
%   (**) xDim, yDim and uDim denote the dimentionality of the state, measurement
%        and control vectors respectively.
%
% KalmanFilterX Methods:
%    KalmanFilterX  - Constructor method
%    predict        - Performs KF prediction step
%    update         - Performs KF update step
%    iterate        - Performs a complete KF iteration (Predict & Update)
%    smooth         - Performs KF smoothing on a provided set of estimates
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    properties
        StateMean 
        StateCovar     
        PredStateMean        
        PredStateCovar  
        PredMeasMean         
        InnovErrCovar   
        CrossCovar 
        KalmanGain        
        Measurement    
        ControlInput     
    end
    
    methods
        function this = KalmanFilterX(varargin)
        % KALMANFILTER Constructor method
        %   
        % DESCRIPTION: 
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
                    if (isfield(varargin{1},'priorStateMean'))
                        this.priorStateMean = varargin{1}.priorStateMean;
                        this.filtStateMean  = this.priorStateMean;
                    end
                    if (isfield(varargin{1},'priorStateCov'))
                        this.priorStateCov  = varargin{1}.priorStateCov;
                        this.filtStateCov   = this.priorStateCov;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addParameter('priorStateMean',NaN);
            parser.addParameter('priorStateCov',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.priorStateMean))
                this.StateMean = parser.Results.priorStateMean;
            end
            
            if(~isnan(parser.Results.priorStateCov))
                this.StateCovar  = parser.Results.priorStateCov;
            end
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the KalmanFilter with a certain set of
        % parameters. 
        %   
        % DESCRIPTION: 
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
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    if (isfield(varargin{1},'priorStateMean'))
                        this.priorStateMean = varargin{1}.priorStateMean;
                        this.filtStateMean  = this.priorStateMean;
                    end
                    if (isfield(varargin{1},'priorStateCov'))
                        this.priorStateCov  = varargin{1}.priorStateCov;
                        this.filtStateCov   = this.priorStateCov;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addParameter('priorStateMean',NaN);
            parser.addParameter('priorStateCov',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.priorStateMean))
                this.StateMean = parser.Results.priorStateMean;
            end
            
            if(~isnan(parser.Results.priorStateCov))
                this.StateCovar  = parser.Results.priorStateCov;
            end
        end
        
        function predict(this)
        % PREDICT Perform Kalman Filter prediction step
        %   
        % DESCRIPTION: 
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % MORE DETAILS:
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
        %  See also update, smooth.
            
            % Extract model parameters
            F = this.Model.Dyn.feval();
            Q = this.Model.Dyn.covariance();
            H = this.Model.Obs.heval();
            R = this.Model.Obs.covariance();
            if(~isempty(this.Model.Ctr))
                B   = this.Model.Ctr.beval();
                Qu  = this.Model.Ctr.covariance();
            else
                this.ControlInput   = 0;
                B   = 0;
                Qu  = 0;
            end
            % Perform prediction
            [this.PredStateMean, this.PredStateCovar, this.PredMeasMean, this.InnovErrCovar, this.CrossCovar] = ...
                KalmanFilterX_Predict(this.StateMean, this.StateCovar, F, Q, H, R, this.ControlInput, Qu); 
        end
        
        
        function update(this)
        % UPDATE Perform Kalman Filter update step
        %   
        % DESCRIPTION: 
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
        
            if(size(this.Measurement,2)>1)
                error('[KF] More than one measurement have been provided for update. Use KalmanFilterX.UpdateMulti() function instead!');
            elseif size(this.Measurement,2)==0
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                this.StateMean = this.PredStateMean;
                this.StateCovar = this.PredStateCovar;
                return;
            end
        
            % Perform single measurement update
            [this.StateMean, this.StateCovar, this.KalmanGain] = ...
                KalmanFilterX_Update(this.PredStateMean,this.PredStateCovar,...
                                     this.Measurement,this.PredMeasMean,this.InnovErrCovar,this.CrossCovar);
        end
        
        function UpdatePDA(this, assocWeights)
        % UpdatePDA - Performs KF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to measurements. Default = [0, ones(1,nData)/nData];
        %
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.y = y_new; % y_new is the new measurement)
        %       kf.UpdateMulti(assocWeights);
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also KalmanFilterX, Predict, Iterate, Smooth, resample.
        
            nData = size(this.Params.y,2);  
            
            if(~nData)
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                this.Params.x = this.Params.xPred;
                this.Params.P = this.Params.PPred;
                return;
            end
            
            if(~exist('assocWeights','var'))
                warning('[KF] No association weights have been supplied to update track! Applying default "assocWeights = [0, ones(1,nData)/nData];"...');
                assocWeights = [0, ones(1,nData)/nData]; % (1 x Nm+1)
            end
            
            [this.Params.x,this.Params.P,this.Params.K] = ...
                KalmanFilterX_UpdatePDA(this.Params.xPred,this.Params.PPred,this.Params.y,...
                                        assocWeights,this.Params.yPred,this.Params.S,this.Params.Pxy);
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
    end
end