classdef KalmanFilterX < FilterX % Handle class with copy functionality
% KalmanFilterX class
%
% Summary of KalmanFilterX:
% This is a class implementation of a vanilla Kalman Filter.
%
% IMPORTANT NOTE: The Kalman Filter only supports linear-Gaussian models. Thus, the model functions and noise covariances called within the filter ONLY depend on time, but NOT the state.
%                 (i.e. Calls made within the filter are of the form this.DynModel.sys(k), this.DynModel.cov(k), this.ObsModel.obs(k), this.ObsModel.cov(k) etc., where k is the time index/interval)
%                 In general, all of the pre-defined DynamicModelX, ObservationModelX and ControlModelX classes comply with this rule, but is an important factor to consider when designing your own classes.
%                 See the DynamicModelX, ObservationModelX and ControlModelX templates for more details one how to make your models compatible with all ProjectX filters.
%
% KalmanFilterX Properties:
%    - Params   = structure with fields:
%       .k          = time index. Can also act as a time interval (Dt), depending on the underlying models. 
%       .x (*)      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state
%       .P (*)      = Estimated state covariance (P_{k|k}) - (nx x nx) matrix 
%       .xPred     = Predicted state mean (x_{k|k-1}) - (nx x 1) column vector
%       .PPred     = Predicted state mean (P_{k|k-1}) - (nx x nx) matrix
%       .y          = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
%       .y_Pred     = Predicted measurement mean (H*x_{k|k-1}) - (ny x 1) column vector
%       .S          = Innovation covariance (S_k) - (ny x ny) column vector
%       .K          = Kalman Gain (K_k) - (nx x ny) column vector
%   
%   - DynModel (*)  = Object handle to DynamicModelX SubClass     | (TO DO: LinearGaussDynModelX) 
%   - ObsModel (*)  = Object handle to ObservationModelX SubClass | (TO DO: LinearGaussObsModelX)
%   - CtrModel (*)  = Object handle to ControlModelX SubClass     | (TO DO: LinearCtrModelX)
%
%   (*) Signifies properties necessary to instantiate a class object
%
% KalmanFilterX Methods:
%    KalmanFilterX  - Constructor method
%    Predict        - Performs KF prediction step
%    Update         - Performs KF update step
%    Iterate        - Performs a complete KF iteration (Predict & Update)
%    Smooth         - Performs KF smoothing on a provided set of estimates
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
    properties
    end
    methods
        function this = KalmanFilterX(Init)
        % KalmanFilterX - Constructor method
        %   
        %   Inputs:
        %       Required
        %       ========
        %       DynModel  => DynamicModelX SubClass instance     
        %       ObsModel  => ObservationModelX SubClass instance   
        %               
        %       Optional
        %       ========
        %       CtrModel  => ControlModelX SubClass instance           
        %       x_init    => Prior state mean    
        %       P_init    => Prior state covariance
        %       u_init    => Prior control input 
        %       k         => Time variable k. If constant, k can be set 
        %                    once on initialisation and then left the same
        %                    for as long as it remains constant.
        %       
        %   
        %   Usage:
        %       kf = KalmanFilterX(DynModel, ObsModel, CtrModel, x_init, P_init, u_init, k);
        %       kf = KalmanFilterX(DynModel, ObsModel, CtrModel, 'x_init', x_init, 'P_init', P_init,  'u_init;, u_init, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
            
            % Call SuperClass method
            this@FilterX(Init);
            
%             % Add DynModel
%             if(~isfield(Init,'DynModel'))
%                 error('[KF] No DynModel provided!');
%             else
%                 this.DynModel = Init.DynModel;
%             end
%             
%             % Add ObsModel
%             if(~isfield(Init,'ObsModel'))
%                 error('[KF] No ObsModel provided!');
%             else
%                 this.ObsModel = Init.ObsModel;
%             end
%             
%             % Validate CtrModel
%             if(isfield(Init,'CtrModel'))
%                 this.CtrModel = Init.CtrModel;
%             end
            
            % Validate u_init
            if isfield(Init,'u_init')&&isfield(this,'CtrModel')
                this.Params.u = Init.u_init;
            end
        
            % Validate k
            if isfield(Init,'k') 
                this.Params.k = Init.k;
            end
            
            % Validate x_init
            if isfield(Init,'x_init') 
                this.Params.x = Init.x_init;
            end
            
            % Validate P_init
            if isfield(Init,'P_init') 
                this.Params.P = Init.P_init;
            end
        end
        
        function Predict(this)
        % Predict - Performs KF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and the control input "this.Params.u"
        %           need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.k = 1; % 1 sec)
        %       kf.Predict();
        %
        %   See also KalmanFilterX, Update, Iterate, Smooth.
            
            % Extract model parameters
            this.Params.F       = this.DynModel.sys(this.Params.k);
            this.Params.Q       = this.DynModel.sys_cov(this.Params.k);
            this.Params.H       = this.ObsModel.obs(this.Params.k);
            this.Params.R       = this.ObsModel.obs_cov(this.Params.k);
            if(~isempty(this.CtrModel))
                this.Params.B   = this.CtrModel.ctr(this.Params.k);
                this.Params.Qu  = this.CtrModel.ctr_cov(this.Params.k);
            else
                this.Params.u   = 0;
                this.Params.B   = 0;
                this.Params.Qu  = 0;
            end
            
            % Perform prediction
            [this.Params.xPred, this.Params.PPred, this.Params.yPred, this.Params.S, this.Params.Pxy] = ...
                KalmanFilterX_Predict(this.Params.x, this.Params.P, this.Params.F, this.Params.Q, ...
                                      this.Params.H, this.Params.R, this.Params.u, this.Params.B, this.Params.Qu);             
        end
        
        
        function Update(this)
        % Update - Performs KF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.y = y_new; % y_new is the new measurement)
        %       kf.Update(); 
        %
        %   See also KalmanFilterX, Predict, Iterate, Smooth.
        
            if(size(this.Params.y,2)>1)
                error('[KF] More than one measurement have been provided for update. Use KalmanFilterX.UpdateMulti() function instead!');
            elseif size(this.Params.y,2)==0
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                this.Params.x = this.Params.xPred;
                this.Params.P = this.Params.PPred;
                return;
            end
        
            % Perform single measurement update
            [this.Params.x, this.Params.P, this.Params.K] = ...
                KalmanFilterX_Update(this.Params.xPred,this.Params.PPred,...
                                     this.Params.y,this.Params.yPred,this.Params.S,this.Params.Pxy);
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
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Systems, vol. 29, no. 6, pp. 82-100, Dec. 2009.
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
        
        function Iterate(this)
        % Iterate - Performs a complete KF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.k = 1; % 1 sec)
        %       (kf.Params.y = y_new; % y_new is the new measurement)
        %       kf.Iterate();
        %
        %   See also KalmanFilterX, Predict, Update, Smooth.
        
            this.Predict();  % Predict         
            this.Update();   % Update
        end
        
        function smoothedEstimates = Smooth(this, filteredEstimates, interval)
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
        
        % ===============================>
        % ACCESS METHODS
        % ===============================>
        
        function [xPred, PPred] = getXpred(this)
        % getXpred - Returns the predicted state mean and covariance
        %   
        %   Inputs:
        %       N/A
        %   
        %   Outputs:   
        %       xPred : Last predicted state mean
        %       PPred : Last predicted state covariance
        %
        %   Usage:
        %       [xPred, PPred] = kf.getXpred()
        %   
        %   NOTE: At least one kf.Predict() call must have preceded.
        %
        %   See also KalmanFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract state prediction values from UKF 
            xPred = this.Params.xPred;
            PPred = this.Params.PPred;    
        end
        
        function [yPred, S] = getYpred(this)
        % getYpred - Returns the predicted measurement mean and innovation covariance
        %   
        %   Inputs:
        %       N/A
        %   
        %   Outputs:   
        %       yPred : Last predicted measurement mean
        %       S      : Last innovation covariance
        %
        %   Usage:
        %       [xPred, PPred] = kf.getXpred()
        %   
        %   NOTE: At least one kf.Predict() call must have preceded.
        %
        %   See also KalmanFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract measurement prediction values from UKF
            yPred = this.Params.yPred;
            S      = this.Params.S;    
        end
        
        function likelihoods = getObsLikelihoods(this,k,DataList,xPred)
        % getObsLikelihoods - Computes and returns the likelihood of a given series of observations, given a predicted state.
        %   
        %   Inputs:
        %       k        : time index/interval 
        %                  (Optional, default = this.Params.k)
        %       DataList : a (ny x Nm) matrix, where ny is the measurement dimensions and Nm is the number of measurements 
        %                  (Optional, default = this.Params.y)
        %       xPred   : a (nx x Ns) matrix, where nx is the state dimensionality and Ns is the number of predicted state samples
        %   
        %   Outputs:   
        %       Likelihoods : a (Nm x Ns) row vector, where each column corresponds to the likelihood of the respective measurement index
        %
        %   Usage:
        %       Likelihoods = getObsLikelihoods() returns the likelihoods of the internal this.Params.y measurements, using this.Params.k time index (Applicable only to time-variant measurement models. For time-invariant models, the time index has no effect)
        %       Likelihoods = getObsLikelihoods(k) or Likelihoods = getObsLikelihoods('k',k) returns the likelihoods of the internal this.Params.y measurements, using the provided k time index (applicable only to time-variant measurement models. For time-invariant models, the time index has no effect) 
        %       Likelihoods = getObsLikelihoods(k, DataList) or Likelihoods = getObsLikelihoods('k',k, 'DataList', DataList) returns the likelihoods of the DataList measurements, using the provided k time index (applicable only to time-variant measurement models. For time-invariant models, the time index has no effect)
        %       Likelihoods = getObsLikelihoods(k, DataList, xPred) or Likelihoods = getObsLikelihoods('k',k, 'DataList', DataList, 'xPred', xPred) returns the likelihoods of the DataList measurements, given the provided state mean(s) and using the provided k time index (applicable only to time-variant measurement models. For time-invariant models, the time index has no effect)
        %   NOTE: At least one kf.Predict() call must have preceded.
        %
        %   See also KalmanFilterX, Predict, Update, Smooth, Iterate.
            
            % Compute and return likelihoods
            likelihoods = this.ObsModel.eval_likelihood(k, DataList, xPred)';
        end
            
    end
end