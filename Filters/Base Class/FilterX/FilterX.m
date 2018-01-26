classdef FilterX < matlab.mixin.Copyable % Handle class with copy functionality
% FilterX Base class
%
% Summary of FilterX:
% This is the Base class for all TrackingX filters.
% Any custom defined Filters should be derived from this FilterX base class. 
%
% FilterX Properties:
%   - Params   = structure with fields
%   - DynModel  = Object handle to DynamicModelX SubClass     
%   - ObsModel  = Object handle to ObservationModelX SubClass 
%   - CtrModel  = Object handle to ControlModelX SubClass     
%
%   (*) Signifies properties necessary to instantiate a class object
%
% KalmanFilterX Methods:
%    FilterX  - Constructor method
%    Predict        - Performs filter prediction step
%    Update         - Performs filter update/correction step
%    Iterate        - Performs a complete filter iteration (Predict & Update)
%    Smooth         - Performs smoothing on a provided set of estimates
% 
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    
     properties
        Params
        DynModel
        ObsModel
        CtrModel
    end
    methods
        function this = FilterX(Init)
        % FilterX - Constructor method
        %   
        %   Inputs:
        %       Init => Structure with fields:
        %
        %       Required
        %       ========
        %       .DynModel  => DynamicModelX SubClass instance     
        %       .ObsModel  => ObservationModelX SubClass instance   
        %               
        %       Optional
        %       ========
        %       .CtrModel  => ControlModelX SubClass instance
        %
        %   Usage:
        %       f = FilterX(Init);
        %
        %   See also Predict, Update, Iterate, Smooth.
            
            % Add DynModel
            if(~isfield(Init,'DynModel'))
                error('No DynModel provided!');
            else
                this.DynModel = Init.DynModel;
            end
            
            % Add ObsModel
            if(~isfield(Init,'ObsModel'))
                error('No ObsModel provided!');
            else
                this.ObsModel = Init.ObsModel;
            end
            
            % Validate CtrModel
            if(isfield(Init,'CtrModel'))
                this.CtrModel = Init.CtrModel;
            end
           
        end
        
        function Predict(this)
        % Predict - Performs prediction step
        %   
        %   Inputs:
        %       N/A  
        %   
        %   Usage:
        %       f.Predict();
        %
        %   See also FilterX, Update, Iterate, Smooth.
                        
        end
        
        
        function Update(this)
        % Update - Performs KF update step
        %   
        %   Inputs:
        %       N/A 
        %
        %   Usage:
        %       f.Update(); 
        %
        %   See also FilterX, Predict, Iterate, Smooth.
        
        end
        
        function Iterate(this)
        % Iterate - Performs a complete KF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %
        %   Usage:
        %       kf.Iterate();
        %
        %   See also FilterX, Predict, Update, Smooth.
        
            this.Predict();  % Predict         
            this.Update();   % Update
        end
        
        function smoothedEstimates = Smooth(this, filteredEstimates, interval)
        % Smooth - Performs KF smoothing on a provided set of estimates
        %   
        %   Inputs:
        %       filteredEstimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %          
        %   Usage:
        %       kf.Smooth(filteredEstimates);
        %
        %   See also FilterX, Predict, Update, Iterate.
                
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
        %       [xPred, PPred] = f.getXpred()
        %   
        %   NOTE: At least one f.Predict() call must have preceded.
        %
        %   See also FilterX, Predict, Update, Smooth, Iterate.
        
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
            S     = this.Params.S;    
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