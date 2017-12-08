classdef EKalmanFilterX<KalmanFilterX
    % EKalmanFilterX class
    %
    % Summary of EKalmanFilterX:
    % This is a class implementation of a first-order Extended Kalman Filter.
    %
    % EKalmanFilterX Properties:
    %    - Params       = structurem with fields:
    %       .k          = time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .x (*)      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state
    %       .P (*)      = Estimated state covariance (P_{k|k}) - (nx x nx) matrix 
    %       .x_pred     = Predicted state mean (x_{k|k-1}) - (nx x 1) column vector
    %       .P_pred     = Predicted state mean (P_{k|k-1}) - (nx x nx) matrix
    %       .y          = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .y_pred     = Predicted measurement mean (H*x_{k|k-1}) - (ny x 1) column vector
    %       .S          = Innovation covariance (S_k) - (ny x ny) column vector
    %       .K          = Kalman Gain (K_k) - (ny x ny) column vector
    %       .Fjac       = Transition matrix Jacobian - (nx x nx) matrix (Computed internally)
    %       .Hjac       = Observation matrix Jacobian - (ny x ny) matrix (Computed internally)
    %   
    %   - DynModel (*)   = Object handle to Dynamic Model Class
    %   - ObsModel (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies properties necessary to instantiate a class object
    %
    % EKalmanFilterX Methods:
    %    EKalmanFilterX  - Constructor method
    %    Predict         - Performs EKF prediction step
    %    Update          - Performs EKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs EKF smoothing on a provided set of estimates
    %
    
    properties
    end
    
    methods
        function this = EKalmanFilterX(Init)
        % EKalmanFilterX - Constructor method
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
        %       kf = EKalmanFilterX(DynModel, ObsModel, CtrModel, x_init, P_init, u_init, k);
        %       kf = EKalmanFilterX(DynModel, ObsModel, CtrModel, 'x_init', x_init, 'P_init', P_init,  'u_init;, u_init, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
        
           % Call SuperClass method
            this@KalmanFilterX(Init);
        end
        
        function Predict(this)
        % Predict - Performs EKF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.Params.k = 1; % 1 sec)
        %       ekf.Predict();
        %
        %   See also EKalmanFilterX, Update, Iterate, Smooth.
        
            % Extract model parameters
            this.Params.f       = @(x) this.DynModel.sys(this.Params.k,x);
            this.Params.Q       = this.DynModel.sys_cov(this.Params.k);
            this.Params.h       = @(x) this.ObsModel.obs(this.Params.k,x);
            this.Params.R       = this.ObsModel.obs_cov(this.Params.k);
            if(~isempty(this.CtrModel))
                this.Params.b   = @(x) this.CtrModel.ctr(this.Params.k,x);
                this.Params.Qu  = this.CtrModel.ctr_cov(this.Params.k);
            else
                this.Params.u   = 0;
                this.Params.b   = @(x) 0;
                this.Params.Qu  = 0;
            end
            
            % Perform prediction
            [this.Params.xPred, this.Params.PPred, this.Params.yPred, this.Params.S, this.Params.Pxy, this.Params.F, this.Params.H, this.Params.B] = ...
                EKalmanFilterX_Predict(this.Params.x, this.Params.P, this.Params.f, this.Params.Q, ...
                                      this.Params.h, this.Params.R, this.Params.u, this.Params.b, this.Params.Qu);                                                                  
        end   
        
        function Update(this)
        % Update - Performs EKF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.Params.y = y_new; % y_new is the new measurement)
        %       ekf.Update(); 
        %
        %   See also EKalmanFilterX, Predict, Iterate, Smooth.
        
            % Call SuperClass method
            Update@KalmanFilterX(this);
        
        end
        
        function UpdatePDA(this, assocWeights)
        % UpdateMulti - Performs EKF update step, for multiple measurements
        %               Update is performed according to the generic (J)PDAF equations [1] 
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.Params.y = y_new; % y_new is the new measurement)
        %       ekf.UpdateMulti(assocWeights);
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Systems, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also EKalmanFilterX, Predict, Iterate, Smooth, resample.
        
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
            UpdatePDA@KalmanFilterX(this, assocWeights);
        end
        
        function Iterate(this)
        % Iterate - Performs a complete EKF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ekf.Params.k = 1; % 1 sec)
        %       (ekf.Params.y = y_new; % y_new is the new measurement)
        %       ekf.Iterate();
        %
        %   See also EKalmanFilterX, Predict, Update, Smooth.
        
           % Call SuperClass method
            Iterate@KalmanFilterX(this);
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
        %   See also EKalmanFilterX, Predict, Update, Iterate.
        
            % Call SuperClass method
            if nargin==2
                smoothedEstimates = Smooth@KalmanFilterX(this, filteredEstimates);
            else
                smoothedEstimates = Smooth@KalmanFilterX(this, filteredEstimates, interval); 
            end
        end
        
        % ===============================>
        % ACCESS METHODS
        % ===============================>
        
        function [x_pred, P_pred] = getXpred(this)
        % getXpred - Returns the predicted state mean and covariance
        %   
        %   Inputs:
        %       N/A
        %   
        %   Outputs:   
        %       x_pred : Last predicted state mean
        %       P_pred : Last predicted state covariance
        %
        %   Usage:
        %       [x_pred, P_pred] = ekf.getXpred()
        %   
        %   NOTE: At least one ekf.Predict() call must have preceded.
        %
        %   See also EKalmanFilterX, Predict, Update, Smooth, Iterate.
        
            % Call SuperClass method
            [x_pred, P_pred] = getXpred@KalmanFilterX(this);  
        end
        
        function [y_pred, S] = getYpred(this)
        % getYpred - Returns the predicted measurement mean and innovation covariance
        %   
        %   Inputs:
        %       N/A
        %   
        %   Outputs:   
        %       y_pred : Last predicted measurement mean
        %       S      : Last innovation covariance
        %
        %   Usage:
        %       [x_pred, P_pred] = ekf.getXpred()
        %   
        %   NOTE: At least one ekf.Predict() call must have preceded.
        %
        %   See also EKalmanFilterX, Predict, Update, Smooth, Iterate.
        
            % Call SuperClass method
            [y_pred, S] = getYpred@KalmanFilterX(this);    
        end
    end
    
end

