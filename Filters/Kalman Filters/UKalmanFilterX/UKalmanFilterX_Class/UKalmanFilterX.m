classdef UKalmanFilterX < KalmanFilterX
    % UKalmanFilterX class
    %
    % Summary of UKalmanFilterX:
    % This is a class implementation of a scaled Unscented Kalman Filter.
    %
    % UKalmanFilterX Properties:
    %    - Params       = structure, with fields:
    %       .k          = time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .x (*)      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state
    %       .P (*)      = Estimated state covariance (P_{k|k}) - (nx x nx) matrix 
    %       .x_pred     = Predicted state mean (x_{k|k-1}) - (nx x 1) column vector
    %       .P_pred     = Predicted state mean (P_{k|k-1}) - (nx x nx) matrix
    %       .y          = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .y_pred     = Predicted measurement mean (H*x_{k|k-1}) - (ny x 1) column vector
    %       .S          = Innovation covariance (S_k) - (ny x ny) column vector
    %       .K          = Kalman Gain (K_k) - (ny x ny) column vector
    %       .alpha      = Default 0.5 |
    %       .kappa      = Default 0   |=> UKF scaling parameters
    %       .beta       = Default 2   |
    %   
    %   - DynModel (*)   = Object handle to Dynamic Model Class
    %   - ObsModel (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies properties necessary to instantiate a class thisect
    %
    % UKalmanFilterX Methods:
    %    UKalmanFilterX  - Constructor method
    %    Predict         - Performs UKF prediction step
    %    Update          - Performs UKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs UKF smoothing on a provided set of estimates
    % 
    % UKalmanFilterX Example:
  
    properties
    end
    
    methods
        function this = UKalmanFilterX(Init)
        % UKalmanFilterX - Constructor method
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
        %       ukf = UKalmanFilterX(DynModel, ObsModel, CtrModel, x_init, P_init, u_init, k);
        %       ukf = UKalmanFilterX(DynModel, ObsModel, CtrModel, 'x_init', x_init, 'P_init', P_init,  'u_init;, u_init, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
                 
            % Call SuperClass method
            this@KalmanFilterX(Init);
            
            % Validate alpha, kappa, betta
            if ~isfield(Init,'alpha')
                disp('[UKF] No alpha provided.. Setting "alpha=0.5"...'); this.Params.alpha = 0.5; 
            else
                this.Params.alpha = Init.alpha;
            end
            if ~isfield(Init, 'kappa')
                disp('[UKF] No kappa provided.. Setting "kappa=0"...'); this.Params.kappa = 0; 
            else
                this.Params.kappa = Init.kappa;
            end
            if ~isfield(Init,'beta')
                disp('[UKF] No beta provided.. Setting "beta=2"...'); this.Params.beta = 2; 
            else
                this.Params.beta = Init.beta;
            end
        end
        
        function Predict(this)
        % Predict - Performs UKF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ukf.Params.k = 1; % 1 sec)
        %       ukf.Predict();
        %
        %   See also UKalmanFilterX, Update, Iterate, Smooth.
        
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
            [this.Params.xPred, this.Params.PPred, this.Params.yPred, this.Params.S, this.Params.Pxy] = ...
                UKalmanFilterX_Predict(this.Params.alpha, this.Params.kappa, this.Params.beta,...
                                       this.Params.x, this.Params.P, this.Params.f, this.Params.Q, ...
                                       this.Params.h, this.Params.R, this.Params.u, this.Params.b, this.Params.Qu);         
        end
        
        function Update(this)
        % Update - Performs UKF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ukf.Params.y = y_new; % y_new is the new measurement)
        %       ukf.Update(); 
        %
        %   See also UKalmanFilterX, Predict, Iterate, Smooth.
        
            % Call SuperClass method
            Update@KalmanFilterX(this)
  
        end
        
        function UpdatePDA(this, assocWeights)
        % UpdatePDA - Performs UKF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ukf.Params.y = y_new; % y_new is the new measurement)
        %       ukf.UpdateMulti(assocWeights);
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Systems, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also UKalmanFilterX, Predict, Iterate, Smooth, resample.
        
            % Call SuperClass method
            UpdatePDA@KalmanFilterX(this, assocWeights);
        end
        
        function Iterate(this)
        % Iterate - Performs a complete UKF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (ukf.Params.k = 1; % 1 sec)
        %       (ukf.Params.y = y_new; % y_new is the new measurement)
        %       ukf.Iterate();
        %
        %   See also UKalmanFilterX, Predict, Update, Smooth.
        
           % Call SuperClass method
            Iterate@KalmanFilterX(this);
        end
        
        function smoothedEstimates = Smooth(this, filteredEstimates)
        % Smooth - Performs UKF smoothing on a provided set of estimates
        %          (Based on [1])
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %   
        %   Outputs:
        %       smoothed_estimates: a copy of the input (1 x N) cell array filtered_estimates, where the .x and .P fields have been replaced with the smoothed estimates   
        %
        %   (Virtual inputs at each iteration)        
        %           -> filtered_estimates{k}.x          : Filtered state mean estimate at timestep k
        %           -> filtered_estimates{k}.P          : Filtered state covariance estimate at each timestep
        %           -> filtered_estimates{k+1}.x_pred   : Predicted state at timestep k+1
        %           -> filtered_estimates{k+1}.P_pred   : Predicted covariance at timestep k+1
        %           -> smoothed_estimates{k+1}.x        : Smoothed state mean estimate at timestep k+1
        %           -> smoothed_estimates{k+1}.P        : Smoothed state covariance estimate at timestep k+1 
        %       where, smoothed_estimates{N} = filtered_estimates{N} on initialisation
        %
        %   (NOTE: The filtered_estimates array can be accumulated by running "filtered_estimates{k} = ukf.Params" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       ukf.Smooth(filtered_estimates);
        %
        %   [1] S. SÄrkkÄ, "Unscented Rauch-Tung-Striebel Smoother," in IEEE Transactions on Automatic Control, vol. 53, no. 3, pp. 845-849, April 2008.
        %
        %   See also UKalmanFilterX, Predict, Update, Iterate.
        
            if(nargin==2)
                smoothedEstimates = UKalmanFilterX_SmoothRTS(filteredEstimates);
            else
                smoothedEstimates = UKalmanFilterX_SmoothRTS(filteredEstimates,interval);
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
        %       [x_pred, P_pred] = ukf.getXpred()
        %   
        %   NOTE: At least one ukf.Predict() call must have preceded.
        %
        %   See also UKalmanFilterX, Predict, Update, Smooth, Iterate.
        
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
        %       [x_pred, P_pred] = ukf.getXpred()
        %   
        %   NOTE: At least one ukf.Predict() call must have preceded.
        %
        %   See also UKalmanFilterX, Predict, Update, Smooth, Iterate.
        
            % Call SuperClass method
            [y_pred, S] = getYpred@KalmanFilterX(this);    
        end

    end
end