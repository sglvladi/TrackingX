classdef KalmanFilterX2 < matlab.mixin.Copyable % Handle class with copy functionality
    % KalmanFilterX class
    %
    % Summary of KalmanFilterX:
    % This is a class implementation of a vanilla Kalman Filter.
    %
    % KalmanFilterX Properties:
    %   
    %   - DynModel          = Object handle to a DynamicModelX object
    %   - ObsModel          = Object handle to an ObservationModelX object
    %   - CtrModel (*)      = Object handle to a ControlModelX object
    %   - Init (*)          = Initial settings structure with fields
    %       .timeIndex      = Time index. Can also act as a time interval (Dt), depending on the underlying models. Default = 1
    %       .initStateMean  = Prior distribution mean. Default = Empty (Copied in Config.stateMean on instantiation)
    %       .initStateCov   = Prior distribution mean. Default = Empty (Copied in Config.stateMean on instantiation)
    %       .initCtrInput   = Inital control input. Only applies if CtrModel has been provided. Default = 0;
    %
    %   Runtime parameters - Computed and modified internally, but can be directly accessed at any time
    %   ~~~~~~~~~~~~~~~~~~
    %   - Params   = structure with fields:
    %       .timeIndex          = Time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .stateMean          = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state
    %       .stateCov           = Estimated state covariance (P_{k|k}) - (nx x nx) matrix        
    %       .predictedStateMean = Predicted state mean (x_{k|k-1}) - (nx x 1) column vector
    %       .predictedStateCov  = Predicted state mean (P_{k|k-1}) - (nx x nx) matrix
    %       .obsList            = Measurement list (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .predictedObsMean   = Predicted measurement mean (H*x_{k|k-1}) - (ny x 1) column vector
    %       .innovationCov      = Innovation covariance (S_k) - (ny x ny) column vector
    %       .kamlanGain         = Kalman Gain (K_k) - (nx x ny) column vector
    %
    %   (*) Optional
    %
    % KalmanFilterX Methods:
    %    KalmanFilterX  - Constructor method
    %    predict        - Performs KF prediction step
    %    update         - Performs KF update step
    %    iterate        - Performs a complete KF iteration (Predict & Update)
    %    smooth         - Performs KF smoothing on a provided set of estimates
    % 
    % KalmanFilterX Example:
    
     properties
        Params
        DynModel
        ObsModel
        CtrModel
    end
    methods
        function obj = KalmanFilterX2(varargin)
        % KalmanFilterX - Constructor method
        %   
        %   Inputs:
        %       DynModel | 
        %       ObsModel | => Check class help for more details
        %       CtrModel | (Optional)
        %       Init     | (Optional)
        %   
        %   Usage:
        %       kf = KalmanFilterX(DynModel, ObsModel);
        %       kf = KalmanFilterX(DynModel, ObsModel, CtrModel);
        %       kf = KalmanFilterX(DynModel, ObsModel, Init); 
        %       kf = KalmanFilterX(DynModel, ObsModel, CtrModel, Init);
        %
        %   See also predict, update, iterate, smooth.
            
            if (~isa(varargin{1},'DynamicModelX'))
                error("[KF] Invalid dynamic model provided. Only DynamicModelX subclasses are supported. To develop your own model, type 'help DynamicModelX' for genereal guidelines.");
            end
            obj.DynModel = varargin{1}; % First argument is always the dynamic model
            if (~isa(varargin{2},'ObservationModelX'))
                error("[KF] Invalid observation model provided. Only ObservationModelX subclasses are supported. To develop your own model, type 'help ObservationModelX' for genereal guidelines.");
            end
            obj.ObsModel = varargin{2}; % Second argument is always the observation model
            switch nargin
                case 3
                    if (isa(varargin{3},'ControlModelX'))
                        obj.CtrModel = varargin{3};
                    else
                        Init = varargin{3};
                        if(isfield(Init, 'timeIndex'))
                            obj.Params.timeIndex = Init.timeIndex;
                        else
                            obj.Params.timeIndex = 1;
                        end
                        if(isfield(Init, 'initStateMean'))
                            obj.Params.stateMean = Init.initStateMean;
                        else
                            warning("[KF] No inital mean provided. It needs to be provided before running the filter!");
                        end
                        if(isfield(Init, 'initStateCov'))
                            obj.Params.stateMean = Init.initStateCov;
                        else
                            warning("[KF] No inital covariance provided. It needs to be provided before running the filter!");
                        end
                    end
                case 4
                    if (~isa(varargin{3},'ControlModelX'))
                        error("[KF] Invalid control model provided. Only ControlModelX subclasses are supported. To develop your own model, type 'help ControlModelX' for genereal guidelines.");
                    end
                    obj.CtrModel = varargin{3};
                    Init = varargin{4};
                    if(isfield(Init, 'timeIndex'))
                        obj.Params.timeIndex = Init.timeIndex;
                    else
                        obj.Params.timeIndex = 1;
                    end
                    if(isfield(Init, 'initStateMean'))
                        obj.Params.stateMean = Init.initStateMean;
                    else
                        warning("[KF] No inital mean provided. It needs to be provided before running the filter!");
                    end
                    if(isfield(Init, 'initStateCov'))
                        obj.Params.stateMean = Init.initStateCov;
                    else
                        warning("[KF] No inital covariance provided. It needs to be provided before running the filter!");
                    end
            end
        end
        
        function Predict(obj)
        % predict - Performs KF prediction step
        %   
        %   Inputs:
        %       N/A
        %
        %   Properties required:
        %       DynModel (sys, sys_cov functions)
        %       ObsModel (obs, obs_cov functions)
        %       Params:
        %           .stateMean
        %           .stateCov
        %
        %   Outputs:
        %       N/A
        %   
        %   Properties modified:
        %       Params:
        %           .predictedStateMean
        %           .predictedStateCov
        %           .predictedObsMean
        %           .innovationCov
        %
        %   (NOTE: The time index/interval "obj.Params.timeIndex" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.timeIndex = 1; % 1 sec)
        %       kf.predict();
        %
        %   See also KalmanFilterX, update, iterate, smooth.
            
            % Compute predicted state mean and covariance
            obj.Params.predictedStateMean = obj.DynModel.sys(obj.Params.timeIndex, obj.Params.stateMean) + obj.Params.B * obj.Params.u;
            obj.Params.predictedStateCov = obj.DynModel.sys_cov(obj.Params.timeIndex, obj.Params.stateCov);
            
            % Compute predicted measurement mean and covariance
            obj.Params.predictedObsMean  = obj.ObsModel.obs(obj.Params.timeIndex, obj.Params.predictedStateMean);
            obj.Params.innovationCov     = obj.ObsModel.obs_cov(obj.Params.k, obj.Params.predictedStateCov);
            
        end
        
        
        function update(obj)
        % update - Performs KF update step for a single measurement
        %   
        %   Inputs:
        %       N/A
        %
        %   Properties required:
        %       DynModel (sys, sys_cov functions)
        %       ObsModel (obs, obs_cov functions)
        %       Params:
        %           .predictedStateMean
        %           .predictedStateCov
        %           .predictedObsMean
        %           .innovationCov
        %           .obsList
        %
        %   Outputs:
        %       N/A
        %   
        %   Properties modified:
        %       Params:
        %           .kalmanGain
        %           .stateMean
        %           .stateCov
        %
        %   (NOTE: The measurement "obj.Params.obsList" needs to be updated, when necessary, before calling this method)        
        %   
        %   Usage:
        %       (kf.Params.obsList = y_new; % y_new is the new measurement)
        %       kf.update(); 
        %
        %   See also KalmanFilterX, predict, iterate, smooth.
        
            if(size(obj.Params.obsList,2)>1)
                error('[KF] More than one measurement have been provided for update. Use KalmanFilterX.UpdateMulti() function instead!');
            elseif size(obj.Params.obsList,2)==0
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                obj.Params.stateMean = obj.Params.predictedStateMean;
                obj.Params.stateCov = obj.Params.predictedStateCov;
                return;
            end
        
            % Compute Kalman gain
            obj.Params.kalmanGain = obj.Params.predictedStateCov * obj.ObsModel.obs(obj.Params.k)' / (obj.Params.innovationCov);

            % Compute filtered estimates
            obj.Params.stateMean = obj.Params.predictedStateMean + obj.Params.kalmanGain * (obj.Params.obsList - obj.Params.predictedObsMean);
            obj.Params.stateCov = obj.Params.P_pred - obj.Params.K * obj.ObsModel.Params.h(obj.Params.k) * obj.Params.P_pred;
        end
        
        function UpdateMulti(obj, assocWeights)
        % UpdateMulti - Performs KF update step, for multiple measurements
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %       LikelihoodMatrix: a (Nm x Np) likelihood matrix, where Nm is the number of measurements and Np is the number of particles.
        %
        %   (NOTE: The measurement "obj.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.y = y_new; % y_new is the new measurement)
        %       pf.Update(); 
        %
        %   See also ParticleFilterX, Predict, Iterate, Smooth, resample.
            ObsNum = size(obj.Params.y,2);  
            ObsDim = size(obj.Params.y,1); 
            
            if(~ObsNum)
                warning('[KF] No measurements have been supplied to update track! Skipping Update step...');
                obj.Params.x = obj.Params.x_pred;
                obj.Params.P = obj.Params.P_pred;
                return;
            end
            
            if(~exist('assocWeights','var'))
                warning('[KF] No association weights have been supplied to update track! Applying default "assocWeights = [0, ones(1,ObsNum)/ObsNum];"...');
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            
            % Compute Kalman gain
            innov_err      = obj.Params.y - obj.Params.y_pred(:,ones(1,ObsNum)); % error (innovation) for each sample
            obj.Params.K   = obj.Params.P_pred*obj.ObsModel.Params.h(obj.Params.k)'/obj.Params.S;  

            % update
            %Pc              = (eye(size(obj.Params.x,1)) - obj.Params.K*obj.ObsModel.Params.h(obj.Params.k))*obj.Params.P_pred;
            Pc              = obj.Params.P_pred - obj.Params.K*obj.Params.S*obj.Params.K';
            tot_innov_err   = innov_err*assocWeights(2:end)';
            Pgag            = obj.Params.K*((innov_err.*assocWeights(ones(ObsDim,1),2:end))*innov_err' - tot_innov_err*tot_innov_err')*obj.Params.K';
            
            obj.Params.x    = obj.Params.x_pred + obj.Params.K*tot_innov_err;  
            obj.Params.P    = assocWeights(1)*obj.Params.P_pred + (1-assocWeights(1))*Pc + Pgag;
        end
        
        function Iterate(obj)
        % Iterate - Performs a complete KF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "obj.Params.k" and measurement "obj.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (kf.Params.k = 1; % 1 sec)
        %       (kf.Params.y = y_new; % y_new is the new measurement)
        %       kf.Iterate();
        %
        %   See also KalmanFilterX, Predict, Update, Smooth.
        
            obj.Predict();  % Predict         
            obj.Update();   % Update
        end
        
        function smoothed_estimates = Smooth(obj, filtered_estimates)
        % Smooth - Performs KF smoothing on a provided set of estimates
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of obj.Params after each iteration
        %                            
        %   (NOTE: The filtered_estimates array can be computed by running "filtered_estimates{k} = kf.Params" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       kf.Smooth(filtered_estimates);
        %
        %   See also KalmanFilterX, Predict, Update, Iterate.
        
            % Allocate memory
            N                           = length(filtered_estimates);
            smoothed_estimates          = cell(1,N);
            smoothed_estimates{N}       = filtered_estimates{N}; 
            
            % Perform Rauch–Tung–Striebel Backward Recursion
            for k = N-1:-1:1
                smoothed_estimates{k}.C     = filtered_estimates{k}.P * obj.DynModel.Params.F(filtered_estimates{k+1}.k)' / filtered_estimates{k+1}.P_pred;
                smoothed_estimates{k}.x     = filtered_estimates{k}.x + smoothed_estimates{k}.C * (smoothed_estimates{k+1}.x - filtered_estimates{k+1}.x_pred);
                smoothed_estimates{k}.P     = filtered_estimates{k}.P + smoothed_estimates{k}.C * (smoothed_estimates{k+1}.P - filtered_estimates{k+1}.P_pred) * smoothed_estimates{k}.C';                            
            end
        end
            
    end
end