classdef EParticleFilterX < ParticleFilterX
    % EParticleFilterX class
    %
    % Summary of EParticleFilterX:
    % This is a class implementation of an Extended Particle Filter.
    % An ExtendedKalmanFilterX instance is used to generate the Optimal Proposal.
    %
    % EParticleFilterX Properties:
    %    - Params       = structure, with fields:
    %       .k                      = time index. Can also act as a time interval (Dt), depending on the underlying models. 
    %       .Np                (*1) = number of particles. Default = 1000
    %       .particles         (*2) = particles - (nx x Np [x T]) matrix (3rd dimension T is optional)
    %       .w                      = weights - (1 x Np [x T]) vector. Default = repmat(1/Np,Np) (3rd dimension T is optional)  
    %       .x                      = Estimated state mean (x_{k|k}) - (nx x 1) column vector, where nx is the dimensionality of the state (Weighted average avg(w*.particles))
    %       .y                      = Measurement (y_k) - (ny x 1) column vector, where ny is the dimensionality of the measurement
    %       .gen_x0            (*1) = function handle of a procedure that samples from the initial pdf p_x0
    %       .resample_strategy (*)  = resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
    %   
    %   - DynModel (*)   = Object handle to Dynamic Model Class
    %   - ObsModel (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies compulsory properties to instantiate a class thisect
    %   (*x) Signify valid combinations of optional properties required to instantiate a class thisect
    %
    % EParticleFilterX Methods:
    %    EParticleFilterX - Constructor method
    %    Predict         - Performs EKF prediction step (Particles are not propagated!)
    %    Update          - Performs EKF and PF update steps
    %    Iterate         - Performs a complete EPF iteration (Predict & Update)
    %    Smooth          - Performs EPF smoothing on a provided set of estimates
    % 
    % EParticleFilterX Example:

    properties
        ekf   % Instance of a EKalmanFilter, used to generate optimal proposal
    end
    
    methods
        function this = EParticleFilterX(Init)
        % EParticleFilterX - Constructor method
        %   
        %   Inputs:
        %       Required
        %       ========
        %       DynModel  => DynamicModelX SubClass instance     
        %       ObsModel  => ObservationModelX SubClass instance   
        %               
        %       Optional
        %       ========
        %       CtrModel            => ControlModelX SubClass instance 
        %       Np                  => Number of particles   | ===================================>
        %       gen_x0              => Prior state pdf       | Supply either a prior pdf with Np, or a 
        %       particles_init      => Prior state particles | combination of particles and weights
        %       w_init              => Prior state weights   | ===================================>
        %       u_init              => Prior control input
        %       resample_strategy   => Resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
        %       k                   => Time variable k. If constant, k can be set once on initialisation and then left the same
        %                              for as long as it remains constant.
        %   
        %   Usage:
        %       epf = EParticleFilterX(DynModel, ObsModel, CtrModel, particles_init, w_init, u_init, reampling_strategy, k);
        %       epf = EParticleFilterX(DynModel, ObsModel, CtrModel, 'particles_init', particles_init, 'w_init', w_init,
        %                              'u_init;, u_init, 'reampling_strategy', reampling_strategy, 'k', k);
        %           or
        %       epf = EParticleFilterX(DynModel, ObsModel, CtrModel, Np, gen_x0, u_init, reampling_strategy, k);
        %       epf = EParticleFilterX(DynModel, ObsModel, CtrModel, 'Np', Np, 'gen_x0', gen_x0, 'u_init;, u_init, 
        %                              'resampling_strategy', resampling_strategy, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
        
           % Call SuperClass method
           this@ParticleFilterX(Init);
           
           % Instantiate EKF
           this.ekf = EKalmanFilterX(Init);
        end
        
        function Predict(this)
        % Predict - Performs EPF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       Usage:
        %       (epf.Params.k = 1; % 1 sec)
        %       (epf.Params.y = y_new; % y_new is the new measurement)
        %       epf.Iterate();
        %
        %   See also EParticleFilterX, Update, Iterate, Smooth, resample.
        
            % Update EKF time
            this.ekf.Params.k = this.Params.k;
            
            % Compute EKF prior mean and covariance
            this.ekf.Params.x = sum(repmat(this.Params.w,size(this.Params.particles,1),1).*this.Params.particles,2);
            this.ekf.Params.P = weightedcov(this.Params.particles',this.Params.w');
            
            % Iterate EKF to obtain Optimal Proposal
            this.ekf.Predict();
            this.Params.xPred = this.ekf.Params.xPred;
            this.Params.PPred = this.ekf.Params.PPred;
        end
        
        function Update(this)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (epf.Params.y = y_new; % y_new is the new measurement)
        %       epf.Update(); 
        %
        %   See also EParticleFilterX, Predict, Iterate, Smooth, resample.
            
            % Update EKF measurement
            this.ekf.Params.y = this.Params.y;
            
            % Perform EKF update to obtain Optimal Proposal
            this.ekf.Update();
            
            % Sample from EKF proposal
            this.Params.particles = mvnrnd(this.ekf.Params.x', this.ekf.Params.P,this.Params.Np)'; % Sample from optimal proposal
            
            % Call SuperClass method
            Update@ParticleFilterX(this);             
        end
        
        function UpdatePDA(this, assocWeights, LikelihoodMatrix)
        % UpdatePDA - Performs EPF update step, for multiple measurements
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to
        %                       measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %       LikelihoodMatrix: a (Nm x Np) likelihood matrix, where Nm is the number of measurements and Np is the number of particles.
        %
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.y = y_new; % y_new is the new measurement)
        %       pf.Update(); 
        %
        %   See also EParticleFilterX, Predict, Iterate, Smooth, resample.
            % Update EKF measurement
            this.ekf.Params.y = this.Params.y;
            
            % Perform EKF update to obtain Optimal Proposal
            this.ekf.UpdatePDA(assocWeights);
            
            % Sample from EKF proposal
            this.Params.particles = mvnrnd(this.ekf.Params.x', this.ekf.Params.P,this.Params.Np)'; % Sample from optimal proposal
            
            if(exist('LikelihoodMatrix','var'))
                UpdatePDA@ParticleFilterX(this, assocWeights, LikelihoodMatrix); 
            else
                UpdatePDA@ParticleFilterX(this, assocWeights);
            end
        end
        
        function Iterate(this)
        % Iterate - Performs a complete bootstrap PF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (epf.Params.k = 1; % 1 sec)
        %       (epf.Params.y = y_new; % y_new is the new measurement)
        %       epf.Iterate();
        %
        %   See also EParticleFilterX, Predict, Update, Smooth, resample.
        
           % Call SuperClass method
           Iterate@ParticleFilterX(this);
        end
        
        function smoothed_estimates = Smooth(this, filtered_estimates)
        % Smooth - Performs FBS smoothing on a provided set of estimates
        %   
        %    [=====================]
        %    [!!! NOT DEVELOPED !!!] 
        %    [=====================]
        %
        %   See also EParticleFilterX, Predict, Update, Iterate.
            
            % [=====================]
            % [!!! NOT DEVELOPED !!!] 
            % [=====================]
            error("[EPF] EParticleFilterX Smooth function has not been developed yet!");
            
            % Call SuperClass method 
            % smoothed_estimates = Smooth@ParticleFilterX(this);
        end        
        
        function [xk, wk, idx] = resample(this, xk, wk, resampling_strategy)
        % resample - Performs particle resampling
        %   
        %   Inputs:
        %       xk : Initial particle set (nx x Np) matrix, where nx is the dimensionality of the state and Np is the number of particles
        %       wk : Initial weights (1 x Np) matrix
        %       resampling_strategy : resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
        %   
        %   Outputs:   
        %       xk : Resampled particle set (nx x Np) matrix, where nx is the dimensionality of the state and Np is the number of particles
        %       wk : Resampled (Uniform) weights (1 x Np) matrix
        %       idx : Indices of resampled particles (Map to old particles)
        %
        %   Usage:
        %       [xk, wk] = resample(this, xk, wk, resampling_strategy)
        %       [xk, wk, idx] = resample(this, xk, wk, resampling_strategy)
        %
        %   See also EParticleFilterX, Predict, Update, Smooth, Iterate.    
        
            % Call SuperClass method
            [xk, wk, idx] = resample@ParticleFilterX(this, xk, wk, resampling_strategy);
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
        %       [xPred, PPred] = ekf.getXpred()
        %   
        %   NOTE: At least one ekf.Predict() call must have preceded.
        %
        %   See also EParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract state prediction values from UKF 
            xPred = this.ekf.Params.xPred;
            PPred = this.ekf.Params.PPred;   
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
        %       [xPred, PPred] = ekf.getXpred()
        %   
        %   NOTE: At least one ekf.Predict() call must have preceded.
        %
        %   See also EParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract measurement prediction values from UKF
            yPred = this.ekf.Params.yPred;
            S      = this.ekf.Params.S;  
        end
    end
end