classdef UParticleFilterX < ParticleFilterX
    % UParticleFilterX class
    %
    % Summary of UParticleFilterX:
    % This is a class implementation of an Unscented Particle Filter.
    % An UnscentedKalmanFilterX instance is used to generate the Optimal Proposal.
    %
    % UParticleFilterX Properties:
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
    % UParticleFilterX Methods:
    %    UParticleFilterX - Constructor method
    %    Predict         - Performs UKF prediction step (Particles are not propagated!)
    %    Update          - Performs UKF and PF update steps (Particles are resampled)
    %    Iterate         - Performs a complete UPF iteration (Predict & Update)
    %    Smooth          - Performs UPF smoothing on a provided set of estimates
    % 
    % UParticleFilterX Example:

    properties
        ukf   % Instance of a UKalmanFilter, used to generate optimal proposal
    end
    
    methods
        function this = UParticleFilterX(Init)
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
        %       upf = UParticleFilterX(DynModel, ObsModel, CtrModel, particles_init, w_init, u_init, reampling_strategy, k);
        %       upf = UParticleFilterX(DynModel, ObsModel, CtrModel, 'particles_init', particles_init, 'w_init', w_init,
        %                              'u_init;, u_init, 'reampling_strategy', reampling_strategy, 'k', k);
        %           or
        %       upf = UParticleFilterX(DynModel, ObsModel, CtrModel, Np, gen_x0, u_init, reampling_strategy, k);
        %       upf = UParticleFilterX(DynModel, ObsModel, CtrModel, 'Np', Np, 'gen_x0', gen_x0, 'u_init;, u_init, 
        %                              'resampling_strategy', resampling_strategy, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
        
           % Call SuperClass method
           this@ParticleFilterX(Init);
           
           % Instantiate EKF
           this.ukf = UKalmanFilterX(Init);
        end
        
        function Predict(this)
        % Predict - Performs UPF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       Usage:
        %       (upf.Params.k = 1; % 1 sec)
        %       (upf.Params.y = y_new; % y_new is the new measurement)
        %       upf.Iterate();
        %
        %   See also UParticleFilterX, Update, Iterate, Smooth, resample.
        
            % Update UKF time 
            this.ukf.Params.k = this.Params.k;
                        
            % Compute UKF prior mean and covariance
            this.ukf.Params.x = sum(repmat(this.Params.w,size(this.Params.particles,1),1).*this.Params.particles,2);
            this.ukf.Params.P = weightedcov(this.Params.particles',this.Params.w');
            
            % Perform UKF prediction
            this.ukf.Predict();
            this.Params.xPred = this.ukf.Params.xPred;
            this.Params.PPred = this.ukf.Params.PPred;
        end
        
        function Update(this)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (upf.Params.y = y_new; % y_new is the new measurement)
        %       upf.Update(); 
        %
        %   See also UParticleFilterX, Predict, Iterate, Smooth, resample.
            
            % Update UKF measurement
            this.ukf.Params.y = this.Params.y;
            
            % Perform UKF update to obtain Optimal Proposal
            this.ukf.Update();
            
            % Sample from UKF proposal
            this.Params.particles = mvnrnd(this.ukf.Params.x', this.ukf.Params.P, this.Params.Np)'; % Sample from optimal proposal
            
            % Call SuperClass method
            Update@ParticleFilterX(this);         
        end
        
        function UpdatePDA(this, assocWeights, LikelihoodMatrix)
        % UpdatePDA - Performs UPF update step, for multiple measurements
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
        %       upf.UpdateMulti(assocWeights, LikelihoodMatrix); 
        %
        %   See also UParticleFilterX, Predict, Iterate, Smooth, resample.
        
            % Update EKF measurement
            this.ukf.Params.y = this.Params.y;
            
            % Perform EKF update to obtain Optimal Proposal
            this.ukf.UpdatePDA(assocWeights);
            
            % Sample from EKF proposal
            this.Params.particles = mvnrnd(this.ukf.Params.x', this.ukf.Params.P,this.Params.Np)'; % Sample from optimal proposal
            
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
        %       (upf.Params.k = 1; % 1 sec)
        %       (upf.Params.y = y_new; % y_new is the new measurement)
        %       upf.Iterate();
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, resample.
         
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
        %   See also UParticleFilterX, Predict, Update, Iterate.
            
            % [=====================]
            % [!!! NOT DEVELOPED !!!] 
            % [=====================]
            error("[UPF] UParticleFilterX Smooth function has not been developed yet!");
            
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
        %       [xk, wk] = upf.resample(xk, wk, resampling_strategy)
        %       [xk, wk, idx] = upf.resample(xk, wk, resampling_strategy)
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
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
        %       [xPred, PPred] = upf.getXpred()
        %   
        %   NOTE: At least one upf.Predict() call must have preceded.
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract state prediction values from UKF 
            xPred = this.ukf.Params.xPred;
            PPred = this.ukf.Params.PPred;    
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
        %       [xPred, PPred] = upf.getXpred()
        %   
        %   NOTE: At least one upf.Predict() call must have preceded.
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Extract measurement prediction values from UKF
            yPred = this.ukf.Params.yPred;
            S = this.ukf.Params.S;    
        end
    end
end