classdef ParticleFilterX < matlab.mixin.Copyable
    % ParticleFilterX class
    %
    % Summary of ParticleFilterX:
    % This is a class implementation of a bootstrap Particle Filter.
    %
    % ParticleFilterX Properties:
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
    % ParticleFilterX Methods:
    %    ParticleFilterX - Constructor method
    %    Predict         - Performs UKF prediction step
    %    Update          - Performs UKF update step
    %    Iterate         - Performs a complete EKF iteration (Predict & Update)
    %    Smooth          - Performs UKF smoothing on a provided set of estimates
    % 
    % ParticleFilterX Example:

    properties
        Params
        DynModel
        ObsModel
        CtrModel
    end
    
    methods
        function this = ParticleFilterX(Init)
        % ParticleFilterX - Constructor method
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
        %       pf = ParticleFilterX(DynModel, ObsModel, CtrModel, particles_init, w_init, u_init, reampling_strategy, k);
        %       pf = ParticleFilterX(DynModel, ObsModel, CtrModel, 'particles_init', particles_init, 'w_init', w_init,
        %                            'u_init;, u_init, 'reampling_strategy', reampling_strategy, 'k', k);
        %           or
        %       pf = ParticleFilterX(DynModel, ObsModel, CtrModel, gen_x0, u_init, reampling_strategy, k);
        %       pf = ParticleFilterX(DynModel, ObsModel, CtrModel, 'gen_x0', gen_x0, 'u_init;, u_init, 
        %                              'resampling_strategy', resampling_strategy, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.
            
            % Add DynModel
            if(~isfield(Init,'DynModel'))
                error('[KF] No DynModel provided!');
            else
                this.DynModel = Init.DynModel;
            end
            
            % Add ObsModel
            if(~isfield(Init,'ObsModel'))
                error('[KF] No ObsModel provided!');
            else
                this.ObsModel = Init.ObsModel;
            end
            
            % Validate CtrModel
            if(isfield(Init,'CtrModel'))
                this.CtrModel = Init.CtrModel;
            end
            
            % Validate u_init
            if isfield(Init,'u_init')&&isfield(this, 'CtrModel')
                this.Params.u = p.Results.u_init;
            end
        
            % Validate k
            if isfield(Init,'k') 
                this.Params.k = Init.k;
            else
                this.Params.k = 1;
            end
            
            % Validate .Np
            if ~isfield(Init,'Np') && isfield(Init,'particles_init')
                fprintf('[PF] Number of particles missing... Assuming "Np = size(particles,2)"..\n');
                this.Params.Np = size(Init.particles_init,2);
            elseif ~isfield(Init,'Np')
                fprintf('[PF] Number of particles missing... Assuming "Np = 1000"..\n');
                this.Params.Np = 1000;
            else
                this.Params.Np = Init.Np;
            end
            
            % Validate .gen_x0
            if isfield(Init,'gen_x0')
                this.Params.gen_x0 = Init.gen_x0;
            end
            
            % Validate .particles
            if (~isfield(Init,'particles_init'))
                fprintf('[PF] Particles not given... Proceeding to generation of initial particles..\n');
                if ~isfield(Init,'gen_x0')
                    fprintf('[PF] Function handle to sample from initial pdf not given... Cannot proceed..\n');
                    error('[PF] Please supply either an initial set of particles, or a function handle (gen_0) to allow for generation of initial ones!\n');
                else
                    this.Params.particles = Init.gen_x0(this.Params.Np)'; % at time k=1
                    this.Params.w = repmat(1/Init.Np, 1, Init.Np);
                    fprintf('[PF] Generated %d particles with uniform weights\n',this.Params.Np);
                end
            else
                if size(Init.particles_init,2)~=this.Params.Np
                    error('[PF] Given number of particles (Np) is different that the size of supplied particle list! Aborting..\n');
                end
                this.Params.particles = Init.particles_init;
            end
            
            % Validate .w
            if ~isfield(Init,'w_init')
                fprintf('[PF] Initial set of weights not given... Proceeding to auto initialisation!\n');
                this.Params.w = repmat(1/this.Params.Np, 1, this.Params.Np);
                fprintf('[PF] Uniform weights for %d particles have been created\n', this.Params.Np);
            else
                if (all(Init.w_init ==0))
                    fprintf('[PF] Initial set of weights given as all zeros... Proceeding to auto initialisation!\n');
                    this.Params.w = repmat(1/this.Params.Np, 1, this.Params.Np);
                    fprintf('[PF] Uniform weights for %d particles have been created\n', this.Params.Np); 
                else
                    this.Params.w = Init.w_init;
                end   
            end
            
            % Validate .resample_strategy
            if ~isfield(Init,'resampling_strategy')
                fprintf('[PF] Resampling strategy not given... Assuming "resampling_strategy = systematic_resampling"..\n');
                this.Params.resampling_strategy = 'systematic_resampling';
            else
                this.Params.resampling_strategy = Init.resampling_strategy;
            end 
             
            % State estimate based on particles
            this.Params.x(:,1) = sum(this.Params.w.*this.Params.particles,2);      
        end
        
        function Predict(this)
        % Predict - Performs bootstrap PF prediction step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.k = 1; % 1 sec)
        %       pf.Predict();
        %
        %   See also ParticleFilterX, Update, Iterate, Smooth, resample.
        
            % Extract model parameters
            this.Params.f  = @(x,wk) this.DynModel.sys(this.Params.k, x, wk); % Transition function
            this.Params.wk = this.DynModel.sys_noise(this.Params.k, this.Params.Np); % Process noise
        
            % Propagate particles through the dynamic model
            this.Params.particles = ParticleFilterX_Predict(this.Params.f,this.Params.particles,this.Params.wk);
        end
        
        function Update(this)
        % Update - Performs bootstrap PF update step
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.y = y_new; % y_new is the new measurement)
        %       pf.Update(); 
        %
        %   See also ParticleFilterX, Predict, Iterate, Smooth, resample.
        
            if(size(this.Params.y,2)>1)
                error('[PF] More than one measurement have been provided for update. Use ParticleFilterX.UpdateMulti() function instead!');
            elseif size(this.Params.y,2)==0
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
            end
            
            % Perform update
            [this.Params.particles, this.Params.w, this.Params.x] = ...
                ParticleFilterX_Update(@(y,x) this.ObsModel.eval(this.Params.k,y,x),this.Params.y,...
                                       this.Params.particles,this.Params.w, this.Params.resampling_strategy);            
        end
        
        function UpdatePDA(this, assocWeights, LikelihoodMatrix)
        % UpdatePDA - Performs bootstrap PF update step, for multiple measurements
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
        %   See also ParticleFilterX, Predict, Iterate, Smooth, resample.
            ObsNum = size(this.Params.y,2);  
            if(~ObsNum)
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
                return;
            end
            
            if(~exist('assocWeights','var'))
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            if(~exist('LikelihoodMatrix','var') && isfield(this.Params, 'LikelihoodMatrix'))
                LikelihoodMatrix = this.Params.LikelihoodMatrix;  
            elseif ~exist('LikelihoodMatrix','var')
                LikelihoodMatrix = this.ObsModel.eval(this.Params.k, this.Params.y , this.Params.particles);
            end
            
            % Perform update
            [this.Params.particles, this.Params.w, this.Params.x] = ...
                ParticleFilterX_UpdatePDA(@(y,x) this.ObsModel.eval(this.Params.k,y,x),this.Params.y,this.Params.particles,...
                                          this.Params.w, this.Params.resampling_strategy, assocWeights, LikelihoodMatrix);  
            clear this.Params.LikelihoodMatrix;
        end
        
        function Iterate(this)
        % Iterate - Performs a complete bootstrap PF iteration (Predict & Update)
        %   
        %   Inputs:
        %       N/A 
        %   (NOTE: The time index/interval "this.Params.k" and measurement "this.Params.y" need to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.k = 1; % 1 sec)
        %       (pf.Params.y = y_new; % y_new is the new measurement)
        %       pf.Iterate();
        %
        %   See also ParticleFilterX, Predict, Update, Smooth, resample.
           this.Predict();
           this.Update();
        end
        
        function smoothed_estimates = Smooth(this, filtered_estimates)
        % Smooth - Performs FBS smoothing on a provided set of estimates
        %           (Based on [1])
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %   
        %   Outputs:
        %       smoothed_estimates: a copy of the input (1 x N) cell array filtered_estimates, where the .x and .P fields have been replaced with the smoothed estimates   
        %
        %   (Virtual inputs at each iteration)        
        %           -> filtered_estimates{k}.particles          : Filtered state mean estimate at timestep k
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
        %   [1] Mike Klaas, Mark Briers, Nando de Freitas, Arnaud Doucet, Simon Maskell, and Dustin Lang. 2006. Fast particle smoothing: if I had a million particles. In Proceedings of the 23rd international conference on Machine learning (ICML '06). ACM, New York, NY, USA, 481-488.
        %
        %   See also ParticleFilterX, Predict, Update, Iterate.
        
            % Allocate memory
            N                           = length(filtered_estimates);
            smoothed_estimates          = cell(1,N);
            smoothed_estimates{N}       = filtered_estimates{N}; 
            
            % Perform Rauch–Tung–Striebel Backward Recursion
            for k = N-1:-1:1
                lik = this.DynModel.eval(filtered_estimates{k}.k, filtered_estimates{k+1}.particles, filtered_estimates{k}.particles);
                denom = sum(filtered_estimates{k}.w(ones(this.Params.Np,1),:).*lik,2)'; % denom(1,j)
                smoothed_estimates{k}.w = filtered_estimates{k}.w(1,:) .* sum(smoothed_estimates{k+1}.w(ones(this.Params.Np,1),:).*lik'./denom(ones(this.Params.Np,1),:),2)';
                smoothed_estimates{k}.particles =  filtered_estimates{k}.particles;
                smoothed_estimates{k}.x = sum(smoothed_estimates{k}.w.*smoothed_estimates{k}.particles,2);
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
        %       [x_pred, P_pred] = upf.getXpred()
        %   
        %   NOTE: At least one upf.Predict() call must have preceded.
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Compute mean and covariance of particles 
            x_pred = sum(this.Params.w.*this.Params.particles,2);
            P_pred = (this.Params.particles-x_pred)*diag(this.Params.w)*(this.Params.particles-x_pred)';
            P_pred = .5*(P_pred+P_pred'); % Force symmetry
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
        %       [x_pred, P_pred] = upf.getXpred()
        %   
        %   NOTE: At least one upf.Predict() call must have preceded.
        %
        %   See also UParticleFilterX, Predict, Update, Smooth, Iterate.
        
            % Project particles into measurement space and compute mean and covariance
            trans_parts = this.ObsModel.obs(this.Params.k, this.Params.particles, this.ObsModel.obs_noise(this.Params.k,this.Params.Np));
            y_pred      = sum(this.Params.w.*trans_parts,2);
            S = weightedcov(trans_parts', this.Params.w);
        end
        
    end
end