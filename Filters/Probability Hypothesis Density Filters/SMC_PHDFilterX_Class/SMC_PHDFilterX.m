classdef SMC_PHDFilterX < handle
    % SMC_PHDFilterX class
    %
    % Summary of SMC_PHDFilterX:
    % This is a class implementation of a scaled Unscented Kalman Filter.
    %
    % SMC_PHDFilterX Properties:
    %    - Params       = structure, with fields:
    %       .type               = Type of PHD filter to use. (Currently only 'standard' and 'search' are implemented)
    %                             -> 'standard' is equivalent to generic SMC PHD [1]
    %                             -> 'search' is equivalent to [2], i.e. parameterised to detect, initialise and confirm targets
    %       .k                  = Time index (k) or time since last iteration (Dt)
    %       .Np                 = Number of particles
    %       .particles          = Particles (n_x by Np matrix, n_x being the dimensionality of the state vector)
    %       .w                  = Weights (1 by Np matrix)
    %       .pBirth             = Probability of birth (Only required if .birth_strategy='mixture')
    %       .pDeath             = Probability of death (1-e_k|k-1 in [1])
    %       .pConf              = Confirmation probability [2] (Only required if .type="search")
    %       .NpConf             = Number of particles for confirmed targets (Only required if .type="search")
    %       .NewTracksList      = List of tracks to be initiated [2] (Computed internally and only present if .type="search")
    %       .pDetect            = Probability of Detection 
    %       .Jk                 = Number of birth particles (Only required if .birth_strategy='expansion')
    %       .y                  = Measurements (n_y by Nm, Nm being the number of measurements)
    %       .lambda             = Clutter rate per unit volume (equivalent to ?_k(z) in [1] under the assumption of Poisson distributed clutter number)
    %       .resample_strategy  = Resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
    %       .birth_strategy     = Stategy for generating birth particles. Options are: i) 'expansion'; ii) 'mixture'; iii) 'obs_oriented'
    %       .gen_x0             = Birth intensity
    %       .gen_x1             = Observation oriented birth intensity (Only required if 'obs_oriented' birth_strategy has been selected)
    %   
    %   - DynModel (*)   = Object handle to Dynamic Model Class
    %   - ObsModel (*)   = Object handle to Observation Model Class
    %
    %   (*) Signifies properties necessary to instantiate a class thisect
    %
    % SMC_PHDFilterX Methods:
    %    SMC_PHDFilterX  - Constructor method
    %    Predict         - Performs SMC_PHD prediction step
    %    Update          - Performs SMC_PHD update step
    %    Iterate         - Performs a complete SMC_PHD iteration (Predict & Update)
    %    Smooth          - Performs SMC_PHD smoothing on a provided set of estimates

%
%       * Function Handles
%       -------------------
%       .sys;               = Transition function (Dynamics)
%       .sys_noise;         = Dynamics noise generator function
%       .likelihood         = Likelihood pdf p(z|x)
%       .obs_model          = Observation model (without noise), used to project state (particles) in measurement space
%       .gen_x0             = Birth particle sampling function 
%                            (Assumed to have a single input parameter, which is the number of particles to be generated)
%       .gen_x1 (*)         = Birth particle sampling function around measurements (Only required if .birth_strategy='obs_oriented')
%                            (Assumed to have a single input parameter, which is the number of particles to be generated)
%
%       Author: Lyudmil Vladimirov, University of Liverpool
%
% [1]  B. N. Vo, S. Singh and A. Doucet, "Sequential Monte Carlo methods for multitarget filtering with random finite sets," in IEEE Transactions on Aerospace and Electronic Systems, vol. 41, no. 4, pp. 1224-1245, Oct. 2005.
% [2]  P. Horridge and S. Maskell,  “Using a probabilistic hypothesis density filter to confirm tracks in a multi-target environment,” in2011 Jahrestagung der Gesellschaft fr Informatik, October 2011.
% [3]  B. ngu Vo and S. Singh, “Sequential monte carlo implementation of the phd filter for multi-targettracking,” inIn Proceedings of the Sixth International Conference on Information Fusion, pp. 792–799, 2003.
% =====================================================================================>
    properties
        Params
        DynModel
        ObsModel
    end
    
    methods
        function this = SMC_PHDFilterX(Init)
        % SMC_PHDFilterX - Constructor method
        %   
        %   Inputs:
        %       Required
        %       ========
        %       DynModel  => DynamicModelX SubClass instance     
        %       ObsModel  => ObservationModelX SubClass instance   
        %               
        %       Optional
        %       ========
        %       type               = Type of PHD filter to use. (Currently only 'standard' and 'search' are implemented)
        %                             -> 'standard' is equivalent to generic SMC PHD [1]
        %                             -> 'search' is equivalent to [2], i.e. parameterised to detect, initialise and confirm targets
        %       k                  = Time index (k) or time since last iteration (Dt)
        %       Np                 = Number of particles
        %       Jk                 = Number of birth particles (Only required if .birth_strategy='expansion')
        %       particles_init     = Particles (n_x by Np matrix, n_x being the dimensionality of the state vector)
        %       w_init             = Weights (1 by Np matrix)
        %       pBirth             = Probability of birth (Only required if .birth_strategy='mixture')
        %       pDeath             = Probability of death (1-e_k|k-1 in [1])
        %       pConf              = Confirmation probability [2] (Only required if .type="search")
        %       pDetect            = Probability of Detection 
        %       lambda             = Clutter rate per unit volume (equivalent to k(z) in [1] under the assumption of Poisson distributed clutter number)
        %       resample_strategy  = Resampling strategy. Set it either to 'multinomial_resampling' or 'systematic_resampling'
        %       birth_strategy     = Stategy for generating birth particles. Options are: i) 'expansion'; ii) 'mixture'; iii) 'obs_oriented'
        %       gen_x0             = Birth intensity
        %       gen_x1             = Observation oriented birth intensity (Only required if 'obs_oriented' birth_strategy has been selected)
        %       
        %   
        %   Usage:
        %       phd = UKalmanFilterX(DynModel, ObsModel, CtrModel, x_init, P_init, u_init, k);
        %       ukf = UKalmanFilterX(DynModel, ObsModel, CtrModel, 'x_init', x_init, 'P_init', P_init,  'u_init;, u_init, 'k', k);
        %
        %   See also Predict, Update, Iterate, Smooth.    
                   
            % Add DynModel
            if ~isfield(Init, 'DynModel')
                error('No dynamic model provided!');
            end
            this.DynModel = Init.DynModel;
            
            % Add ObsModel
            if ~isfield(Init, 'ObsModel')
                error('No observation model provided!');
            end
            this.ObsModel = Init.ObsModel;
            
            % Validate type
            if ~isfield(Init, 'type') 
                this.Params.type = 'standard';
            else
                this.Params.type = Init.type;
            end
            
            if(strcmp(this.Params.type, 'search'))
                % Validate .pConf
                if ~isfield(Init, 'pConf')
                    fprintf('[PHD] Probability pConf missing... Assumming pConf = 0.9.\n');
                    this.Params.pConf = 0.9;
                else
                    this.Params.pConf = Init.pConf;
                end
            
                % Validate .NpConf
                if ~isfield(Init, 'NpConf')
                    fprintf('[PHD] Number of confirm particles missing... Assumming NpConf = 1000.\n');
                    this.Params.NpConf = 1000;
                else
                    this.Params.NpConf = Init.NpConf;
                end
            end
        
            % Validate k
            if ~isfield(Init, 'k')
                this.Params.k = 1;
            else
                this.Params.k = Init.k;
            end
            
            % Validate .Np
            if ~isfield(Init, 'Np') && isfield(Init, 'particles_init')
                fprintf('[PHD] Number of particles missing... Assuming "Np = size(particles,2)"..\n');
                this.Params.Np = size(Init.particles_init,2);
            elseif ~isfield(Init, 'Np')
                fprintf('[PHD] Number of particles missing... Assuming "Np = 1000"..\n');
                this.Params.Np = 1000;
            else
                this.Params.Np = Init.Np;
            end
            
            % Validate .particles
            if (~isfield(Init, 'particles_init'))
                fprintf('[PHD] Particles not given... Proceeding to generation of initial particles..\n');
                if ~isfield(Init, 'gen_x0')
                    fprintf('[PHD] Function handle to sample from initial pdf not given... Cannot proceed..\n');
                    error('[PHD] Please supply either an initial set of particles, or a function handle (gen_0) to allow for generation of initial ones!\n');
                else
                    this.Params.particles = Init.gen_x0(this.Params.Np)'; % at time k=1
                    this.Params.w = repmat(1/Init.Np, 1, Init.Np);
                    fprintf('[PHD] Generated %d particles with uniform weights\n',this.Params.Np);
                end
            else
                if size(Init.particles_init,2)~=this.Params.Np
                    error('[PHD] Given number of particles (Np) is different that the size of supplied particle list! Aborting..\n');
                end
                this.Params.particles = Init.particles_init;
            end
            
            % Validate .w
            if ~isfield(Init, 'w_init')
                fprintf('[PHD] Initial set of weights not given... Proceeding to auto initialisation!\n');
                this.Params.w = repmat(1/this.Params.Np, 1, this.Params.Np);
                fprintf('[PHD] Uniform weights for %d particles have been created\n', this.Params.Np);
            else
                if (all(Init.w_init ==0))
                    fprintf('[PHD] Initial set of weights given as all zeros... Proceeding to auto initialisation!\n');
                    this.Params.w = repmat(1/this.Params.Np, 1, this.Params.Np);
                    fprintf('[PHD] Uniform weights for %d particles have been created\n', this.Params.Np); 
                else
                    this.Params.w = Init.w_init;
                end   
            end
            
            % Validate .resample_strategy
            if ~isfield(Init, 'resampling_strategy')
                fprintf('[PHD] Resampling strategy not given... Assuming "resampling_strategy = systematic_resampling"..\n');
                this.Params.resampling_strategy = 'systematic_resampling';
            else
                this.Params.resampling_strategy = Init.resampling_strategy;
            end 
            
            % Validate .birth_strategy
            if ~isfield(Init, 'birth_strategy')
                fprintf('[PHD] Birth strategy not given... Assuming "birth_strategy = mixture"..\n');
                this.Params.birth_strategy = 'mixture';
            else
                this.Params.birth_strategy = Init.birth_strategy;
            end 
            
            % Validate .Jk
            if (strcmp(this.Params.birth_strategy,'expansion') && ~isfield(Init, 'Jk'))
                error('[PHD] Birth strategy is set to "expansion" but number of birth particles (Jk) is not supplied... Aborting...');
            elseif strcmp(this.Params.birth_strategy,'expansion')
                this.Params.Jk = Init.Jk;             
            end
            
            % Validate .pBirth
            if (strcmp(this.Params.birth_strategy,'mixture') && ~isfield(Init, 'pBirth'))
                error('[PHD] Birth strategy is set to "mixture" but birth probability (pBirth) is not supplied... Aborting...');
            else
                this.Params.pBirth = Init.pBirth;
            end
            
            % Validate .pDeath
            if ~isfield(Init, 'pDeath')
                fprintf('[PHD] Probability pDeath missing... Assumming pDeath = 0.005.\n');
                this.Params.pDeath = 0.005;
            else
                this.Params.pDeath = Init.pDeath;
            end
            
            % Validate .pDetect
            if ~isfield(Init, 'pDetect')
                fprintf('[PHD] Probability pDetect missing... Assumming pDetect = 0.9.\n');
                this.Params.pDetect = 0.9;
            else
                this.Params.pDetect = Init.pDetect;
            end
            
                       
            % Validate .lambda
            if ~isfield(Init, 'lambda')
                fprintf('[PHD] Clutter density (lambda) missing... Assumming lambda = 1.\n');
                this.Params.lambda = 1;
            else
                this.Params.lambda = Init.lambda;
            end
            
            % Validate .gen_x0
            if ~isfield(Init, 'gen_x0')
                error('[PHD] Birth density (gen_x0) missing!');
            else
                this.Params.gen_x0 = Init.gen_x0;
            end
            
            % Validate .gen_x1
            if ~isfield(Init, 'gen_x1') && strcmp(this.Params.birth_strategy,'obs_oriented')
                error('[PHD] Birth strategy is set to "mixture" but Birth density (gen_x1) is not supplied... Aborting...');
            elseif strcmp(this.Params.birth_strategy,'obs_oriented')
                this.Params.gen_x1 = Init.gen_x1;
            end
            
        end
        
        % Predict function
        % ----------------
        % Performs the relevant SMC PHD prediction algo, based on the selected .type
        function Predict(this)
            switch this.Params.type
                case 'standard' % [1]
                    this.Predict_Standard();
                case 'search' % [2]
                    this.Predict_Search();
            end
        end
        
        % Update function
        % ----------------
        % Performs the relevant SMC PHD update algo, based on the selected .type
        function Update(this)
            switch this.Params.type
                case 'standard' % [1]
                    this.Update_Standard();
                case 'search' % [2]
                    this.Update_Search();
            end
        end
        
        % Predict_Standard function
        % ----------------
        % Performs the "standard" prediction step, as dictated by [1]
        function Predict_Standard(this)
                        
            % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.Params.birth_strategy, 'expansion')) 
                % Expansion method is equivalent to Eqs. (25-26) of [1]
                % assuming that:
                %  -> e_k|k-1(?) = 1-pDeath
                %  -> b_k|k-1(x|?) = 0 (No spawned targets)
                %  -> q_k(x_k|x_k-1,Z_k) = f_k|k-1(x|?)  
                %  i.e. (25) = (1-pDeath)*w_k-1^i
                % and
                %  -> ?_k(x_k) = 0.2*p_k(x_k|Z_k)
                %  i.e  (26) = 0.2/Jk
            
                % Expand number of particles to accomodate for births
                this.Params.particles = [this.Params.particles, zeros(size(this.Params.particles, 1), this.Params.Jk)]; 
                this.Params.w = [this.Params.w, zeros(1, this.Params.Jk)];
                this.Params.Np_total = this.Params.Np + this.Params.Jk;  

                % Generate Np normally predicted particles
                this.Params.particles(:,1:this.Params.Np) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:this.Params.Np), this.DynModel.sys_noise(this.Params.k, this.Params.Np)); % Simply propagate all particles
                this.Params.w(:,1:this.Params.Np) = (1-this.Params.pDeath)* this.Params.w(:,1:this.Params.Np);

                % Generate birth particles 
                this.Params.particles(:,this.Params.Np+1:end) = this.Params.gen_x0(this.Params.Jk)';
                this.Params.w(:,this.Params.Np+1:end) = 0.2/this.Params.Jk;
                
            elseif(strcmp(this.Params.birth_strategy, 'mixture'))
                % Mixture method is equivalent to the one proposed in Section 5 of [2]
                
                % Compute mixture components
                a = this.Params.pBirth;
                b = (1-this.Params.pDeath);
                Np_n = binornd(this.Params.Np,b/(a+b)); % Number of normally predicted particles
                Np_b = this.Params.Np - Np_n; % Number of birth particles

                % Generate normally predicted particles 
                if(Np_n)
                    this.Params.particles(:,1:Np_n) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:Np_n), this.DynModel.sys_noise(this.Params.k, Np_n)); % Simply propagate all particles
                end

                % Generate birth particles 
                if(Np_b>0)
                    this.Params.particles(:,Np_n+1:end) = this.Params.gen_x0(Np_b)';
                end
                this.Params.w(:, Np_n+1:end) = 0.2/(this.Params.Np-Np_n); % Assign weights to birth particles
                
                this.Params.Np_total = this.Params.Np;
                
            elseif(strcmp(this.Params.birth_strategy, 'obs_oriented'))
                % =========================================================>
                % UNDER DEVELOPMENT: Still in beta version. Not recommended
                % =========================================================>
                this.Params.Jk = size(this.Params.z,2)*100; % 100 particles are assigned per measurement
                
                % Expand number of particles to accomodate for births
                this.Params.particles = [this.Params.particles, zeros(4, this.Params.Jk)]; 
                this.Params.w = [this.Params.w, zeros(1, this.Params.Jk)];
                this.Params.Np_total = this.Params.Np + this.Params.Jk;  

                % Generate Np normally predicted particles
                this.Params.particles(:,1:this.Params.Np) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:this.Params.Np), this.DynModel.sys_noise(this.Params.k,this.Params.Np)); % Simply propagate all particles
                this.Params.w(:,1:this.Params.Np) = (1-this.Params.pDeath)* this.Params.w(:,1:this.Params.Np);

                % Generate birth particles 
                for i=1:size(this.Params.z,2)
                    this.Params.particles(:,this.Params.Np+1+(i-1)*100:this.Params.Np+(i-1)*100+100) = this.Params.gen_x0(this.Params.z(:,i), 100)';
                end
                %this.Params.particles(:,this.Params.Np+1:end) = this.Params.gen_x0(this.Params.Jk)';
                this.Params.w(:,this.Params.Np+1:end) = 0.2/this.Params.Jk;
            else
                error('Birth strategy "%s" not defined', this.Params.birth_strategy); 
            end
            
        end
        
        % Update_Standard function
        % ----------------
        % Performs the "standard" update step, as dictated by [1]
        function Update_Standard(this)
            
            
            % Compute g(z|x) matrix as in [1] 
            %this.Params.g = zeros(this.Params.Np_total,size(this.Params.y, 2));
            %for i = 1:size(this.Params.y, 2)
            this.Params.g = this.ObsModel.eval(this.Params.k, this.Params.y, this.Params.particles)';
            %end
            
            % Compute C_k(z) Eq. (27) of [1]  
            Ck = zeros(1,size(this.Params.y,2));
            for i = 1:size(this.Params.y,2)   % for all measurements
                Ck(i) = sum(this.Params.pDetect*this.Params.g(:,i)'.*this.Params.w,2);
            end
            this.Params.Ck = Ck;
            
            % Update weights Eq. (28) of [1]
            this.Params.w = (1-this.Params.pDetect + sum(this.Params.pDetect*this.Params.g./(ones(this.Params.Np_total,1)*(this.Params.lambda+this.Params.Ck)),2))'.*this.Params.w;
            
            % Resample (equivalent to Step 3 of [1]
            this.Params.Nk = sum(this.Params.w,2); % Compute total mass
            [this.Params.particles, this.Params.w] = this.resample(this.Params.particles, (this.Params.w/this.Params.Nk), this.Params.resampling_strategy, this.Params.Np); % Resample
            this.Params.w = this.Params.w*this.Params.Nk; % Rescale

        end
        
        % Predict_Search function
        % ----------------
        % Performs the "search" prediction step, as dictated by [2]
        function Predict_Search(this)
                        
             % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.Params.birth_strategy, 'expansion')) 
                % Expansion method is equivalent to Eqs. (25-26) of [1]
                % assuming that:
                %  -> e_k|k-1(?) = 1-pDeath
                %  -> b_k|k-1(x|?) = 0 (No spawned targets)
                %  -> q_k(x_k|x_k-1,Z_k) = f_k|k-1(x|?)  
                %  i.e. (25) = (1-pDeath)*w_k-1^i
                % and
                %  -> ?_k(x_k) = 0.2*p_k(x_k|Z_k)
                %  i.e  (26) = 0.2/Jk
            
                % Expand number of particles to accomodate for births
                this.Params.particles = [this.Params.particles, zeros(size(this.Params.particles, 1), this.Params.Jk)]; 
                this.Params.w = [this.Params.w, zeros(1, this.Params.Jk)];
                this.Params.Np_total = this.Params.Np + this.Params.Jk;  

                % Generate Np normally predicted particles
                this.Params.particles(:,1:this.Params.Np) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:this.Params.Np), this.DynModel.sys_noise(this.Params.k, this.Params.Np)); % Simply propagate all particles
                this.Params.w(:,1:this.Params.Np) = (1-this.Params.pDeath)* this.Params.w(:,1:this.Params.Np);

                % Generate birth particles 
                this.Params.particles(:,this.Params.Np+1:end) = this.Params.gen_x0(this.Params.Jk)';
                this.Params.w(:,this.Params.Np+1:end) = 0.2/this.Params.Jk;
                
            elseif(strcmp(this.Params.birth_strategy, 'mixture'))
                % Mixture method is equivalent to the one proposed in Section 5 of [2]
                
                % Compute mixture components
                a = this.Params.pBirth;
                b = (1-this.Params.pDeath);
                Np_n = binornd(this.Params.Np,b/(a+b)); % Number of normally predicted particles
                Np_b = this.Params.Np - Np_n; % Number of birth particles

                % Generate normally predicted particles 
                if(Np_n)
                    this.Params.particles(:,1:Np_n) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:Np_n), this.DynModel.sys_noise(this.Params.k, Np_n)); % Simply propagate all particles
                end

                % Generate birth particles 
                if(Np_b>0)
                    this.Params.particles(:,Np_n+1:end) = this.Params.gen_x0(Np_b)';
                end
                this.Params.w(:, Np_n+1:end) = 0.2/(this.Params.Np-Np_n); % Assign weights to birth particles
                
                this.Params.Np_total = this.Params.Np;
                
           elseif(strcmp(this.Params.birth_strategy, 'obs_oriented'))
                % =========================================================>
                % UNDER DEVELOPMENT: Still in beta version. Not recommended
                % =========================================================>
                this.Params.Jk2 = size(this.Params.y,2)*10; % 10 particles are assigned per measurement
                this.Params.Np_total = this.Params.Np +this. Params.Jk + this.Params.Jk2; 
                
                % Expand number of particles to accomodate for births
                this.Params.particles = [this.Params.particles, zeros(4,  this.Params.Jk + this.Params.Jk2)]; 
                this.Params.w = [this.Params.w, zeros(1,  this.Params.Jk + this.Params.Jk2)]; 

                % Generate Np normally predicted particles
                this.Params.particles(:,1:this.Params.Np) = this.DynModel.sys(this.Params.k, this.Params.particles(:,1:this.Params.Np), this.DynModel.sys_noise(this.Params.Np)); % Simply propagate all particles
                this.Params.w(:,1:this.Params.Np) = (1-this.Params.pDeath)* this.Params.w(:,1:this.Params.Np);

                % Generate birth particles 
                for i=1:size(this.Params.y,2)
                    this.Params.particles(:,this.Params.Np+1+(i-1)*10:this.Params.Np+(i-1)*10+10) = this.Params.gen_x1(this.Params.y(:,i), 10)';
                end
                %Params.particles(:,Params.Np+1:end) = Params.gen_x0(Params.Jk)';
                this.Params.w(:,this.Params.Np+1:this.Params.Np+this.Params.Jk2) = 0.2/(this.Params.Jk+this.Params.Jk2);
                
                % Generate birth particles 
                this.Params.particles(:,this.Params.Np+this.Params.Jk2+1:end) = this.Params.gen_x0(this.Params.Jk)';
                this.Params.w(:,this.Params.Np+1:end) = 0.2/(this.Params.Jk+this.Params.Jk2);
            else
                error('Birth strategy "%s" not defined.. Choose between "expansion" or "mixture" strategies!', this.Params.birth_strategy); 
            end               
        end
        
        % Update_Search function
        % ----------------
        % Performs the "search" update step, as dictated by [2]
        function Update_Search(this)
            
            % Tranform particles to measurement space
            %trans_particles = this.ObsModel.obs(this.Params.particles(:,:)); 
            
            % Get rhi measurement weights (computed externally as in Eq. (16) in [2])
            %this.Params.rhi = this.Params.rhi==1;  %  ones(1,size(this.Params.y,2)) Assume all measurements are unused
            
            % Perform particle gating
            % =========================================================>
            % UNDER DEVELOPMENT: Still in beta version. Not recommended
            % =========================================================>
            %DistM = ones(Params.Np_total, size(Params.z, 2))*1000;
            %for i = 1:size(Params.z, 2)
            %    DistM(:,i) = mahalDist(trans_particles', Params.z(:,i), Params.R, 2);
            %    valid_particles(:,i) = DistM(:,i)<10;
            %end
            % % Compute g(z|x) matrix as in [1] 
            %Params.g = zeros(size(trans_particles,2),size(Params.z, 2));
            %for i = 1:size(Params.z, 2)
            %    Params.g(find(valid_particles(:,i)),i) = Params.likelihood(Params.k, trans_particles(1:2,find(valid_particles(:,i)))', Params.z(:,i)');
            %end
            
            % Compute g(z|x) matrix as in [1] 
            this.Params.g = this.ObsModel.eval(this.Params.k, this.Params.y, this.Params.particles)';

            % Compute C_k(z) Eq. (27) of [1]  
            Ck = zeros(1,size(this.Params.y,2));
            for i = 1:size(this.Params.y,2)   % for all measurements
                Ck(i) = sum(this.Params.pDetect*this.Params.g(:,i)'.*this.Params.w,2);
            end
            this.Params.Ck = Ck;

            % Calculate pi Eq. (21) of [2]
            this.Params.pi = zeros(1, size(this.Params.y,2));
            for j = 1:size(this.Params.y,2)
                this.Params.pi(j) = sum((this.Params.pDetect*this.Params.rhi(j)*this.Params.g(:,j)'/(this.Params.lambda+this.Params.Ck(j))).*this.Params.w,2);
            end
            
            % Update weights Eq. (28) of [1]
            w = zeros(size(this.Params.y,2)+1, this.Params.Np_total);
            w(1,:) = (1-this.Params.pDetect)*this.Params.w;
            for j = 1:size(this.Params.y,2)
                w(j+1,:) = (this.Params.pDetect*this.Params.rhi(j)*this.Params.g(:,j)'/(this.Params.lambda+this.Params.Ck(j))).*this.Params.w; 
            end
            %this.Params.w = (1-this.Params.pDetect + sum(this.Params.pDetect*this.Params.g./(ones(this.Params.Np_total,1)*(this.Params.lambda+this.Params.Ck)),2))'.*this.Params.w;
            
            % Select measurements to be used for spawning new tracks
            CritMeasurements = find(this.Params.pi>this.Params.pConf);
            this.Params.NewTracks = [];
            
            % Initiate new tracks
            for j = 1:size(CritMeasurements,2)
                MeasInd = CritMeasurements(j); % Index of measurement
                
                % Get particles and weights
                NewTrack.particles = this.Params.particles;
                NewTrack.w         = w(MeasInd+1,:);
                % Resample particles to ensure they are correctly localised
                % around the measurement
                Nk = sum(NewTrack.w,2);
                [NewTrack.particles, NewTrack.w] = this.resample(NewTrack.particles, (NewTrack.w/Nk), this.Params.resampling_strategy, this.Params.NpConf);
                NewTrack.ExistProb = this.Params.pi(MeasInd);
                this.Params.NewTracks{end+1} = NewTrack; 
            end
            
            % Select measurements which are not to be used for new tracks
            NonCritMeasurements = setdiff([1:size(this.Params.y, 2)], CritMeasurements);
            
            w(w(CritMeasurements,:)>0) = 0;
            
            % Rescale new particle weights, considering only non critical measurements
            this.Params.w = sum(w([1,NonCritMeasurements+1],:),1);
            %Params.w = (1-Params.pDetect + sum(Params.pDetect*Params.g./(ones(Params.Np_total,1)*(Params.lambda+Params.Ck)),2))'.*Params.w;
            
            % Resample (equivalent to Step 3 of [1]
            this.Params.Nk = sum(this.Params.w,2); % Compute total mass
            [this.Params.particles, this.Params.w] = this.resample(this.Params.particles, (this.Params.w/this.Params.Nk), this.Params.resampling_strategy, this.Params.Np);
            this.Params.w = this.Params.w*this.Params.Nk; % Rescale
          
        end
        
       % Resampling function
       % -------------------
        function [xk, wk, idx] = resample(this, xk, wk, resampling_strategy, Np_new)
            Np = length(wk);  % Np = number of particles
            switch resampling_strategy
               case 'multinomial_resampling'
                  with_replacement = true;
                  idx = randsample(1:Np, Np_new, with_replacement, wk);
                %{
                  THIS IS EQUIVALENT TO:
                  edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
                  edges(end) = 1;                 % get the upper edge exact
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, idx] = histc(sort(rand(Np,1)), edges);
                %}
               case 'systematic_resampling'
                  % this is performing latin hypercube sampling on wk
                  edges = min([0 cumsum(wk)],1); % protect against accumulated round-off
                  edges(end) = 1;                 % get the upper edge exact
                  u1 = rand/Np_new;
                  % this works like the inverse of the empirical distribution and returns
                  % the interval where the sample is to be found
                  [~, ~, idx] = histcounts(u1:1/Np_new:1, edges);
               otherwise
                  error('Resampling strategy not implemented\n')
            end
            xk = xk(:,idx);                    % extract new particles
            wk = repmat(1/Np_new, 1, Np_new);          % now all particles have the same weight
        end
    end
end