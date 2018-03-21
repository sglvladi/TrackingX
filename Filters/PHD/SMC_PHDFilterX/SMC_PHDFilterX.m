classdef SMC_PHDFilterX < FilterX
% SMC_PHDFilterX class
%
% Summary of SMC_PHDFilterX:
% This is a class implementation of a Sequential Monte Carlo (SMC) Probabilistic
% Hypothesis Density (PHD) Filter.
%
% SMC_PHDFilterX Properties:
%   + NumParticles      The number of particles employed by the PHD Filter
%   + Particles         A (NumStateDims x NumParticles) matrix used to store 
%                       the last computed/set filtered particles  
%   + Weights           A (1 x NumParticles) vector used to store the weights
%                       of the last computed/set filtered particles
%   + PredParticles     A (NumStateDims x NumParticles) matrix used to store 
%                       the last computed/set predicted particles  
%   + PredWeights       A (1 x NumParticles) vector used to store the weights
%                       of the last computed/set predicted particles
%   + MeasurementList   A (NumObsDims x NumObs) matrix used to store the received 
%                       measurements
%   + ResamplingScheme  Method used for particle resampling, specified as 
%                       'multinomial', 'systematic'. Default = 'systematic'
%   + ResamplingPolicy  A (1 x 2) cell array, specifying the resampling trigger
%                       conditions. ReamplingPolicy{1} should be a string
%                       which can be either "TimeInterval", in which case 
%                       ReamplingPolicy{2} should specify the number of 
%                       iterations after which resampling should be done,
%                       or "EffectiveRatio", in which case Resampling{2}
%                       should specify the minimum ratio of effective particles
%                       which when reached will trigger the resampling process.
%                       Default ResamplingPolicy = {"TimeInterval",1}, meaning
%                       that resampling is performed on every iteration of
%                       the Particle Filter (upon update).                       
%   + Resampler         An object handle to a ResamplerX subclass. If a 
%                       Resampler is provided, then it will override any choice
%                       specified within the ResamplingScheme. ResamplingPolicy
%                       will not be affected.
%   + BirthScheme       A (1 x 3) cell array, specifying the particle birth
%                       scheme. BirthScheme{1} should be a string which can be
%                       set to either "Mixture", in which case BirthScheme{2}
%                       should specify the probability of birth of new particles,
%                       or "Expansion", in which case BirthScheme{2} should
%                       specify the number of particles to be "birthed" at
%                       each iteration of the filter.
%                       Default BirthScheme = {"Mixture",0.5} meaning that
%                       particles are birthed using the mixture scheme, with
%                       a birth probability of 50%.
%   + BirthIntFcn       A function handle, which when called generates a set 
%                       of initial particles and weights.
%   + ProbOfDeath       The probability that a target may cease to exist
%                       between consecutive iterations of the filter.
%   + ProbOfDetection   The probablity that a target will be detected in
%                       a given measurement scan.
%   + NumTargets        The estimated number of targets following an update step.
%   + Model             An object handle to StateSpaceModelX object
%       + Dyn = Object handle to DynamicModelX SubClass      
%       + Obs = Object handle to ObservationModelX SubClass 
%       + Ctr = Object handle to ControlModelX SubClass 
%
% SMC_PHDFilterX Methods:
%   + SMC_PHDFilterX  - Constructor method
%   + predict         - Performs SMC_PHD prediction step
%   + update          - Performs SMC_PHD update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1]  B. N. Vo, S. Singh and A. Doucet, "Sequential Monte Carlo methods for multitarget filtering with random finite sets," in IEEE Transactions on Aerospace and Electronic Systems, vol. 41, no. 4, pp. 1224-1245, Oct. 2005.
% [2]  P. Horridge and S. Maskell,  “Using a probabilistic hypothesis density filter to confirm tracks in a multi-target environment,” in2011 Jahrestagung der Gesellschaft fr Informatik, October 2011.
% [3]  B. ngu Vo and S. Singh, “Sequential monte carlo implementation of the phd filter for multi-targettracking,” inIn Proceedings of the Sixth International Conference on Information Fusion, pp. 792–799, 2003.
%
% See also ParticleFilterX, KalmanFilerX.

    properties (Access = private, Hidden)
        
    end
    properties
        NumParticles = 10000
        NumParticlesTotal
        NumMeasurements
        NumTargets = 0
        Particles     
        Weights     
        PredParticles   
        PredWeights    
        MeasurementList
        ResamplingScheme = 'Systematic'
        ResamplingPolicy = {'TimeInterval',1}                    
        Resampler        
        BirthScheme = {'Mixture',0.5};
        BirthIntFcn 
        ProbOfDeath  
        ProbOfDetection 
        ClutterRate 
        MeasLikelihood
        MeasWeights = 1;
        WeightsPerMeasurement
        e_k
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
    end
    
    methods
        function this = SMC_PHDFilterX(varargin)
        % SMC_PHDFilterX - Constructor method
        %   
        % DESCRIPTION: 
        % * phd = SMC_PHDFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * phd = SMC_PHDFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * phd = SMC_PHDFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * phd = SMC_PHDFilterX(ssm,priorDistFcn) returns an object
        %   handle, preconfigured with the provided StateSpaceModel object 
        %   handle ssm and the prior information about the state, provided  
        %   in the form of the priorDistFcn function.
        % * phd = SMC_PHDFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * NumParticles        (Scalar) The number of particles to be employed by the  
        %                       SMC PHD Filter. [Default = 10000]
        % * PriorParticles      (NumStateDims x NumParticles) The initial set of particles
        %                       to be used by the Particle Filter. These are copied into
        %                       the Particles property by the constructor.
        % * PriorWeights        (1 x NumParticles matrix) The initial set of weights to be used 
        %                       by the Particle Filter. These are copied into the Weights 
        %                       property by the constructor. [Default = 1/NumParticles]
        % * PriorDistFcn        (function handle) A function handle, which when called
        %                       generates a set of initial particles and weights, which
        %                       are consecutively copied into the Particles and Weights
        %                       properties respectively. The function should accept exactly ONE
        %                       argument, which is the number of particles to be generated and
        %                       return 2 outputs. If a PriorDistFcn is specified, then any values
        %                       provided for the PriorParticles and PriorWeights arguments are ignored.
        % * ResamplingScheme    (String) Method used for particle resampling, specified as 
        %                       'Multinomial', 'Systematic'. [Default = 'Systematic']
        % * ResamplingPolicy    (1 x 2 cell array) specifying the resampling trigger
        %                       conditions. ReamplingPolicy{1} should be a (String)
        %                       which can be either "TimeInterval", in which case 
        %                       ReamplingPolicy{2} should be a (Scalar) specifying the number  
        %                       of iterations after which resampling should be performed,
        %                       or "EffectiveRatio", in which case Resampling{2} should be
        %                       a (Scalar) specifying the minimum ratio of effective particles 
        %                       which, when reached, will trigger the resampling process
        %                       [Default ResamplingPolicy = {'TimeInterval',1}], meaning
        %                       that resampling is performed on every iteration of
        %                       the Particle Filter (upon update).                       
        % * Resampler           An object handle to a ResamplerX subclass. If a 
        %                       Resampler is provided, then it will override any choice
        %                       specified within the ResamplingScheme. ResamplingPolicy
        %                       will not be affected.
        %  * BirthScheme        (1 x 2 cell array) specifying the particle birth
        %                       scheme. BirthScheme{1} should be a string which can be
        %                       set to either "Mixture", in which case BirthScheme{2}
        %                       should specify the probability of birth of new particles,
        %                       or "Expansion", in which case BirthScheme{2} should
        %                       specify the number of particles to be "birthed" at
        %                       each iteration of the filter.
        %                       Default BirthScheme = {"Mixture",0.5} meaning that
        %                       particles are birthed using the mixture scheme, with
        %                       a birth probability of 50%.
        % * BirthIntFcn         (function handle) A function handle, which when called 
        %                       generates a set of initial particles and weights.
        % * ProbOfDeath         (Scalar) The probability that a target may cease to exist
        %                       between consecutive iterations of the filter.
        % * ProbOfDetection     (Scalar) The probablity that a target will be detected in
        %                       a given measurement scan.
        % * NumTargets          (Scalar) The estimated number of targets following an 
        %                       update step.
        % * Model               An object handle to StateSpaceModelX object
        %
        %  See also predict, update.   
            
            % Call SuperClass method
            this@FilterX(varargin{:});
            
            this.e_k = 0.5*ones(1,this.NumParticles);
            if(nargin==0)
                this.Resampler = SystematicResamplerX();
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'NumParticles'))
                        this.NumParticles  = config.NumParticles;
                    end
                    if (isfield(config,'PriorDistFcn'))
                        [this.Particles,this.Weights] = ...
                            config.PriorDistFcn(this.NumParticles);
                    elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                        this.Particles = config.priorParticles;
                        this.Weights = config.priorWeights;
                    end
                     if (isfield(config,'Resampler'))
                        this.Resampler = config.Resampler;
                     elseif (isfield(config,'ResamplingScheme'))
                        if(strcmp(config.ResamplingScheme,'Systematic'))
                            this.Resampler = SystematicResamplerX();
                        elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                            this.Resampler = MultinomialResamplerX();
                        end
                     else
                         this.Resampler = SystematicResamplerX();
                     end
                     if (isfield(config,'ResamplingPolicy'))
                         this.ResamplingPolicy = config.ResamplingPolicy;
                     end
                     if (isfield(config,'BirthScheme'))
                         this.BirthScheme = config.BirthScheme;
                     end
                     if (isfield(config,'BirthIntFcn'))
                         this.BirthIntFcn = config.BirthIntFcn;
                     end
                     if (isfield(config,'ProbOfDeath'))
                         this.ProbOfDeath = config.ProbOfDeath;
                     end
                     if (isfield(config,'ProbOfDetection'))
                         this.ProbOfDetection = config.ProbOfDetection;
                     end
                     if (isfield(config,'ClutterRate'))
                         this.ClutterRate = config.ClutterRate;
                     end
                end
                this.e_k = 0.5*ones(1,this.NumParticles);
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'NumParticles'))
                this.NumParticles  = config.NumParticles;
            end
            if (isfield(config,'PriorDistFcn'))
                [this.Particles,this.Weights] = ...
                    config.PriorDistFcn(this.NumParticles);
            elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                 this.Particles = config.PriorParticles;
                 this.Weights = config.PriotWeights;
            end
             if (isfield(config,'Resampler'))
                 this.Resampler = config.Resampler;
             elseif (isfield(config,'ResamplingScheme'))
                 if(strcmp(config.ResamplingScheme,'Systematic'))
                     this.Resampler = SystematicResamplerX();
                 elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                     this.Resampler = MultinomialResamplerX();
                 end
             else
                 this.Resampler = SystematicResamplerX();
             end
             if (isfield(config,'ResamplingPolicy'))
                 this.ResamplingPolicy = config.ResamplingPolicy;
             end
             if (isfield(config,'BirthScheme'))
                 this.BirthScheme = config.BirthScheme;
             end
             if (isfield(config,'BirthIntFcn'))
                 this.BirthIntFcn = config.BirthIntFcn;
             end
             if (isfield(config,'ProbOfDeath'))
                 this.ProbOfDeath = config.ProbOfDeath;
             end
             if (isfield(config,'ProbOfDetection'))
                 this.ProbOfDetection = config.ProbOfDetection;
             end
             if (isfield(config,'ClutterRate'))
                 this.ClutterRate = config.ClutterRate;
             end
             this.e_k = 0.5*ones(1,this.NumParticles);
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the SMC PHD Filter with a certain 
        % set of parameters.  
        %   
        % DESCRIPTION: 
        % * initialise(phd,ssm) initialises the SMC_PHDFilterX object 
        %   phd with the provided StateSpaceModelX object ssm.
        % * initialise(phd,priorParticles,priorWeights)initialises the 
        %   SMC_PHDFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(phd,ssm,priorDistFcn) initialises the SMC_PHDFilterX
        %   object pf with the provided StateSpaceModel object handle ssm
        %   and the prior information about the state, provided in the form 
        %   of the priorDistFcn function.
        % * initialise(phd,___,Name,Value,___) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value
        %   pair arguments.
        %
        % INPUT ARGUMENTS:
        % * NumParticles        (Scalar) The number of particles to be employed by the  
        %                       SMC PHD Filter. [Default = 10000]
        % * PriorParticles      (NumStateDims x NumParticles) The initial set of particles
        %                       to be used by the Particle Filter. These are copied into
        %                       the Particles property by the constructor.
        % * PriorWeights        (1 x NumParticles matrix) The initial set of weights to be used 
        %                       by the Particle Filter. These are copied into the Weights 
        %                       property by the constructor. [Default = 1/NumParticles]
        % * PriorDistFcn        (function handle) A function handle, which when called
        %                       generates a set of initial particles and weights, which
        %                       are consecutively copied into the Particles and Weights
        %                       properties respectively. The function should accept exactly ONE
        %                       argument, which is the number of particles to be generated and
        %                       return 2 outputs. If a PriorDistFcn is specified, then any values
        %                       provided for the PriorParticles and PriorWeights arguments are ignored.
        % * ResamplingScheme    (String) Method used for particle resampling, specified as 
        %                       'Multinomial', 'Systematic'. [Default = 'Systematic']
        % * ResamplingPolicy    (1 x 2 cell array) specifying the resampling trigger
        %                       conditions. ReamplingPolicy{1} should be a (String)
        %                       which can be either "TimeInterval", in which case 
        %                       ReamplingPolicy{2} should be a (Scalar) specifying the number  
        %                       of iterations after which resampling should be performed,
        %                       or "EffectiveRatio", in which case Resampling{2} should be
        %                       a (Scalar) specifying the minimum ratio of effective particles 
        %                       which, when reached, will trigger the resampling process
        %                       [Default ResamplingPolicy = {'TimeInterval',1}], meaning
        %                       that resampling is performed on every iteration of
        %                       the Particle Filter (upon update).                       
        % * Resampler           An object handle to a ResamplerX subclass. If a 
        %                       Resampler is provided, then it will override any choice
        %                       specified within the ResamplingScheme. ResamplingPolicy
        %                       will not be affected.
        %  * BirthScheme        (1 x 2 cell array) specifying the particle birth
        %                       scheme. BirthScheme{1} should be a string which can be
        %                       set to either "Mixture", in which case BirthScheme{2}
        %                       should specify the probability of birth of new particles,
        %                       or "Expansion", in which case BirthScheme{2} should
        %                       specify the number of particles to be "birthed" at
        %                       each iteration of the filter.
        %                       Default BirthScheme = {"Mixture",0.5} meaning that
        %                       particles are birthed using the mixture scheme, with
        %                       a birth probability of 50%.
        % * BirthIntFcn         (function handle) A function handle, which when called 
        %                       generates a set of initial particles and weights.
        % * ProbOfDeath         (Scalar) The probability that a target may cease to exist
        %                       between consecutive iterations of the filter.
        % * ProbOfDetection     (Scalar) The probablity that a target will be detected in
        %                       a given measurement scan.
        % * NumTargets          (Scalar) The estimated number of targets following an 
        %                       update step.
        % * Model               An object handle to StateSpaceModelX object
        %
        %  See also predict, update.   
                        
             this.e_k = 0.5*ones(1,this.NumParticles);
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'Model'))
                        this.Model = config.Model;
                    end
                    if (isfield(config,'NumParticles'))
                        this.NumParticles  = config.NumParticles;
                    end
                    if (isfield(config,'PriorDistFcn'))
                        [this.Particles,this.Weights] = ...
                            config.PriorDistFcn(this.NumParticles);
                    elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                         this.Particles = config.PriorParticles;
                         this.Weights = config.PriotWeights;
                    end
                    if (isfield(config,'Resampler'))
                        this.Resampler = config.Resampler;
                    elseif (isfield(config,'ResamplingScheme'))
                        if(strcmp(config.ResamplingScheme,'Systematic'))
                            this.Resampler = SystematicResamplerX();
                        elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                            this.Resampler = MultinomialResamplerX();
                        end
                    else
                        this.Resampler = SystematicResamplerX();
                    end
                    if (isfield(config,'ResamplingPolicy'))
                        this.ResamplingPolicy = config.ResamplingPolicy;
                    end
                    if (isfield(config,'BirthScheme'))
                        this.BirthScheme = config.BirthScheme;
                    end
                    if (isfield(config,'BirthIntFcn'))
                        this.BirthIntFcn = config.BirthIntFcn;
                    end
                    if (isfield(config,'ProbOfDeath'))
                        this.ProbOfDeath = config.ProbOfDeath;
                    end
                    if (isfield(config,'ProbOfDetection'))
                        this.ProbOfDetection = config.ProbOfDetection;
                    end
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'NumParticles'))
                this.NumParticles  = config.NumParticles;
            end
            if (isfield(config,'PriorDistFcn'))
                [this.Particles,this.Weights] = ...
                    config.PriorDistFcn(this.NumParticles);
            elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                 this.Particles = config.PriorParticles;
                 this.Weights = config.PriotWeights;
            end
             if (isfield(config,'Resampler'))
                 this.Resampler = config.Resampler;
             elseif (isfield(config,'ResamplingScheme'))
                 if(strcmp(config.ResamplingScheme,'Systematic'))
                     this.Resampler = SystematicResamplerX();
                 elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                     this.Resampler = MultinomialResamplerX();
                 end
             else
                 this.Resampler = SystematicResamplerX();
             end
             if (isfield(config,'ResamplingPolicy'))
                 this.ResamplingPolicy = config.ResamplingPolicy;
             end
             if (isfield(config,'BirthScheme'))
                 this.BirthScheme = config.BirthScheme;
             end
             if (isfield(config,'BirthIntFcn'))
                 this.BirthIntFcn = config.BirthIntFcn;
             end
             if (isfield(config,'ProbOfDeath'))
                 this.ProbOfDeath = config.ProbOfDeath;
             end
             if (isfield(config,'ProbOfDetection'))
                 this.ProbOfDetection = config.ProbOfDetection;
             end
        end
        
        function predict(this)
        % PREDICT Perform SMC PHD Filter prediction step
        %   
        % DESCRIPTION: 
        % * predict(this) calculates the predicted PHD
        %
        % MORE DETAILS:
        % * SMC_PHDFilterX uses the Model class property, which should 
        %   be an instance of the TrackingX.Models.StateSpaceModel class, in order
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
            
            % reset measLikelihood matrix
            this.pmeasLikelihood_ = [];
        
            % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.BirthScheme(1), 'Expansion')) 
                % Expansion method is equivalent to Eqs. (25-26) of [1]
                % assuming that:
                %  -> e_k|k-1(?) = 1-pDeath
                %  -> b_k|k-1(x|?) = 0 (No spawned targets)
                %  -> q_k(x_k|x_k-1,Z_k) = f_k|k-1(x|?)  
                %  i.e. (25) = (1-pDeath)*w_k-1^i
                % and
                %  -> ?_k(x_k) = 0.2*p_k(x_k|Z_k)
                %  i.e  (26) = 0.2/Jk
                
                Jk = this.BirthScheme{2};
            
                % Expand number of particles to accomodate for births
                this.PredParticles = [this.Particles, zeros(size(this.Particles, 1), Jk)]; 
                this.PredWeights = [this.Weights, zeros(1, Jk)];
                this.NumParticlesTotal = this.NumParticles + Jk;  

                % Generate Np normally predicted particles
                this.PredParticles(:,1:this.NumParticles) = ...
                    this.Model.Dyn.feval(this.Particles(:,1:this.NumParticles), ...
                        this.Model.Dyn.random(this.NumParticles)); 
                this.PredWeights(:,1:this.NumParticles) = ...
                    (1-this.ProbOfDeath)* this.Weights(:,1:this.NumParticles);

                % Generate birth particles 
                this.PredParticles(:,this.NumParticles+1:end) = this.BirthIntFcn(Jk);
                this.PredWeights(:,this.NumParticles+1:end) = 0.2/Jk;
                this.e_k = [this.e_k ones(1, Jk)*0.5];
                this.e_k = this.e_k*(1-this.ProbOfDeath);
                
            elseif(strcmp(this.BirthScheme(1), 'Mixture'))
                % Mixture method is equivalent to the one proposed in Section 5 of [2]
                
                % Compute mixture components
                a = this.BirthScheme{2};
                b = (1-this.ProbOfDeath);
                Np_n = binornd(this.NumParticles,b/(a+b)); % Number of normally predicted particles
                Np_b = this.NumParticles - Np_n; % Number of birth particles

                % Generate normally predicted particles 
                if(Np_n)
                    this.PredParticles(:,1:Np_n) = ...
                        this.Model.Dyn.feval(this.Particles(:,1:Np_n), ...
                            this.Model.Dyn.random(Np_n)); % Simply propagate all particles
                end

                % Generate birth particles 
                if(Np_b>0)
                    this.PredParticles(:,Np_n+1:this.NumParticles) = this.BirthIntFcn(Np_b);
                end
                this.PredWeights = this.Weights * (1-this.ProbOfDeath);
                this.PredWeights(:, Np_n+1:this.NumParticles) = 0.2/(Np_b); % Assign weights to birth particles
                %this.PredWeights = this.PredWeights 
                this.e_k(:, Np_n+1:this.NumParticles) = this.BirthScheme{2};
                this.NumParticlesTotal = this.NumParticles;
                this.e_k = this.e_k*(1-this.ProbOfDeath);% + this.BirthScheme{2};
            end
        end
        
        % Update function
        % ----------------
        % Performs the relevant SMC PHD update algo, based on the selected .type
        function update(this)
            
            this.NumMeasurements = size(this.MeasurementList,2);
            
            % Compute g(z|x) matrix as in [1] 
            g = this.MeasLikelihood;
                        
            % Compute C_k(z) Eq. (27) of [1]  
            Ck = zeros(1,this.NumMeasurements);
            for i = 1:this.NumMeasurements   % for all measurements
                Ck(i) = sum(this.ProbOfDetection*g(i,:).*this.PredWeights,2);
            end
            
            % Calculate w^{n,i} Eq. (21) of [2]
            this.WeightsPerMeasurement = zeros(this.NumMeasurements, this.NumParticlesTotal);
            for i = 1:this.NumMeasurements
                this.WeightsPerMeasurement(i,:) = (this.ProbOfDetection*g(i,:)/(this.ClutterRate+Ck(i))).*this.PredWeights;
            end
            this.WeightsPerMeasurement = this.MeasWeights'.*this.WeightsPerMeasurement;
            
            betta = zeros(this.NumMeasurements+1,this.NumParticlesTotal);
            betta(1,:) = (1-this.ProbOfDetection)/this.NumParticlesTotal;
            for j = 1:this.NumMeasurements
                betta(j+1,:) = this.ProbOfDetection*g(j,:)/this.ClutterRate;
            end
            Lj = 1 + (1-this.ProbOfDetection) + sum(betta(2:end,:),1);
%             betta = betta./(sum(betta,1));
%             lambda = betta(1,:)/this.NumMeasurements;
%             e_ij = (lambda+betta(2:end,:)).*this.e_k./(this.e_k + betta(1,:).*(1-this.e_k)/(1-this.ProbOfDetection));
%             m_ij = mean(e_ij,2);      
            try
                this.e_k = this.e_k./(this.e_k + (1-this.e_k)./Lj);
            catch
                
            end
            % Calculate pi Eq. (21) of [2]
            p_i = zeros(1,this.NumMeasurements);
            for i = 1:this.NumMeasurements
                p_i(i)= sum(this.WeightsPerMeasurement(i,:),2);
                if(p_i(i)>0.8)
                    %warning('Detected target %d! %f - %f = %f %%',i,p_i(i),sum(this.WeightsPerMeasurement(i,:).*this.e_k), (p_i(i)-sum(this.WeightsPerMeasurement(i,:).*this.e_k))/p_i(i)*100);
                end
            end
            
            % Update weights Eq. (28) of [1]
            this.Weights = (1-this.ProbOfDetection)*this.PredWeights + sum(this.WeightsPerMeasurement,1);
                 
%             this.e_k = ((1-this.ProbOfDetection)*this.ClutterRate+sum(this.ProbOfDetection*g,1)).*this.e_k...
%                        ./(((1-this.ProbOfDetection)*this.ClutterRate+sum(this.ProbOfDetection*g,1)).*this.e_k ...
%                            + this.ClutterRate*(1-this.e_k));
            
            %e_i = p_i(ones(1,this.NumParticles),:)'.*g./(sum(g,2));
            
            % Resample (equivalent to Step 3 of [1]
            %this.Weights = this.Weights .* this.e_k;
            this.NumTargets = sum(this.Weights,2); % Compute total mass
            this.Particles = this.PredParticles;
            [this.Particles, this.Weights, idx] = ...
                this.Resampler.resample(this.Particles, ...
                    (this.Weights/this.NumTargets),this.NumParticles); % Resample
            this.Weights = this.Weights*this.NumTargets; % Rescale
%            this.e_k = this.e_k(:,idx);
        end
    
        function MeasLikelihood = get.MeasLikelihood(this)
            MeasLikelihood = getMeasLikelihood(this);
        end
    end
    
    methods (Access = protected)
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        function MeasLikelihood = getMeasLikelihood(this)
            if(isempty(this.pmeasLikelihood_))
                this.pmeasLikelihood_ = this.Model.Obs.pdf(this.MeasurementList, this.PredParticles);               
            end
            MeasLikelihood = this.pmeasLikelihood_;
        end
    end
end