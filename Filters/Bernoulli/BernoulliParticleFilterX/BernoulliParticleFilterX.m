classdef BernoulliParticleFilterX < FilterX
% BernoulliParticleFilterX class
%
% Summary of BernoulliParticleFilterX:
% This is a class implementation of a Bernoulli Particle Filter [1].
%
% BernoulliParticleFilterX Properties:
%   + NumParticles - The number of particles employed by the PHD Filter
%   + Particles - A (NumStateDims x NumParticles) matrix used to store 
%                       the last computed/set filtered particles  
%   + Weights - A (1 x NumParticles) vector used to store the weights
%               of the last computed/set filtered particles
%   + PredParticles - A (NumStateDims x NumParticles) matrix used to store 
%                     the last computed/set predicted particles  
%   + PredWeights - A (1 x NumParticles) vector used to store the weights
%                   of the last computed/set predicted particles
%   + MeasurementList - A (NumObsDims x NumObs) matrix used to store the received 
%                       measurements
%   + ResamplingScheme - Method used for particle resampling, specified as 
%                        'multinomial', 'systematic'. Default = 'systematic'
%   + ResamplingPolicy - A (1 x 2) cell array, specifying the resampling trigger
%                        conditions. ReamplingPolicy{1} should be a string
%                        which can be either "TimeInterval", in which case 
%                        ReamplingPolicy{2} should specify the number of 
%                        iterations after which resampling should be done,
%                        or "EffectiveRatio", in which case Resampling{2}
%                        should specify the minimum ratio of effective particles
%                        which when reached will trigger the resampling process.
%                        Default ResamplingPolicy = {"TimeInterval",1}, meaning
%                        that resampling is performed on every iteration of
%                        the Particle Filter (upon update).                       
%   + Resampler - An object handle to a ResamplerX subclass. If a 
%                 Resampler is provided, then it will override any choice
%                 specified within the ResamplingScheme. ResamplingPolicy
%                 will not be affected.
%   + BirthModel - A (1 x 3) cell array, specifying the particle birth
%                   scheme. BirthModel{1} should be a string which can be
%                   set to either "Mixture", in which case BirthModel{2}
%                   should specify the probability of birth of new particles,
%                   or "Expansion", in which case BirthModel{2} should
%                   specify the number of particles to be "birthed" at
%                   each iteration of the filter.
%                   Default BirthModel = {"Mixture",0.5} meaning that
%                   particles are birthed using the mixture scheme, with
%                   a birth probability of 50%.
%   + BirthIntFcn - A function handle, which when called generates a set 
%                   of initial particles and weights.
%   + ProbOfSurvive - The probability that a target may cease to exist
%                   between consecutive iterations of the filter.
%   + ProbOfDetection - The probablity that a target will be detected in
%                       a given measurement scan.
%   + NumTargets - The estimated number of targets following an update step.
%   + Model - An object handle to StateSpaceModelX object
%       + Dyn - Object handle to DynamicModelX SubClass      
%       + Obs - Object handle to ObservationModelX SubClass 
%       + Ctr - Object handle to ControlModelX SubClass 
%
% BernoulliParticleFilterX Methods:
%   + BernoulliParticleFilterX - Constructor method
%   + predict - Performs Bernoulli prediction step
%   + update - Performs Bernoulli update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1]  Risti, Branko & Vo, Ba-Ngu & Vo, Ba-Ngu & Farina, Alfonso. (2013). A Tutorial on Bernoulli Filters: Theory, Implementation and Applications. IEEE Transactions on Signal Processing. 61. 10.1109/TSP.2013.2257765. 
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
        ProbOfExistence
        BirthParticles = []
        BirthWeights = []
        PredParticles   
        PredWeights
        PredProbOfExistence
        MeasurementList
        ResamplingScheme = 'Systematic'
        ResamplingPolicy = {'TimeInterval',1}                    
        Resampler        
        BirthModel
        ProbOfBirth
        ProbOfSurvive  
        ProbOfDetection 
        ClutterRate
        ClutterIntFcn
        MeasLikelihood
        MeasWeights = 1;
        WeightsPerMeasurement
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
    end
    
    methods
        function this = BernoulliParticleFilterX(varargin)
        % BernoulliParticleFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % NumParticles: scalar
        %   The number of particles to be employed by the  SMC PHD Filter. 
        %   (default = 10000)
        % PriorParticles: (NumStateDims x NumParticles) matrix
        %   The initial set of particles to be used by the Particle Filter. 
        %   These are copied into the Particles property by the constructor.
        % PriorWeights: (1 x NumParticles) row vector
        %   The initial set of weights to be used by the Particle Filter. 
        %   These are copied into the Weights property by the constructor. 
        %   (default = 1/NumParticles, which implies uniform weights)
        % PriorDistFcn: function handle
        %   A function handle, of the form [parts, weights] = PriorDistFcn(NumParticles),
        %   which when called generates a set of initial particles and weights, 
        %   that are consecutively copied into the Particles and Weights properties
        %   respectively. The function should accept exactly ONE argument, 
        %   which is the number of particles (NumParticles) to be generated and
        %   return 2 outputs. If a PriorDistFcn is specified, then any values provided for the
        %   PriorParticles and PriorWeights arguments are ignored.
        % ResamplingScheme: string 
        %   Method used for particle resampling, specified as either 'Multinomial'
        %   or 'Systematic'. (default = 'Systematic')
        % ResamplingPolicy: (1 x 2) cell array
        %   Specifies the resampling trigger conditions. ReamplingPolicy{1} 
        %   should be a string which can be either:
        %       1) "TimeInterval", in which case ReamplingPolicy{2} should 
        %          be a scalar specifying the number of iterations after which
        %          resampling should be performed, or
        %       2) "EffectiveRatio", in which case Resampling{2} should be
        %          a scalar specifying the minimum ratio of effective particles 
        %          which, when reached, will trigger the resampling process
        %   (default ResamplingPolicy = {'TimeInterval',1}], implying that 
        %   resampling is performed on every update of the Particle Filter).                       
        % Resampler: ResamplerX object handle, optional
        %   An object handle to a ResamplerX subclass. If a Resampler is provided,
        %   then it will override any choice specified within the ResamplingScheme. 
        %   ResamplingPolicy will not be affected.
        % BirthModel: (1 x 2) cell array
        %   Specifies the particle birth scheme. BirthModel{1} should be 
        %   a string which can be set to either:
        %       1) "Mixture", in which case BirthModel{2} should specify the
        %          probability of birth of new particles, or
        %       2) "Expansion", in which case BirthModel{2} should specify 
        %          the number of particles to be "birthed" at each iteration 
        %          of the filter.
        %   (default BirthModel = {"Mixture",0.5} implying that particles 
        %    are born using the mixture scheme, with a birth probability of 50%)
        % BirthIntFcn: function handle
        %   A function handle, [parts, weights] = BirthIntFcn(NumParticles), 
        %   which when called generates a set of initial particles and weights.
        % ProbOfSurvive: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % ProbOfDetection: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        % NumTargets: scalar
        %   The estimated number of targets following an update step.
        %
        % Usage
        % -----
        % * phd = BernoulliParticleFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * phd = BernoulliParticleFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * phd = BernoulliParticleFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * phd = BernoulliParticleFilterX(ssm,priorDistFcn) returns an object
        %   handle, preconfigured with the provided StateSpaceModel object 
        %   handle ssm and the prior information about the state, provided  
        %   in the form of the priorDistFcn function.
        % * phd = BernoulliParticleFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
            
            % Call SuperClass method
            this@FilterX(varargin{:});
           
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
                     if (isfield(config,'BirthModel'))
                         this.BirthModel = config.BirthModel;
                     end
                     if (isfield(config,'ProbOfSurvive'))
                         this.ProbOfSurvive = config.ProbOfSurvive;
                     end
                     if (isfield(config,'ProbOfDetection'))
                         this.ProbOfDetection = config.ProbOfDetection;
                     end
                     if (isfield(config,'ClutterRate'))
                         this.ClutterRate = config.ClutterRate;
                     end
                     if (isfield(config,'ClutterIntFcn'))
                         this.ClutterIntFcn = config.ClutterIntFcn;
                     end
                     if (isfield(config,'ProbOfExistence'))
                         this.ProbOfExistence = config.ProbOfExistence;
                     end
                end
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
             if (isfield(config,'BirthModel'))
                 this.BirthModel = config.BirthModel;
             end
             if (isfield(config,'ProbOfSurvive'))
                 this.ProbOfSurvive = config.ProbOfSurvive;
             end
             if (isfield(config,'ProbOfDetection'))
                 this.ProbOfDetection = config.ProbOfDetection;
             end
             if (isfield(config,'ClutterRate'))
                 this.ClutterRate = config.ClutterRate;
             end
             if (isfield(config,'ClutterIntFcn'))
                 this.ClutterIntFcn = config.ClutterIntFcn;
             end
             if (isfield(config,'ProbOfExistence'))
                 this.ProbOfExistence = config.ProbOfExistence;
             end
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the Bernoulli Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % NumParticles: scalar
        %   The number of particles to be employed by the  SMC PHD Filter. 
        %   (default = 10000)
        % PriorParticles: (NumStateDims x NumParticles) matrix
        %   The initial set of particles to be used by the Particle Filter. 
        %   These are copied into the Particles property by the constructor.
        % PriorWeights: (1 x NumParticles) row vector
        %   The initial set of weights to be used by the Particle Filter. 
        %   These are copied into the Weights property by the constructor. 
        %   (default = 1/NumParticles, which implies uniform weights)
        % PriorDistFcn: function handle
        %   A function handle, of the form [parts, weights] = PriorDistFcn(NumParticles),
        %   which when called generates a set of initial particles and weights, 
        %   that are consecutively copied into the Particles and Weights properties
        %   respectively. The function should accept exactly ONE argument, 
        %   which is the number of particles (NumParticles) to be generated and
        %   return 2 outputs. If a PriorDistFcn is specified, then any values provided for the
        %   PriorParticles and PriorWeights arguments are ignored.
        % ResamplingScheme: string 
        %   Method used for particle resampling, specified as either 'Multinomial'
        %   or 'Systematic'. (default = 'Systematic')
        % ResamplingPolicy: (1 x 2) cell array
        %   Specifies the resampling trigger conditions. ReamplingPolicy{1} 
        %   should be a string which can be either:
        %       1) "TimeInterval", in which case ReamplingPolicy{2} should 
        %          be a scalar specifying the number of iterations after which
        %          resampling should be performed, or
        %       2) "EffectiveRatio", in which case Resampling{2} should be
        %          a scalar specifying the minimum ratio of effective particles 
        %          which, when reached, will trigger the resampling process
        %   (default ResamplingPolicy = {'TimeInterval',1}], implying that 
        %   resampling is performed on every update of the Particle Filter).                       
        % Resampler: ResamplerX object handle, optional
        %   An object handle to a ResamplerX subclass. If a Resampler is provided,
        %   then it will override any choice specified within the ResamplingScheme. 
        %   ResamplingPolicy will not be affected.
        % BirthModel: (1 x 2) cell array
        %   Specifies the particle birth scheme. BirthModel{1} should be 
        %   a string which can be set to either:
        %       1) "Mixture", in which case BirthModel{2} should specify the
        %          probability of birth of new particles, or
        %       2) "Expansion", in which case BirthModel{2} should specify 
        %          the number of particles to be "birthed" at each iteration 
        %          of the filter.
        %   (default BirthModel = {"Mixture",0.5} implying that particles 
        %    are born using the mixture scheme, with a birth probability of 50%)
        % BirthIntFcn: function handle
        %   A function handle, [parts, weights] = BirthIntFcn(NumParticles), 
        %   which when called generates a set of initial particles and weights.
        % ProbOfSurvive: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % ProbOfDetection: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        % NumTargets: scalar
        %   The estimated number of targets following an update step.
        %
        % Usage
        % ----- 
        % * initialise(phd,ssm) initialises the BernoulliParticleFilterX object 
        %   phd with the provided StateSpaceModelX object ssm.
        % * initialise(phd,priorParticles,priorWeights)initialises the 
        %   BernoulliParticleFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(phd,ssm,priorDistFcn) initialises the BernoulliParticleFilterX
        %   object pf with the provided StateSpaceModel object handle ssm
        %   and the prior information about the state, provided in the form 
        %   of the priorDistFcn function.
        % * initialise(phd,___,Name,Value,___) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value
        %   pair arguments.
        %
        %  See also predict, update.   
                        
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
                    if (isfield(config,'BirthModel'))
                        this.BirthModel = config.BirthModel;
                    end
                    if (isfield(config,'ProbOfSurvive'))
                        this.ProbOfSurvive = config.ProbOfSurvive;
                    end
                    if (isfield(config,'ProbOfDetection'))
                        this.ProbOfDetection = config.ProbOfDetection;
                    end
                    if (isfield(config,'ClutterRate'))
                         this.ClutterRate = config.ClutterRate;
                    end
                    if (isfield(config,'ClutterIntFcn'))
                         this.ClutterIntFcn = config.ClutterIntFcn;
                    end
                    if (isfield(config,'ProbOfExistence'))
                         this.ProbOfExistence = config.ProbOfExistence;
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
             if (isfield(config,'BirthModel'))
                 this.BirthModel = config.BirthModel;
             end
             if (isfield(config,'ProbOfSurvive'))
                 this.ProbOfSurvive = config.ProbOfSurvive;
             end
             if (isfield(config,'ProbOfDetection'))
                 this.ProbOfDetection = config.ProbOfDetection;
             end
             if (isfield(config,'ProbOfExistence'))
                 this.ProbOfExistence = config.ProbOfExistence;
             end
        end
        
        function predict(this)
        % PREDICT Perform Bernoulli Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted Bernoulli PDF
        %
        % More details
        % ------------
        % * BernoulliParticleFilterX uses the Model class property, which should 
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
            
            % Compute predicted probability of Existence - (28) of [1]
            this.PredProbOfExistence = this.BirthModel.ProbOfBirth*(1-this.ProbOfExistence)...
                                        + this.ProbOfSurvive*this.ProbOfExistence;
            
            % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.BirthModel.Schema, 'Expansion')) 
                
                this.NumParticles = size(this.Particles,2);
                if(this.BirthModel.Jk<1)
                    Jk = ceil(this.BirthModel.Jk*this.NumParticles);
                else
                    Jk = this.BirthModel.Jk;
                end
                
                % Expand number of particles to accomodate for births
                this.PredParticles = [this.Particles, zeros(size(this.Particles, 1), Jk)]; 
                this.PredWeights = [this.Weights, zeros(1, Jk)];
                this.NumParticlesTotal = this.NumParticles + Jk;  
                
                % Generate birth particles
                [birthParticles, birthWeights] = this.BirthModel.BirthIntFcn(Jk);
                this.PredParticles(:,this.NumParticles+1:end) = birthParticles;
                
                % Predict all particles - (81) of [1]
                this.PredParticles = this.Model.Dyn.feval(this.PredParticles, this.Model.Dyn.random(this.NumParticlesTotal)); 
                
                % Predict weights - (82) of [1]
                this.PredWeights(:,1:this.NumParticles) = ...
                    this.ProbOfSurvive*(this.ProbOfExistence/this.PredProbOfExistence)*this.Weights(:,1:this.NumParticles);
                this.PredWeights(:,this.NumParticles+1:end) = ...
                    this.BirthModel.ProbOfBirth*(1-this.ProbOfExistence)*birthWeights/(this.PredProbOfExistence*Jk);
                
            else
                % Propagate old (k-1) particles and generate new birth particles
                this.NumParticles = size(this.Particles,2);

                % Copy over particles from last timestep
                this.PredParticles = this.Particles; 
                this.PredWeights = this.Weights;  

                % Append any birth particles
                if(size(this.BirthParticles,2)>0)
                    this.PredParticles = [this.PredParticles, this.BirthParticles];
                    this.PredWeights = [this.PredWeights, this.BirthWeights]; 
                end

                this.NumParticlesTotal = size(this.PredParticles,2); 

                % Predict all particles - (81) of [1]
                this.PredParticles = this.Model.Dyn.feval(this.PredParticles, this.Model.Dyn.random(this.NumParticlesTotal)); 

                % Predict weights - (82) of [1]
                this.PredWeights(:,1:this.NumParticles) = ...
                    this.ProbOfSurvive*(this.ProbOfExistence/this.PredProbOfExistence)*this.Weights(:,1:this.NumParticles);

                if(size(this.BirthParticles,2)>0)
                    this.PredWeights(:,this.NumParticles+1:end) = ...
                        this.BirthModel.ProbOfBirth*(1-this.ProbOfExistence)*this.BirthWeights/(this.PredProbOfExistence);
                end
            end
        end
        
        function update(this)
        % UPDATE Perform Bernoulli Filter update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
        
            this.NumMeasurements = size(this.MeasurementList,2);
            
            % Compute g(z|x) matrix as in [1] 
            g = this.MeasLikelihood;

             % Approximate integral I_1 - (84) of [1]
            I_1 = sum(this.ProbOfDetection*this.PredWeights);
            
            % Approximate integral I_2(z) - (85) of [1]
            I_2 = sum(this.ProbOfDetection*g.*this.PredWeights,2);
            
            % Compute ?k approximation - (86) of [1]
            D_k = I_1 - sum(I_2./(this.ClutterRate.*this.ClutterIntFcn(this.MeasurementList))');
            
            % Update existence probability - (56) of [1] 
            this.ProbOfExistence = (1-D_k)*this.PredProbOfExistence/(1-this.PredProbOfExistence*D_k);
            
            % Update weights - (87) of [1]
            this.Weights = (1-this.ProbOfDetection + this.ProbOfDetection*sum(g./(this.ClutterRate.*this.ClutterIntFcn(this.MeasurementList))')).*this.PredWeights;
            
            % Normalise weights 
            this.Weights = this.Weights./sum(this.Weights);
            
            % Resample
            this.Particles = this.PredParticles;
            [this.Particles, this.Weights, idx] = ...
                this.Resampler.resample(this.Particles, ...
                    this.Weights,this.NumParticles);
                
            if(~strcmp(this.BirthModel.Schema, 'Expansion'))
                % Generate birth particles
                [this.BirthParticles, this.BirthWeights] = this.BirthModel.BirthIntFcn(this.MeasurementList);
            end
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