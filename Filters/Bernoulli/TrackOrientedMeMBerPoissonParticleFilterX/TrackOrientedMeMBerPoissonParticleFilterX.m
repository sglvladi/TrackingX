classdef TrackOrientedMeMBerPoissonParticleFilterX < FilterX
% TrackOrientedMeMBerPoissonParticleFilterX class
%
% Summary of TrackOrientedMeMBerPoissonParticleFilterX:
% This is a class implementation of a Track-Oriented Marginal Hybrid 
% Multi-Bernoulli (MeMBer) Poisson Filter [1].
%
% TrackOrientedMeMBerPoissonParticleFilterX Properties:
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
% TrackOrientedMeMBerPoissonParticleFilterX Methods:
%   + TrackOrientedMeMBerPoissonParticleFilterX - Constructor method
%   + predict - Performs Bernoulli prediction step
%   + update - Performs Bernoulli update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.
%
% See also ParticleFilterX, KalmanFilerX.

    properties (Access = private, Hidden)
        
    end
    properties
        Bernoulli
        Poisson
        NumMeasurements
        NumBernoulliParticles = 1000
        NumPoissonParticles = 10000
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
        AssocLikelihoodMatrix
        MarginalAssociationProbabilities
        Hypothesiser = LoopyBeliefPropagationX();
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
    end
    
    methods
        function this = TrackOrientedMeMBerPoissonParticleFilterX_2(varargin)
        % TrackOrientedMeMBerPoissonParticleFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
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
        % * phd = TrackOrientedMeMBerPoissonParticleFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * phd = TrackOrientedMeMBerPoissonParticleFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * phd = TrackOrientedMeMBerPoissonParticleFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * phd = TrackOrientedMeMBerPoissonParticleFilterX(ssm,priorDistFcn) returns an object
        %   handle, preconfigured with the provided StateSpaceModel object 
        %   handle ssm and the prior information about the state, provided  
        %   in the form of the priorDistFcn function.
        % * phd = TrackOrientedMeMBerPoissonParticleFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
            
            % Call SuperClass method
            this@FilterX(varargin{:});
           
            this.Bernoulli.Predicted = {};
            this.Bernoulli.Updated = {};
            this.Poisson.Predicted = {};
            this.Poisson.Updated = {};
            
            if(nargin==0)
                this.Resampler = SystematicResamplerX();
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
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
                     if (isfield(config,'NumBernoulliParticles'))
                         this.NumBernoulliParticles = config.NumBernoulliParticles;
                     end
                     if (isfield(config,'NumPoissonParticles'))
                         this.NumPoissonParticles = config.NumPoissonParticles;
                     end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
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
             if (isfield(config,'NumBernoulliParticles'))
                 this.NumBernoulliParticles = config.NumBernoulliParticles;
             end
             if (isfield(config,'NumPoissonParticles'))
                 this.NumPoissonParticles = config.NumPoissonParticles;
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
        % * initialise(phd,ssm) initialises the TrackOrientedMeMBerPoissonParticleFilterX object 
        %   phd with the provided StateSpaceModelX object ssm.
        % * initialise(phd,priorParticles,priorWeights)initialises the 
        %   TrackOrientedMeMBerPoissonParticleFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(phd,ssm,priorDistFcn) initialises the TrackOrientedMeMBerPoissonParticleFilterX
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
                    if (isfield(config,'NumBernoulliParticles'))
                        this.NumBernoulliParticles = config.NumBernoulliParticles;
                    end
                    if (isfield(config,'NumPoissonParticles'))
                        this.NumPoissonParticles = config.NumPoissonParticles;
                    end
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
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
            if (isfield(config,'NumBernoulliParticles'))
                this.NumBernoulliParticles = config.NumBernoulliParticles;
            end
            if (isfield(config,'NumPoissonParticles'))
                this.NumPoissonParticles = config.NumPoissonParticles;
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
        % * TrackOrientedMeMBerPoissonParticleFilterX uses the Model class property, which should 
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

            % Predict existing Bernoulli tracks
            this.Bernoulli.Predicted = [];
            for i = 1:numel(this.Bernoulli.Updated)
                
                % Predict state particles
                this.Bernoulli.Predicted{i}.Particles = ...
                  this.Model.Dyn.feval(this.Bernoulli.Updated{i}.Particles, ...
                  this.Model.Dyn.random(size(this.Bernoulli.Updated{i}.Particles,2)));
              
                % Predict state weights
                this.Bernoulli.Predicted{i}.Weights = ...
                    this.ProbOfSurvive*this.Bernoulli.Updated{i}.Weights;
                
                % Predict existence probability
                this.Bernoulli.Predicted{i}.ProbOfExistence = ...
                    sum(this.Bernoulli.Predicted{i}.Weights)*this.Bernoulli.Updated{i}.ProbOfExistence;
                
                % Normalise state weights
                this.Bernoulli.Predicted{i}.Weights = ...
                    this.Bernoulli.Predicted{i}.Weights./sum(this.Bernoulli.Predicted{i}.Weights);
                
                % Copy forward the component trajectory
                this.Bernoulli.Predicted{i}.Trajectory = this.Bernoulli.Updated{i}.Trajectory;
            end

            % Predict existing PPP intensity
            this.Poisson.Predicted = struct;
            this.Poisson.Predicted.Particles = ...
                this.Model.Dyn.feval(this.Poisson.Updated.Particles, ...
                 this.Model.Dyn.random(size(this.Poisson.Updated.Particles,2)));
            this.Poisson.Predicted.Weights = this.ProbOfSurvive*this.Poisson.Updated.Weights;

            % Incorporate birth intensity into PPP
            if(this.BirthModel.Jk<1)
                Jk = ceil(this.BirthModel.Jk*size(this.Poisson.Predicted.Particles,2));
            else
                Jk = this.BirthModel.Jk;
            end
            
            % Generate birth particles
            [birthParticles, birthWeights] = this.BirthModel.BirthIntFcn(Jk);

            % Expand number of particles to accomodate for births
            this.Poisson.Predicted.Particles = [this.Poisson.Predicted.Particles, birthParticles]; 
            this.Poisson.Predicted.Weights = [this.Poisson.Predicted.Weights, birthWeights];
            this.Poisson.Predicted.NumTargets = sum(this.Poisson.Predicted.Weights);
            
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
            
            this.AssocLikelihoodMatrix = zeros(numel(this.Bernoulli.Predicted)+1, this.NumMeasurements+1);

            % Update existing tracks
            this.Bernoulli.Updated = [];
            for i = 1:numel(this.Bernoulli.Predicted)
                
                % Particles remain the same for now
                this.Bernoulli.Updated{i}.Particles = this.Bernoulli.Predicted{i}.Particles;

                % Compute measurement likelihoods for track
                g = this.Model.Obs.pdf(this.MeasurementList, this.Bernoulli.Updated{i}.Particles);
                
                % Store weights per hypothesis (Equivalent to creating multiple hypotheses)
                this.Bernoulli.Updated{i}.WeightsPerHypothesis = zeros(this.NumMeasurements+1,size(this.Bernoulli.Updated{i}.Particles,2));
                
                % Store existence per Hypothesis
                this.Bernoulli.Updated{i}.ExistencePerHypothesis = zeros(1,this.NumMeasurements+1);
                
                % Evaluate missed detection hypothesis
                this.Bernoulli.Updated{i}.WeightsPerHypothesis(1,:) =  (1-this.ProbOfDetection) * this.Bernoulli.Predicted{i}.Weights;
                this.AssocLikelihoodMatrix(i+1,1) = 1-this.Bernoulli.Predicted{i}.ProbOfExistence + this.Bernoulli.Predicted{i}.ProbOfExistence*sum(this.Bernoulli.Updated{i}.WeightsPerHypothesis(1,:));
                this.Bernoulli.Updated{i}.ExistencePerHypothesis(1) = ...
                    this.Bernoulli.Predicted{i}.ProbOfExistence*sum(this.Bernoulli.Updated{i}.WeightsPerHypothesis(1,:))...
                    /this.AssocLikelihoodMatrix(i+1,1);
                
                % Evaluate hypotheses with measurement updates
                this.Bernoulli.Updated{i}.WeightsPerHypothesis(2:end,:) = this.ProbOfDetection * g .* this.Bernoulli.Predicted{i}.Weights;
                
                this.Bernoulli.Updated{i}.ExistencePerHypothesis(2:end) = 1;
                this.AssocLikelihoodMatrix(i+1,2:end) = this.Bernoulli.Predicted{i}.ProbOfExistence*sum(this.Bernoulli.Updated{i}.WeightsPerHypothesis(2:end,:),2)';
                
                % Copy forward the component trajectory
                this.Bernoulli.Updated{i}.Trajectory = this.Bernoulli.Predicted{i}.Trajectory;
            end
            
            %disp(cellfun(@(c) c.ProbOfExistence, this.Bernoulli.Updated, 'UniformOutput', false));
            
            % Create a new track for each measurement by updating PPP with measurement
            % Compute g(z|x) matrix as in [1] 
            g = this.Model.Obs.pdf(this.MeasurementList, this.Poisson.Predicted.Particles);
                        
            % Compute C_k(z) Eq. (27) of [1]  
            Ck = this.ProbOfDetection*g.*this.Poisson.Predicted.Weights;
            
            % Calculate w^{n,i} Eq. (20) of [2]
            PoissonWeightsPerMeasurement = zeros(this.NumMeasurements, size(this.Poisson.Predicted.Particles,2));
            if(this.NumMeasurements>0)
                for j=1:this.NumMeasurements
                    PoissonWeightsPerMeasurement(j,:) = (this.ProbOfDetection*g(j,:).*this.Poisson.Predicted.Weights)/(this.ClutterIntFcn(this.MeasurementList(:,j))+sum(Ck(j),2));
                end
                    %PoissonWeightsPerMeasurement = (Ck+repmat(this.ClutterIntFcn(this.MeasurementList)',1,size(this.Poisson.Predicted.Weights,2))).*this.Poisson.Predicted.Weights;
            end
            
            % Extract new tracks from PPP density
            for j = 1:this.NumMeasurements
                newTrack = struct;
                
                % Get new particles and weights
                newTrack.Particles = this.Poisson.Predicted.Particles;
                newTrack.Weights = (PoissonWeightsPerMeasurement(j,:)./sum(PoissonWeightsPerMeasurement(j,:)));
                
                % Resample to correct number of particles
                [newTrack.Particles,newTrack.Weights,idx] = this.Resampler.resample(newTrack.Particles,newTrack.Weights,this.NumBernoulliParticles);
                
                Ck_j = Ck(j,idx);
                C = sum(Ck_j,2);
                
                % Compute weights per hypothesis
                this.AssocLikelihoodMatrix(1,j+1) = this.ClutterIntFcn(this.MeasurementList(:,j)) + C;
                
                % Compute probability of existence
                newTrack.ProbOfExistence = C/this.AssocLikelihoodMatrix(1,j+1);
                
                % Create new bernoulli component
                this.Bernoulli.Updated{end+1} = newTrack;
            end
            
            % Update (i.e., thin) PPP intensity
            this.Poisson.Updated = struct;
            this.Poisson.Updated.Particles = this.Poisson.Predicted.Particles;
            this.Poisson.Updated.Weights = (1-this.ProbOfDetection)*this.Poisson.Predicted.Weights;
            
            % Resample (equivalent to Step 3 of [1]
            NumPPPTargets = sum(this.Poisson.Updated.Weights,2);
            [this.Poisson.Updated.Particles,this.Poisson.Updated.Weights] = ...
                this.Resampler.resample(this.Poisson.Updated.Particles, ...
                    (this.Poisson.Updated.Weights/NumPPPTargets),this.NumPoissonParticles); % Resample
            this.Poisson.Updated.Weights = this.Poisson.Updated.Weights*NumPPPTargets;
            
            
            % Compute marginal association probabilities for known targets
            [a,b] = lbp(this.AssocLikelihoodMatrix(2:end,:),this.AssocLikelihoodMatrix(1,2:end));
            this.MarginalAssociationProbabilities = [[0,b'];a];
            %this.MarginalAssociationProbabilities = this.Hypothesiser.hypothesise(this.AssocLikelihoodMatrix);
            
            % TOMB/P algorithm for forming new tracks
            % Form continuing tracks
            for i = 1 : numel(this.Bernoulli.Predicted)
                % Update Existence probability
                this.Bernoulli.Updated{i}.ProbOfExistence = sum(this.MarginalAssociationProbabilities(i+1,:).*this.Bernoulli.Updated{i}.ExistencePerHypothesis);
                
                % Update weights
                this.Bernoulli.Updated{i}.Weights = sum(this.MarginalAssociationProbabilities(i+1,:)'.*this.Bernoulli.Updated{i}.WeightsPerHypothesis,1);
                this.Bernoulli.Updated{i}.Weights = this.Bernoulli.Updated{i}.Weights./sum(this.Bernoulli.Updated{i}.Weights,2);
                % Resample
                [this.Bernoulli.Updated{i}.Particles,this.Bernoulli.Updated{i}.Weights] = ...
                    this.Resampler.resample(this.Bernoulli.Updated{i}.Particles, this.Bernoulli.Updated{i}.Weights,this.NumBernoulliParticles); 
                
                this.Bernoulli.Updated{i}.Trajectory.Mean(:,end+1) = mean(this.Bernoulli.Updated{i}.Particles,2);
                this.Bernoulli.Updated{i}.Trajectory.Covariance(:,:,end+1) = cov(this.Bernoulli.Updated{i}.Particles');
                this.Bernoulli.Updated{i}.Trajectory.ProbOfExistence(:,end+1) = this.Bernoulli.Updated{i}.ProbOfExistence;
            end
            % Form new tracks (already single hypothesis)
            for j = numel(this.Bernoulli.Predicted)+ 1 : numel(this.Bernoulli.Updated)
                this.Bernoulli.Updated{j}.ProbOfExistence = this.Bernoulli.Updated{j}.ProbOfExistence*this.MarginalAssociationProbabilities(1,j+1-(numel(this.Bernoulli.Predicted)));
                
                this.Bernoulli.Updated{j}.Trajectory.Mean = mean(this.Bernoulli.Updated{j}.Particles,2);
                this.Bernoulli.Updated{j}.Trajectory.Covariance = cov(this.Bernoulli.Updated{j}.Particles');
                this.Bernoulli.Updated{j}.Trajectory.ProbOfExistence = this.Bernoulli.Updated{j}.ProbOfExistence;
            end
            
            disp(cellfun(@(c) c.ProbOfExistence, this.Bernoulli.Updated, 'UniformOutput', false));
                        
            % Prune bernoulli components
            this.Bernoulli.Updated = this.Bernoulli.Updated(find(cell2mat(cellfun(@(c) c.ProbOfExistence, this.Bernoulli.Updated, 'UniformOutput', false))>5*10^(-2)));
            
            %disp(cellfun(@(c) c.ProbOfExistence, this.Bernoulli.Updated, 'UniformOutput', false));
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