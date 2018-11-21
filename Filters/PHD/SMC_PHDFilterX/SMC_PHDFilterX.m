classdef SMC_PHDFilterX < ParticleFilterX
% SMC_PHDFilterX class
%
% Summary of SMC_PHDFilterX:
% This is a class implementation of a Sequential Monte Carlo (SMC) Probabilistic
% Hypothesis Density (PHD) Filter.
%
% SMC_PHDFilterX Properties:
%   + NumParticles - The number of particles employed by the PHD Filter
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (NumMeasDims x NumMeasurements) matrix used to store
%                       the received measurements
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
%   + BirthScheme - A (1 x 3) cell array, specifying the particle birth
%                   scheme. BirthScheme{1} should be a string which can be
%                   set to either "Mixture", in which case BirthScheme{2}
%                   should specify the probability of birth of new particles,
%                   or "Expansion", in which case BirthScheme{2} should
%                   specify the number of particles to be "birthed" at
%                   each iteration of the filter.
%                   Default BirthScheme = {"Mixture",0.5} meaning that
%                   particles are birthed using the mixture scheme, with
%                   a birth probability of 50%.
%   + BirthIntFcn - A function handle, which when called generates a set 
%                   of initial particles and weights.
%   + ProbOfDeath - The probability that a target may cease to exist
%                   between consecutive iterations of the filter.
%   + DetectionProbability - The probablity that a target will be detected in
%                       a given measurement scan.
%   + NumTargets - The estimated number of targets following an update step.
%   + Model - An object handle to StateSpaceModelX object
%       + Transition (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass
%       + Clutter
%       + Birth
%
% SMC_PHDFilterX Methods:
%   + SMC_PHDFilterX - Constructor method
%   + predict - Performs SMC_PHD prediction step
%   + update - Performs SMC_PHD update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1]  B. N. Vo, S. Singh and A. Doucet, "Sequential Monte Carlo methods for multitarget filtering with random finite sets," in IEEE Transactions on Aerospace and Electronic Systems, vol. 41, no. 4, pp. 1224-1245, Oct. 2005.
% [2]  P. Horridge and S. Maskell,  “Using a probabilistic hypothesis density filter to confirm tracks in a multi-target environment,” in2011 Jahrestagung der Gesellschaft fr Informatik, October 2011.
% [3]  B. ngu Vo and S. Singh, “Sequential monte carlo implementation of the phd filter for multi-targettracking,” inIn Proceedings of the Sixth International Conference on Information Fusion, pp. 792–799, 2003.
%
% See also ParticleFilterX, KalmanFilerX.
    
    properties       
        BirthScheme = {'Mixture',0.5};
        SurvivalProbability  
        DetectionProbability 
        MeasWeights = 1;
        
        
        weightsPerHypothesis_
    end
    
    properties (Dependent)
        IntensityPerHypothesis
    end
    
    properties (Access=protected)
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            initialise_@ParticleFilterX(this,config);
            if (isfield(config,'BirthScheme'))
                this.BirthScheme = config.BirthScheme;
            end
            if (isfield(config,'SurvivalProbability'))
                this.SurvivalProbability = config.SurvivalProbability;
            end
            if (isfield(config,'DetectionProbability'))
                this.DetectionProbability = config.DetectionProbability;
            end
        end
    end
    
    methods
        function this = SMC_PHDFilterX(varargin)
        % SMC_PHDFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % NumParticles: scalar
        %   The number of particles that should be employed by the filter
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a ParticleStateX instance, then it will be converted into
        %   one by sampling from the provided State's distribution.
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
        % BirthScheme: (1 x 2) cell array
        %   Specifies the particle birth scheme. BirthScheme{1} should be 
        %   a string which can be set to either:
        %       1) "Mixture", in which case BirthScheme{2} should specify the
        %          probability of birth of new particles, or
        %       2) "Expansion", in which case BirthScheme{2} should specify 
        %          the number of particles to be "birthed" at each iteration 
        %          of the filter.
        %   (default BirthScheme = {"Mixture",0.5} implying that particles 
        %    are born using the mixture scheme, with a birth probability of 50%)
        % BirthIntFcn: function handle
        %   A function handle, [parts, weights] = BirthIntFcn(NumParticles), 
        %   which when called generates a set of initial particles and weights.
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % DetectionProbability: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        % NumTargets: scalar
        %   The estimated number of targets following an update step.
        %
        % Usage
        % -----
        % * phd = SMC_PHDFilterX(___,Name,Value) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
            
            % Call SuperClass method
            this@ParticleFilterX(varargin{:});
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the SMC PHD Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------
         % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % NumParticles: scalar
        %   The number of particles that should be employed by the filter
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a ParticleStateX instance, then it will be converted into
        %   one by sampling from the provided State's distribution.
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
        % BirthScheme: (1 x 2) cell array
        %   Specifies the particle birth scheme. BirthScheme{1} should be 
        %   a string which can be set to either:
        %       1) "Mixture", in which case BirthScheme{2} should specify the
        %          probability of birth of new particles, or
        %       2) "Expansion", in which case BirthScheme{2} should specify 
        %          the number of particles to be "birthed" at each iteration 
        %          of the filter.
        %   (default BirthScheme = {"Mixture",0.5} implying that particles 
        %    are born using the mixture scheme, with a birth probability of 50%)
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % DetectionProbability: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        %
        % Usage
        % ----- 
        % * initialise(phd,___,Name,Value) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value
        %   pair arguments.
        %
        %  See also predict, update.   
                        
            initialise@ParticleFilterX(this);
            
            % First check to see if a structure was received
            if(nargin==2)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function predict(this)
        % PREDICT Perform SMC PHD Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted PHD
        %
        % More details
        % ------------
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
                    
            % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.BirthScheme(1), 'Expansion')) 
                % Expansion method is equivalent to Eqs. (25-26) of [1]
                
                % Number of birth particles
                Jk = this.BirthScheme{2};

                % Generate Np normally predicted particles
                predParticles = this.Model.Transition.feval(this.StatePosterior.Particles, true); 
                predWeights = this.SurvivalProbability *this.StatePosterior.Weights;

                % Generate birth particles 
                [birthParticles, birthWeights] = this.Model.Birth.random(Jk);
                
                this.StatePrediction = ParticleStateX([predParticles,birthParticles], [predWeights, birthWeights]);
            end
        end
        
        function idx = update(this)
        % UPDATE Perform SMC PHD update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
        
            numMeasurements = size(this.MeasurementList,2);
            
            % Compute g(z|x) matrix as in [1] 
            g = this.MeasurementLikelihoodsPerParticle;
                                    
            % Calculate w^{n,i} Eq. (20) of [2]
            this.weightsPerHypothesis_ = [(1-this.DetectionProbability).*this.StatePrediction.Weights;
                                          zeros(numMeasurements, this.StatePrediction.NumParticles)];
            if(numMeasurements>0)
                % Compute C_k(z) Eq. (27) of [1] 
                Ck = sum(this.DetectionProbability*g.*this.StatePrediction.Weights,2)';
                this.weightsPerHypothesis_(2:end,:) = (this.DetectionProbability*g./(this.Model.Clutter.pdf(this.MeasurementList)+Ck)');
            end
            this.weightsPerHypothesis_(2:end,:) = this.MeasWeights'.*this.weightsPerHypothesis_(2:end,:).*this.StatePrediction.Weights;
            
            % Update weights Eq. (28) of [1]
            postWeights = sum(this.weightsPerHypothesis_,1);
            
            % Resample (equivalent to Step 3 of [1]
            numTargets = sum(postWeights,2); % Compute total mass
            postParticles = this.StatePrediction.Particles;
            [postParticles, postWeights, idx] = ...
                this.Resampler.resample(postParticles, ...
                    (postWeights/numTargets),this.NumParticles); % Resample
            postWeights = postWeights*numTargets; % Rescale
            this.StatePosterior = ParticleStateX(postParticles,postWeights);
        end
        
        function intensityPerHypothesis = get.IntensityPerHypothesis(this)
            intensityPerHypothesis = sum(this.weightsPerHypothesis_,2)';
        end
    
    end
 
end