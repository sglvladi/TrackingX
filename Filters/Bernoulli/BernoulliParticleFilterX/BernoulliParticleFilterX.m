classdef BernoulliParticleFilterX < ParticleFilterX
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
%   + SurvivalProbability - The probability that a target may cease to exist
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

    properties 
        BirthScheme = {'Mixture',0.5}; 
        SurvivalProbability
        MeasWeights = 1;
        
        weightsPerHypothesis_
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
        end
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
        %   The number of particles that should be employed by the filter
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a ParticleStateX instance, then it will be converted into
        %   one by sampling from the provided State's distribution.
        % ResamplingScheme: string , optional
        %   Method used for particle resampling, specified as either 'Multinomial'
        %   or 'Systematic'. (default = 'Systematic')
        % ResamplingPolicy: (1 x 2) cell array, optional
        %   Specifies the resampling trigger conditions. ReamplingPolicy{1} 
        %   should be a string which can be either:
        %       1) "TimeInterval", in which case ReamplingPolicy{2} should 
        %          be a scalar specifying the number of iterations after which
        %          resampling should be performed, or
        %       2) "EffectiveRatio", in which case Resampling{2} should be
        %          a scalar specifying the minimum ratio of effective particles 
        %          which, when reached, will trigger the resampling process
        %  (default ResamplingPolicy = {'TimeInterval',1}], implying that 
        %   resampling is performed on every update of the Particle Filter).                       
        % Resampler: ResamplerX object handle, optional
        %   An object handle to a ResamplerX subclass. If a Resampler is provided,
        %   then it will override any choice specified within the ResamplingScheme. 
        %   ResamplingPolicy will not be affected.
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
        % INITIALISE Initialise the Bernoulli Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a ParticleStateX instance, then it will be converted into
        %   one by sampling from the provided State's distribution.
        % ResamplingScheme: string , optional
        %   Method used for particle resampling, specified as either 'Multinomial'
        %   or 'Systematic'. (default = 'Systematic')
        % ResamplingPolicy: (1 x 2) cell array, optional
        %   Specifies the resampling trigger conditions. ReamplingPolicy{1} 
        %   should be a string which can be either:
        %       1) "TimeInterval", in which case ReamplingPolicy{2} should 
        %          be a scalar specifying the number of iterations after which
        %          resampling should be performed, or
        %       2) "EffectiveRatio", in which case Resampling{2} should be
        %          a scalar specifying the minimum ratio of effective particles 
        %          which, when reached, will trigger the resampling process
        %  (default ResamplingPolicy = {'TimeInterval',1}], implying that 
        %   resampling is performed on every update of the Particle Filter).                       
        % Resampler: ResamplerX object handle, optional
        %   An object handle to a ResamplerX subclass. If a Resampler is provided,
        %   then it will override any choice specified within the ResamplingScheme. 
        %   ResamplingPolicy will not be affected.
        %
        %  See also predict, update, smooth. 
                    
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
            
            % Compute predicted probability of Existence - (28) of [1]
            priorExistence = this.StatePosterior.Metadata.ExistenceProbability;
            predExistence = this.Model.Birth.BirthProbability*(1-priorExistence)...
                            + this.SurvivalProbability*priorExistence;
            
            % Propagate old (k-1) particles and generate new birth particles
            if(strcmp(this.BirthScheme(1), 'Expansion'))
                
                % Expansion method is equivalent to Eqs. (25-26) of [1]
                
                % Number of birth particles
                Jk = this.BirthScheme{2};
                
                % Predict all particles - (81) of [1]
                predParticles = this.Model.Transition.feval(this.StatePosterior.Particles, true);  
                
                % Generate birth particles 
                [birthParticles, birthWeights] = this.Model.Birth.random(Jk);
                
                % Predict weights - (82) of [1]
                predWeights = this.SurvivalProbability*(priorExistence/predExistence)*this.StatePosterior.Weights;
                birthWeights = this.Model.Birth.BirthProbability*(1-priorExistence)*birthWeights/(predExistence*Jk);
                
                this.StatePrediction = ParticleStateX([predParticles,birthParticles], [predWeights, birthWeights]);
                this.StatePrediction.Metadata.ExistenceProbability = predExistence;
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
        
            numMeasurements = this.MeasurementList.NumMeasurements;
            predExistence = this.StatePrediction.Metadata.ExistenceProbability;
            
            % Compute g(z|x) matrix as in [1] 
            g = this.MeasurementLikelihoodsPerParticle;

            P_Ds = this.Model.Detection.pdf(this.StatePrediction.Particles);

            % Approximate integral I_1 - (84) of [1]
            I_1 = sum(P_Ds.*this.StatePrediction.Weights);
                
            if(numMeasurements>0)
                
                % Approximate integral I_2(z) - (85) of [1]
                I_2 = sum(P_Ds.*g.*this.StatePrediction.Weights,2);

                % Compute Dk approximation - (86) of [1]
                D_k = I_1 - sum(I_2./(this.Model.Clutter.ClutterRate.*this.Model.Clutter.pdf(this.MeasurementList.Vectors,'spatial'))');
            
                % Update weights - (87) of [1]
                weights = (1-P_Ds + P_Ds.*sum(g./(this.Model.Clutter.ClutterRate.*this.Model.Clutter.pdf(this.MeasurementList.Vectors,'spatial'))',1)).*this.StatePrediction.Weights;
            
            else
                
                % Compute Dk approximation - (86) of [1]
                D_k = I_1;
                
                % Update weights - (87) of [1]
                weights = (1-P_Ds).*this.StatePrediction.Weights;
            end
            
            % Update existence probability - (56) of [1]
            postExistence = (1-D_k)*predExistence/(1-predExistence*D_k);
            
            % Normalise weights 
            weights = weights./sum(weights);
            
            % Resample
            particles = this.StatePrediction.Particles;
            [particles, weights, idx] = ...
                this.Resampler.resample(particles, weights, this.NumParticles);
            
            % Store posterior
            this.StatePosterior = ParticleStateX(particles,weights);
            this.StatePosterior.Metadata.ExistenceProbability = postExistence;
            
%             if(~strcmp(this.BirthModel.Schema, 'Expansion'))
%                 % Generate birth particles
%                 [this.BirthParticles, this.BirthWeights] = this.BirthModel.BirthIntFcn(this.MeasurementList);
%             end
        end
    end
end