classdef ParticleFilterX < FilterX
% ParticleFilterX class
%
% Summary of ParticleFilterX:
% This is a class implementation of a SIR Particle Filter. (Alg. 4 of [1])
%
% ParticleFilterX Properties: (*)
%   + NumParticles - The number of particles employed by the Particle Filter
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (yDim x 1) matrix used to store the received measurement
%   + ControlInput - A (NumCtrDims x 1) matrix used to store the last received 
%                    control input
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
%   + Model - An object handle to StateSpaceModelX object
%       + Transition (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass   
%
%   (*) NumStateDims, NumObsDims and NumCtrDims denote the dimentionality of 
%       the state, measurement and control vectors respectively.
%
% ParticleFilterX Methods:
%   + ParticleFilterX  - Constructor method
%   + predict        - Performs UKF prediction step
%   + update         - Performs UKF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1] M. S. Arulampalam, S. Maskell, N. Gordon and T. Clapp, "A tutorial on 
%     particle filters for online nonlinear/non-Gaussian Bayesian tracking,"
%     in IEEE Transactions on Signal Processing, vol. 50, no. 2, pp. 174-188, Feb 2002.
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties (Dependent)
        NumParticles
        MeasurementPrediction
        MeasurementLikelihoods
        MeasurementLikelihoodsPerParticle
    end
    
    properties
        StatePrior
        StatePrediction
        StatePosterior
        ControlInput
        ResamplingScheme = 'Systematic'
        ResamplingPolicy = {'TimeInterval',1}
        Resampler = SystematicResamplerX();
    end
    
    properties (Access=protected)
        MeasurementLikelihoods_ = [];
        MeasurementLikelihoodsPerParticle_ = [];
        MeasurementPrediction_ = [];
        NumParticles_ = 1000;
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            initialise_@FilterX(this,config);
            if (isfield(config,'StatePrior'))
                this.StatePrior = config.StatePrior;
                this.StatePosterior = this.StatePrior;
            end
            if (isfield(config,'ResamplingScheme'))
                this.ResamplingScheme = config.ResamplingScheme;
            end
            if (isfield(config,'ResamplingPolicy'))
                this.ResamplingPolicy = config.ResamplingPolicy;
            end
            if (isfield(config,'Resampler'))
                this.Resampler = config.Resampler;
            end
        end
        
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        % MeasurementList
        function measurementList = setMeasurementList(this, newMeasurementList)
            measurementList = newMeasurementList;
            this.MeasurementLikelihoodsPerParticle_ = [];
            this.MeasurementLikelihoods_ = [];
        end
        
        % StatePrior
        function StatePrior = setStatePrior(this,newStatePrior)
            if(isa(newStatePrior,'ParticleStateX'))
                this.NumParticles_ = newStatePrior.NumParticles;
            end
            StatePrior = newStatePrior;
        end
        
        % State Prediction
        function StatePrediction = setStatePrediction(this,newStatePrediction)
            if(isa(newStatePrediction,'ParticleStateX'))
                StatePrediction = newStatePrediction;
            else
                particles = newStatePrediction.Distribution.random(this.NumParticles_);
                weights = repmat(1/this.NumParticles_,1,this.NumParticles_);
                StatePrediction = ParticleStateX(particles,weights);
            end
            this.MeasurementPrediction_ = []; % Reset dependent measurement prediction
            this.MeasurementLikelihoods_ = [];
            this.MeasurementLikelihoodsPerParticle_ = [];
        end
        
        % Measurement Prediction
        function MeasurementPrediction = getMeasurementPrediction(this)
            if(isempty(this.MeasurementPrediction_))
                this.MeasurementPrediction_ =  ...
                    ParticleStateX(this.Model.Measurement.feval(this.StatePrediction.Particles,true),...
                                   this.StatePrediction.Weights);
            end
            MeasurementPrediction = this.MeasurementPrediction_;
        end
        function MeasurementPrediction = setMeasurementPrediction(this,newMeasurementPrediction)
            if(isa(newMeasurementPrediction,'ParticleStateX'))
                MeasurementPrediction = newMeasurementPrediction;
            else
                particles = newMeasurementPrediction.Distribution.random(this.NumParticles_);
                weights = repmat(1/this.NumParticles_,1,this.NumParticles_);
                MeasurementPrediction = ParticleStateX(particles,weights);
            end
            this.MeasurementLikelihood_ = [];
        end
        
        % MeasurementLikelihoods
        function MeasurementLikelihoods = getMeasurementLikelihoods(this)
            if(isempty(this.MeasurementLikelihoods_))
                this.MeasurementLikelihoods_ = mean(this.getMeasurementLikelihoodsPerParticle(),2)';
            end
            MeasurementLikelihoods = this.MeasurementLikelihoods_;
        end
        function MeasurementLikelihoods = setMeasurementLikelihoods(this, newMeasurementLikelihoods)
            MeasurementLikelihoods = newMeasurementLikelihoods;
        end
        
        % MeasurementLikelihoodsPerParticle
        function MeasurementLikelihoodsPerParticle = getMeasurementLikelihoodsPerParticle(this)
            if(isempty(this.MeasurementLikelihoodsPerParticle_))
                this.MeasurementLikelihoodsPerParticle_ = this.Model.Measurement.pdf(this.MeasurementList.Vectors,this.StatePrediction.Particles);
            end
            MeasurementLikelihoodsPerParticle = this.MeasurementLikelihoodsPerParticle_;
        end
        function MeasurementLikelihoodsPerParticle = setMeasurementLikelihoodsPerParticle(this, newMeasurementLikelihoodsPerParticle)
            MeasurementLikelihoodsPerParticle = newMeasurementLikelihoodsPerParticle;
        end
        
        % State Posterior
        function StatePosterior = setStatePosterior(this,newStatePosterior)
            if(isa(newStatePosterior,'ParticleStateX'))
                StatePosterior = newStatePosterior;
                this.NumParticles_ = StatePosterior.NumParticles;
            else
                particles = newStatePosterior.random(this.NumParticles_);
                weights = repmat(1/this.NumParticles_,1,this.NumParticles_);
                StatePosterior = ParticleStateX(particles,weights,newStatePosterior.Timestamp);
            end
        end
    end
    
    methods
        function this = ParticleFilterX(varargin)
        % ParticleFilterX Constructor method
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
        % * pf = ParticleFilterX(___,Name,Value) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update, smooth. 
                 
            % Call SuperClass method
            this@FilterX(varargin{:});
            
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
        % INITIALISE Initialise the Particle Filter with a certain set of
        % parameters.  
        %
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
                    
            initialise@FilterX(this);
            
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
        
        function varargout = predict(this, varargin)
        % predict Perform an SIR Particle Filter prediction step and return
        %   the generated state and measurement predictions.
        % 
        % Parameters
        % ----------
        % prior: StateX, optional
        %   The prior state estimate.
        % timestamp: datetime, optional
        %   A timestamp indicating the time at which prediction is
        %   performed.
        %
        % Returns
        % -------
        % ParticleStateX
        %   The generated state prediction
        % ParticleStateX, optional
        %   The generated measurement prediction
        %
        %  See also update, smooth.
            
            timestamp = [];
            timestamp_old = [];
            for i = 1:min([2,nargin-1])
                if isa(varargin{i},'StateX')
                    this.StatePosterior = varargin{i};
                    timestamp_old = this.StatePosterior.Timestamp;
                elseif isdatetime(varargin{i})
                    timestamp = varargin{i};
                end
            end
            
            if isempty(timestamp)
                dt = this.Model.Transition.TimestepDuration;
                timestamp = this.StatePosterior.Timestamp;
            else
                dt = timestamp - timestamp_old;
            end
                
            % Extract model parameters
            f  = @(x,wk) this.Model.Transition.feval(x, wk, dt); % Transition function
            wk = this.Model.Transition.random(this.StatePosterior.NumParticles, dt); % Process noise
        
            % Propagate particles through the dynamic model
            predParticles = ParticleFilterX_Predict(f,this.StatePosterior.Particles,wk);
            predWeights = this.StatePosterior.Weights;
            
            % Create a state prediction object
            this.StatePrediction = ParticleStateX(predParticles,...
                                                  predWeights,...
                                                  timestamp);
            % Return arguments
            varargout{1} = this.StatePrediction;
            if(nargout==2)
                % Optionally, compute and return measurement prediction
                varargout{2} = this.MeasurementPrediction;
            end
        end
        
        function posterior = update(this, varargin)
        % update Perform SIR Particle Filter update step
        % 
        % Parameters
        % ----------
        % prediction: StateX
        %   An object containing the prediction estimate
        % measurement: MeasurementX
        %   A single measurement/detection object.    
        %
        % Returns
        % -------
        % posterior: StateX
        %   The generated posterior state estimate.   
        %
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        % See also ParticleFilterX, predict, smooth.
            if nargin>1
                if isa(varargin{1},'MeasurementX')
                    this.MeasurementList = varargin{1};
                elseif isa(varargin{1}, 'StateX')
                    this.StatePrediction = varargin{1};
                    this.MeasurementList = varargin{2};
                end
            end
            if(this.MeasurementList.NumMeasurements)
                timestamp = this.MeasurementList.Timestamp;
            else
                timestamp = this.StatePrediction.Timestamp;
            end
            
            % Perform weight update
            weights = ...
                ParticleFilterX_UpdateWeights(@(y,x) this.Model.Measurement.pdf(y,x),...
                        this.MeasurementList.Vectors,this.StatePrediction.Particles,this.StatePrediction.Weights);
            
            % Resampling
            if(strcmp(this.ResamplingPolicy(1),"TimeInterval"))
                [particles,weights] = ...
                    this.Resampler.resample(this.StatePrediction.Particles,weights);
            end
            
            % Store and return posterior
            this.StatePosterior = ParticleStateX(particles,weights,timestamp);
            posterior = this.StatePosterior;
        end
        
        function posterior = updatePDA(this, assocWeights, measurementLikelihoodsPerParticle)
        % UPDATEPDA - Performs PF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % Usage
        % -----
        %  * updatePDA(assocWeights) Performs KF-PDA update step for multiple 
        %    measurements based on the provided (1-by-Nm+1) association weights 
        %    matrix assocWeights.
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        % See also KalmanFilterX, Predict, Iterate, Smooth, resample.
        
            if(nargin<3)
                measurementLikelihoodsPerParticle = this.MeasurementLikelihoodsPerParticle;  
            end
            if(this.MeasurementList.NumMeasurements)
                timestamp = this.MeasurementList.Timestamp;
            else
                timestamp = this.StatePrediction.Timestamp;
            end
            
            % Perform update
            weights = ...
                ParticleFilterX_UpdatePDA(@(y,x) this.Model.Measurement.feval(y,x),...
                                          this.MeasurementList.Vectors,this.StatePrediction.Particles,...
                                          this.StatePrediction.Weights, assocWeights, measurementLikelihoodsPerParticle);
            
            if(strcmp(this.ResamplingPolicy(1),"TimeInterval"))
                [particles,weights] = ...
                    this.Resampler.resample(this.StatePrediction.Particles,weights);
            end
            
            this.StatePosterior = ParticleStateX(particles,weights,timestamp);
            
            posterior = this.StatePosterior;
        end
        
        % ===============================>
        % ACCESS METHODS
        % ===============================>
        
        % StatePrior        
        function statePrior = get.StatePrior(this)
            statePrior = this.StatePrior;
        end
        function set.StatePrior(this, newStatePrior)
            this.StatePrior = setStatePrior(this, newStatePrior);
        end
        
        % StatePrediction
        function statePrediction = get.StatePrediction(this)
            statePrediction = this.StatePrediction;
        end
        function set.StatePrediction(this, newStatePrediction)
            this.StatePrediction = setStatePrediction(this, newStatePrediction);
        end
        
        % MeasurementPrediction
        function measurementPrediction = get.MeasurementPrediction(this)
            measurementPrediction = getMeasurementPrediction(this);
        end
        function set.MeasurementPrediction(this, newMeasurementPrediction)
            this.MeasurementPrediction_ = setMeasurementPrediction(this, newMeasurementPrediction);
        end
        
        % MeasurementLikelihood
        function MeasurementLikelihoods = get.MeasurementLikelihoods(this)
            MeasurementLikelihoods = getMeasurementLikelihoods(this);
        end
        function set.MeasurementLikelihoods(this, newMeasurementLikelihoods)
            this.MeasurementLikelihoods_ = setMeasurementLikelihoods(this, newMeasurementLikelihoods);
        end
        
        % MeasurementLikelihoodsPerParticle
        function MeasurementLikelihoodsPerParticle = get.MeasurementLikelihoodsPerParticle(this)
            MeasurementLikelihoodsPerParticle = getMeasurementLikelihoodsPerParticle(this);
        end
        function set.MeasurementLikelihoodsPerParticle(this, newMeasurementLikelihoodsPerParticle)
            this.MeasurementLikelihoodsPerParticle_ = setMeasurementLikelihoodsPerParticle(this, newMeasurementLikelihoodsPerParticle);
        end
        
        % StatePosterior
        function statePosterior = get.StatePosterior(this)
            statePosterior = this.StatePosterior;
        end
        function set.StatePosterior(this, newStatePosterior)
            this.StatePosterior = setStatePosterior(this, newStatePosterior);
        end
        
        % NumParticles
        function numParticles = get.NumParticles(this)
            numParticles = this.NumParticles_;
        end
        function set.NumParticles(this, numParticles)
            this.NumParticles_ = numParticles;
        end
    end
end