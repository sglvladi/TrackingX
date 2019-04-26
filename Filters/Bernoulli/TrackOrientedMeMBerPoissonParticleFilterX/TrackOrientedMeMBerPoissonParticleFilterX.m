classdef TrackOrientedMeMBerPoissonParticleFilterX < FilterX
% TrackOrientedMeMBerPoissonParticleFilterX class
%
% Summary of TrackOrientedMeMBerPoissonGMFilterX:
% This is a class implementation of a Track-Oriented Marginal Hybrid 
% Multi-Bernoulli (MeMBer) Poisson Filter [1].
%
% TrackOrientedMeMBerPoissonGMFilterX Properties:
%   + MeasurementList - A (NumObsDims x NumObs) matrix used to store the received 
%                       measurements
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
%   + DetectionProbability - The probablity that a target will be detected in
%                       a given measurement scan.
%   + NumTargets - The estimated number of targets following an update step.
%   + Model - An object handle to StateSpaceModelX object
%       + Dyn - Object handle to DynamicModelX SubClass      
%       + Obs - Object handle to ObservationModelX SubClass 
%       + Ctr - Object handle to ControlModelX SubClass 
%
% TrackOrientedMeMBerPoissonGMFilterX Methods:
%   + TrackOrientedMeMBerPoissonGMFilterX - Constructor method
%   + predict - Performs Bernoulli prediction step
%   + update - Performs Bernoulli update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1] J. L. Williams, "Marginal multi-bernoulli filters: RFS derivation of MHT, JIPDA, and association-based member," in IEEE Transactions on Aerospace and Electronic Systems, vol. 51, no. 3, pp. 1664-1687, July 2015.
%
% See also ParticleFilterX, KalmanFilerX.

    properties (Dependent)
        NumComponents
        NumMeasurements
        NumTargets
        MeasurementLikelihoods
        MeasurementLikelihoodsPerComponent
    end
    
    properties
        NumBernoulliParticles = 5000
        NumPoissonParticles = 50000       
        
        StatePrior
        StatePrediction
        MeasurementPrediction
        StatePosterior
        SurvivalProbability  
        DetectionProbability 
        
        Bernoulli
        Poisson
        
        % Component management
        Filter
        Hypothesiser = LoopyBeliefPropagationX();
        
        ResamplingScheme = 'Systematic'
        ResamplingPolicy = {'TimeInterval',1}                    
        Resampler
        
        StoreTrajectories = true
        AssocLikelihoodMatrix
        MarginalAssociationProbabilities
    end
    
    properties (Access=protected)
        MeasurementLikelihoods_ = [];
        MeasurementLikelihoodsPerParticle_ = [];
        MeasurementPrediction_ = [];
        NumComponents_ = 1;
    end
    
    methods (Access=protected)
        function initialise_(this, config)
            initialise_@FilterX(this,config);
            
            this.Resampler = SystematicResamplerX();
            
            if (isfield(config,'StatePrior'))
                this.StatePrior = config.StatePrior;
                this.StatePosterior = this.StatePrior;
            end
            if (isfield(config,'SurvivalProbability'))
                this.SurvivalProbability = config.SurvivalProbability;
            end
            if (isfield(config,'DetectionProbability'))
                this.DetectionProbability = config.DetectionProbability;
            end
            if (isfield(config,'Thresholds'))
                thresholds = config.Thresholds;
                if (isfield(thresholds,'Deletion'))
                    if (isfield(thresholds.Deletion,'Bernoulli'))
                        this.Thresholds.Deletion.Bernoulli = thresholds.Deletion.Bernoulli;
                    end
                    if (isfield(thresholds.Deletion,'Poisson'))
                        this.Thresholds.Deletion.Poisson = thresholds.Deletion.Poisson;
                    end
                end
            end
            if (isfield(config,'Filter'))
                this.Filter = config.Filter;
            end
        end
        
        % StatePrior
        function StatePrior = setStatePrior(this,newStatePrior)
            StatePrior = newStatePrior;
        end
        
        % State Prediction
        function StatePrediction = setStatePrediction(this,newStatePrediction)
            if(isa(newStatePrediction,'GaussianMixtureStateX'))
                StatePrediction = newStatePrediction;
            else
                error('State type must be GaussianMixtureStateX');
            end
            this.MeasurementPrediction_ = []; % Reset dependent measurement prediction
            this.MeasurementLikelihoods_ = [];
            this.MeasurementLikelihoodsPerComponent_ = [];
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
            if(isa(newMeasurementPrediction,'GaussianMixtureStateX'))
                MeasurementPrediction = newMeasurementPrediction;
            else
                error('State type must be GaussianMixtureStateX');
            end
            this.MeasurementLikelihoods_ = [];
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
        
        % MeasurementLikelihoodsPerComponent
        function MeasurementLikelihoodsPerComponent = getMeasurementLikelihoodsPerComponent(this)
            if(isempty(this.MeasurementLikelihoodsPerComponent_))
                this.MeasurementLikelihoodsPerComponent_ = this.Model.Measurement.pdf(this.MeasurementList,this.StatePrediction.Means,this.StatePrediction.Covars);
            end
            MeasurementLikelihoodsPerComponent = this.MeasurementLikelihoodsPerComponent_;
        end
        function MeasurementLikelihoodsPerComponent = setMeasurementLikelihoodsPerComponent(this, newMeasurementLikelihoodsPerComponent)
            MeasurementLikelihoodsPerComponent = newMeasurementLikelihoodsPerComponent;
        end
        
        % State Posterior
        function StatePosterior = setStatePosterior(this,newStatePosterior)
            if(isa(newStatePosterior,'GaussianMixtureStateX'))
                StatePosterior = newStatePosterior;
            else
                error('State type must be GaussianMixtureStateX');
            end
        end
    end
    
    methods
        function this = TrackOrientedMeMBerPoissonParticleFilterX(varargin)
        % TrackOrientedMeMBerPoissonParticleFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
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
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % DetectionProbability: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        %
        % Usage
        % -----
        % * phd = TrackOrientedMeMBerPoissonGMFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * phd = TrackOrientedMeMBerPoissonGMFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * phd = TrackOrientedMeMBerPoissonGMFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * phd = TrackOrientedMeMBerPoissonGMFilterX(ssm,priorDistFcn) returns an object
        %   handle, preconfigured with the provided StateSpaceModel object 
        %   handle ssm and the prior information about the state, provided  
        %   in the form of the priorDistFcn function.
        % * phd = TrackOrientedMeMBerPoissonGMFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
            
           % Call SuperClass method
            this@FilterX(varargin{:});
            
            this.Filter = ParticleFilterX(varargin{:});
            
            
            this.Bernoulli.StatePosterior.Particles = zeros(this.Model.Transition.NumStateDims,this.NumBernoulliParticles,0);
            this.Bernoulli.StatePosterior.Weights = zeros(1,this.NumBernoulliParticles,0);
            this.Bernoulli.StatePosterior.ExistenceProbabilities = zeros(1,0);
            
            this.Poisson.StatePosterior.Particles = zeros(this.Model.Transition.NumStateDims,this.NumPoissonParticles);
            this.Poisson.StatePosterior.Weights = zeros(1,this.NumPoissonParticles);
            
            if(nargin==0)
                return;
            end
            
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
            config = parser.Results;
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
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % DetectionProbability: scalar
        %   The probablity that a target will be detected in a given measurement scan.
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
                    this.initialise_(config);
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
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
            
            % Interpret length of inputs
            bernoulli.StatePrediction = this.Bernoulli.StatePosterior;
            numBernoulli = numel(this.Bernoulli.StatePosterior.ExistenceProbabilities);
            
            % Predict existing Bernoulli tracks
            for i = 1:numBernoulli
                
                % Predict state particles
                bernoulli.StatePrediction.Particles(:,:,i) = ...
                  this.Model.Transition.feval(bernoulli.StatePrediction.Particles(:,:,i), true);
              
                % Predict state weights
                bernoulli.StatePrediction.Weights(:,:,i) = bernoulli.StatePrediction.Weights(:,:,i);
                
                % Predict existence probability
                bernoulli.StatePrediction.ExistenceProbabilities(i) = ...
                    this.SurvivalProbability*bernoulli.StatePrediction.ExistenceProbabilities(i);
                
                % Normalise state weights
                bernoulli.StatePrediction.Weights(:,:,i) = ...
                    bernoulli.StatePrediction.Weights(:,:,i)./sum(bernoulli.StatePrediction.Weights(:,:,i));
            end
            this.Bernoulli.StatePrediction = bernoulli.StatePrediction;

            % Predict existing PPP intensity
            poisson.StatePrediction = copy(this.Poisson.StatePosterior);
            
            % Number of birth particles
            Jk = 5000;

            % Generate Np normally predicted particles
            predParticles = this.Model.Transition.feval(poisson.StatePrediction.Particles, true); 
            predWeights = this.SurvivalProbability *poisson.StatePrediction.Weights;

            % Generate birth particles 
            [birthParticles, birthWeights] = this.Model.Birth.random(Jk);

            this.Poisson.StatePrediction = ParticleStateX([predParticles,birthParticles], [predWeights, birthWeights]);
            
        end
        
        function update(this)
        % UPDATE Perform Bernoulli Filter update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
        
            numMeasurements = size(this.MeasurementList,2);
            numBernoulli = numel(this.Bernoulli.StatePrediction.ExistenceProbabilities);
            
            % Compute Association Likelihood Matrix for known (Bernoulli) tracks
            this.AssocLikelihoodMatrix = zeros(numBernoulli+1, numMeasurements+1);
            
            % Compute measurement likelihoods
            g = zeros(numMeasurements,numBernoulli);
            for i = 1:numBernoulli
                
                g(:,i) = sum(this.Model.Measurement.pdf(this.MeasurementList,this.Bernoulli.StatePrediction.Particles(:,:,i)),2);
                
                % Missed detection hypothesis
                this.AssocLikelihoodMatrix(i+1,1) = ...
                    1 - this.Bernoulli.StatePrediction.ExistenceProbabilities(i) + this.Bernoulli.StatePrediction.ExistenceProbabilities(i)*(1-this.DetectionProbability);
                % True detection hypotheses
                for j = 1:numMeasurements
                    this.AssocLikelihoodMatrix(i+1,j+1) = this.Bernoulli.StatePrediction.ExistenceProbabilities(i)*this.DetectionProbability*g(j,i);
                end
            end
            
            % Create a new track for each measurement by updating PPP with measurement
            new_Bernoulli.Particles = zeros(this.Model.Transition.NumStateDims,this.NumBernoulliParticles,numMeasurements);
            new_Bernoulli.Weights = zeros(1,this.NumBernoulliParticles,numMeasurements);
            new_Bernoulli.ExistenceProbabilities = zeros(1,numMeasurements);
            
            % Compute measurement likelihoods
            g = this.Model.Measurement.pdf(this.MeasurementList,this.Poisson.StatePrediction.Particles);
                        
            % Compute C_k(z) Eq. (27) of [1]  
            Ck = this.DetectionProbability*g.*this.Poisson.StatePrediction.Weights;
            C = sum(Ck,2);
            C_plus = C + this.Model.Clutter.pdf(this.MeasurementList)';
            
            weightsPerHypothesis = [(1-this.DetectionProbability).*this.Poisson.StatePrediction.Weights;
                                          zeros(numMeasurements, this.Poisson.StatePrediction.NumParticles)];
            if(numMeasurements>0)
                % Compute C_k(z) Eq. (27) of [1]
                weightsPerHypothesis(2:end,:) = Ck./C_plus;
            end
            
            % For each measurement, create a new track as a mixture over
            % all PPP GM components
            for j = 1:numMeasurements
                
                % Association likelihood for new track
                this.AssocLikelihoodMatrix(1,j+1) = C_plus(j);
                
                % Weight for new track
                particles = this.Poisson.StatePrediction.Particles;
                weights = weightsPerHypothesis(j,:); %Ck(j)./C_plus(j);
                
                % Resample to correct number of particles
                % =======================================
                [new_Bernoulli.Particles(:,:,j),new_Bernoulli.Weights(1,:,j)] = this.Resampler.resample(particles,weights./sum(weights),this.NumBernoulliParticles);
                new_Bernoulli.ExistenceProbabilities(j) = C(j)/C_plus(j);
                
            end
            
            % Update (i.e., thin) intensity of unknown targets
            poisson.Posterior = copy(this.Poisson.StatePrediction);
            poisson.Posterior.Weights = (1-this.DetectionProbability)*poisson.Posterior.Weights; 
            this.Poisson.StatePosterior = poisson.Posterior;
            
            % Perform Data Association
            [a,b] = lbp(this.AssocLikelihoodMatrix(2:end,:),this.AssocLikelihoodMatrix(1,2:end));
            this.MarginalAssociationProbabilities = [[0,b'];a];
            
            bernoulli = this.formTracksTOMB(a,this.Bernoulli,b,new_Bernoulli);
            [this.Bernoulli,poisson] = this.recycle(bernoulli,this.Poisson,0.1);
            
            % Resample (equivalent to Step 3 of [1]
            numTargets = sum(poisson.StatePosterior.Weights,2); % Compute total mass
            [postParticles, postWeights, idx] = ...
                this.Resampler.resample(poisson.StatePosterior.Particles, ...
                    (poisson.StatePosterior.Weights/numTargets),this.NumPoissonParticles); % Resample
            postWeights = postWeights*numTargets; % Rescale
            this.Poisson.StatePosterior = ParticleStateX(postParticles,postWeights);
            
%             ss = bernoulli.StatePosterior.ExistenceProbabilities > 0.1;
%             bernoulli.StatePosterior.Particles = bernoulli.StatePosterior.Particles(:,:,ss);
%             bernoulli.StatePosterior.Weights = bernoulli.StatePosterior.Weights(:,:,ss);
%             bernoulli.StatePosterior.ExistenceProbabilities = bernoulli.StatePosterior.ExistenceProbabilities(ss);
%             this.Bernoulli = bernoulli;
%             this.Bernoulli.StatePosterior.Trajectories = bernoulli.StatePosterior.Trajectories(ss);
            
        end
        
        function  bernoulli = formTracksTOMB(this,pupd,bernoulli,pnew,bernoulli_new)
                        
            % Infer sizes
            nold = numel(bernoulli.StatePrediction.ExistenceProbabilities);
            nnew = numel(bernoulli_new.ExistenceProbabilities);
            
            bernoulli.StatePosterior.Particles = zeros(this.Model.Transition.NumStateDims,this.NumBernoulliParticles,nold+nnew);
            bernoulli.StatePosterior.Weights = zeros(1,this.NumBernoulliParticles,nold+nnew);
            bernoulli.StatePosterior.ExistenceProbabilities = zeros(1,nold+nnew);

            % Form continuing tracks
            for i = 1:nold
                
                % Existence
                pr = [bernoulli.StatePrediction.ExistenceProbabilities(i)*(1-this.DetectionProbability)/this.AssocLikelihoodMatrix(i+1,1)*pupd(i,1), pupd(i,2:end)];
                bernoulli.StatePosterior.ExistenceProbabilities(i) = sum(pr);
                assocWeigths = pr/bernoulli.StatePosterior.ExistenceProbabilities(i);
                
                % State Update
                this.Filter.StatePrediction = ParticleStateX(bernoulli.StatePrediction.Particles(:,:,i), ...
                                                             bernoulli.StatePrediction.Weights(:,:,i));
                this.Filter.MeasurementList = this.MeasurementList;
                this.Filter.updatePDA(assocWeigths);
                
                bernoulli.StatePosterior.Particles(:,:,i) = this.Filter.StatePosterior.Particles;
                bernoulli.StatePosterior.Weights(:,:,i) = this.Filter.StatePosterior.Weights;
                
                bernoulli.StatePosterior.Trajectories{i} = bernoulli.StatePrediction.Trajectories{i};
                bernoulli.StatePosterior.Trajectories{i}.StateMean(:,end+1) = bernoulli.StatePosterior.Particles(:,:,i)*bernoulli.StatePosterior.Weights(:,:,i)';
                bernoulli.StatePosterior.Trajectories{i}.StateCovar(:,:,end+1) = weightedcov(bernoulli.StatePosterior.Particles(:,:,i),bernoulli.StatePosterior.Weights(:,:,i));
            end

            % Form new tracks (already single hypothesis)
            for j = 1:nnew
                bernoulli.StatePosterior.ExistenceProbabilities(nold+j) = pnew(j) .* bernoulli_new.ExistenceProbabilities(j);
                bernoulli.StatePosterior.Particles(:,:,nold+j) = bernoulli_new.Particles(:,:,j);
                bernoulli.StatePosterior.Weights(:,:,nold+j) = bernoulli_new.Weights(:,:,j);
                bernoulli.StatePosterior.Trajectories{nold+j}.StateMean = bernoulli_new.Particles(:,:,j)*bernoulli_new.Weights(:,:,j)';
                bernoulli.StatePosterior.Trajectories{nold+j}.StateCovar = weightedcov(bernoulli_new.Particles(:,:,j),bernoulli_new.Weights(:,:,j));
            end
        end
        
        function [bernoulli, poisson] = recycle(this, bernoulli, poisson, threshold)
            
            recycleInds = find(bernoulli.StatePosterior.ExistenceProbabilities< threshold);
            numRecycle = numel(recycleInds);
            
            if(numRecycle>0)
                                
                % Append recycled bernoulli components
                for i=1:numRecycle
                    poisson.StatePosterior.Particles = [poisson.StatePosterior.Particles, bernoulli.StatePosterior.Particles(:,:,recycleInds(i))];
                    poisson.StatePosterior.Weights = ...
                        [poisson.StatePosterior.Weights, bernoulli.StatePosterior.Weights(:,:,recycleInds(i))*bernoulli.StatePosterior.ExistenceProbabilities(recycleInds(i))];
                end
                    
                % Remove old bernoulli components
                bernoulli.StatePosterior.Particles(:,:,recycleInds) = [];
                bernoulli.StatePosterior.Weights(:,:,recycleInds) = [];
                bernoulli.StatePosterior.ExistenceProbabilities(recycleInds) = [];
                bernoulli.StatePosterior.Trajectories(recycleInds) = [];
                
                
            end
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