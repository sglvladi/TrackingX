classdef TrackOrientedMeMBerPoissonGMFilterX < FilterX
% TrackOrientedMeMBerPoissonGMFilterX class
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
        Thresholds
        
        StoreTrajectories = true
        AssocLikelihoodMatrix
        MarginalAssociationProbabilities
    end
    
    properties (Access=protected)
        MeasurementLikelihoods_ = [];
        MeasurementLikelihoodsPerComponent_ = [];
        MeasurementPrediction_ = [];
        NumComponents_ = 1;
        KalmanGains_ = [];
    end
    
     methods (Access=protected)
        function initialise_(this, config)
            initialise_@FilterX(this,config);
            
            % Apply defaults
            this.Thresholds.Deletion.Bernoulli = 1e-4;
            this.Thresholds.Deletion.Poisson = 1e-4;
            
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
        function this = TrackOrientedMeMBerPoissonGMFilterX(varargin)
        % TrackOrientedMeMBerPoissonGMFilterX - Constructor method
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
           
            this.Bernoulli.StatePosterior = GaussianMixtureStateX();
            this.Poisson.StatePosterior = GaussianMixtureStateX();
            
            this.Filter = KalmanFilterX(varargin{:});
            
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
        % * initialise(phd,ssm) initialises the TrackOrientedMeMBerPoissonGMFilterX object 
        %   phd with the provided StateSpaceModelX object ssm.
        % * initialise(phd,priorParticles,priorWeights)initialises the 
        %   TrackOrientedMeMBerPoissonGMFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(phd,ssm,priorDistFcn) initialises the TrackOrientedMeMBerPoissonGMFilterX
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
        % * TrackOrientedMeMBerPoissonGMFilterX uses the Model class property, which should 
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
            bernoulli.StatePrediction = copy(this.Bernoulli.StatePosterior);
            poisson.StatePrediction = copy(this.Poisson.StatePosterior);
            
            numBernoulli = bernoulli.StatePrediction.NumComponents;
            numPoisson = poisson.StatePrediction.NumComponents;
            numStateDims = this.Model.Transition.NumStateDims;
            numMeasDims = this.Model.Measurement.NumMeasDims;
            
            bernoulli.MeasurementPrediction = GaussianMixtureStateX('empty',numMeasDims,numBernoulli);
            bernoulli.MeasurementPrediction.addprop('KalmanGains');
            bernoulli.MeasurementPrediction.KalmanGains = zeros(numStateDims,numMeasDims,numBernoulli);
            
            % Predict state pdf of existing tracks
            for i = 1:numBernoulli
                
                % PDF
                this.Filter.StatePosterior = GaussianStateX(bernoulli.StatePrediction.Means(:,i),...
                                                            bernoulli.StatePrediction.Covars(:,:,i));
                [statePrediction_i, measurementPrediction_i] = this.Filter.predict();
                bernoulli.StatePrediction.Means(:,i) = statePrediction_i.Mean;
                bernoulli.StatePrediction.Covars(:,:,i) = statePrediction_i.Covar;
                bernoulli.MeasurementPrediction.Means(:,i) = measurementPrediction_i.Mean;
                bernoulli.MeasurementPrediction.Covars(:,:,i) = measurementPrediction_i.Covar;
                bernoulli.MeasurementPrediction.KalmanGains(:,:,i) = this.Filter.KalmanGain;
                
            end
            % Existence
            bernoulli.StatePrediction.Weights = this.SurvivalProbability.*bernoulli.StatePrediction.Weights;
            this.Bernoulli.StatePrediction = bernoulli.StatePrediction;
            
            bernoulli.MeasurementPrediction.Weights = bernoulli.StatePrediction.Weights;
            this.Bernoulli.MeasurementPrediction = bernoulli.MeasurementPrediction;
            
            
            % Predict existing PPP intensity
            for k = 1:numPoisson
                % Use underlying filter to predict components
                this.Filter.StatePosterior = GaussianStateX(poisson.StatePrediction.Means(:,k),...
                                                            poisson.StatePrediction.Covars(:,:,k));
                statePrediction_k = this.Filter.predictState();
                
                % Assing predicted state to components
                poisson.StatePrediction.Means(:,k) = statePrediction_k.Mean;
                poisson.StatePrediction.Covars(:,:,k) = statePrediction_k.Covar;
            end
            % Intensity
            poisson.StatePrediction.Weights = this.SurvivalProbability.*poisson.StatePrediction.Weights;
                        
            % Add birth components to PPP
            poisson.StatePrediction = poisson.StatePrediction + this.Model.Birth.Distribution;

            % Not shown in paper--truncate low weight components
            ss = poisson.StatePrediction.Weights > this.Thresholds.Deletion.Poisson;
            this.Poisson.StatePrediction =  GaussianMixtureStateX(poisson.StatePrediction.Means(:,ss),...
                                                                  poisson.StatePrediction.Covars(:,:,ss), ...
                                                                  poisson.StatePrediction.Weights(ss));
            
            % Predict measurement for poisson (surviving and birth)
            numPoisson = this.Poisson.StatePrediction.NumComponents;
            poisson.MeasurementPrediction = GaussianMixtureStateX('empty',numMeasDims,numPoisson);
            poisson.MeasurementPrediction.addprop('KalmanGains');
            poisson.MeasurementPrediction.KalmanGains = zeros(numStateDims,numMeasDims,numPoisson);
            
            % Perform measurement prediction
            for k = 1:poisson.StatePrediction.NumComponents
                % Use underlying filter to predict measurement
                this.Filter.StatePrediction = GaussianStateX(poisson.StatePrediction.Means(:,k),poisson.StatePrediction.Covars(:,:,k));
                measurementPrediction_k = this.Filter.predictMeasurement();
                
                poisson.MeasurementPrediction.Means(:,k) = measurementPrediction_k.Mean;
                poisson.MeasurementPrediction.Covars(:,:,k) = measurementPrediction_k.Covar;
                poisson.MeasurementPrediction.KalmanGains(:,:,k) = this.Filter.KalmanGain;
            end
            poisson.MeasurementPrediction.Weights = poisson.StatePrediction.Weights;
            this.Poisson.MeasurementPrediction = copy(poisson.MeasurementPrediction);
        end
        
        function update(this)
        % UPDATE Perform Bernoulli Filter update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
                   
            % Interpret sizes from inputs
            numMeasurements = size(this.MeasurementList,2);
            numBernoulli = this.Bernoulli.StatePrediction.NumComponents;
            numPoisson = this.Poisson.StatePrediction.NumComponents;
            numStateDims = this.Model.Transition.NumStateDims;
            
            % Compute Association Likelihood Matrix for known (Bernoulli) tracks
            this.AssocLikelihoodMatrix = zeros(numBernoulli+1, numMeasurements+1);
            
            % Compute measurement likelihoods
            g = this.Model.Measurement.pdf(this.MeasurementList,this.Bernoulli.StatePrediction.Means,this.Bernoulli.StatePrediction.Covars);
            for i = 1:numBernoulli
                % Missed detection hypothesis
                this.AssocLikelihoodMatrix(i+1,1) = ...
                    1 - this.Bernoulli.StatePrediction.Weights(i) + this.Bernoulli.StatePrediction.Weights(i)*(1-this.DetectionProbability);
                % True detection hypotheses
                for j = 1:numMeasurements
                    this.AssocLikelihoodMatrix(i+1,j+1) = this.Bernoulli.StatePrediction.Weights(i)*this.DetectionProbability*g(j,i);
                end
            end
            
            % Create a new track for each measurement by updating PPP with measurement
            new_Bernoulli = GaussianMixtureStateX('empty',numStateDims,numMeasurements);

            % Compute measurement likelihoods
            g = this.Model.Measurement.pdf(this.MeasurementList,this.Poisson.StatePrediction.Means,this.Poisson.StatePrediction.Covars);

            % Compute normalising constant(s)
            Ck = this.DetectionProbability*g.*this.Poisson.StatePrediction.Weights;
            C = sum(Ck,2);
            C_plus = C + this.Model.Clutter.pdf(this.MeasurementList)';
            
            % Pre-allocate memory
            Means = zeros(numStateDims,numPoisson,numMeasurements);
            Covars = zeros(numStateDims,numStateDims,numPoisson);
            Weights = Ck./C;
            
            % Standard Kalman Update for all PPP components and mesurements
            this.Filter.MeasurementList = this.MeasurementList;
            for k = 1:numPoisson
                
                % Use underlying filter to update components
                this.Filter.StatePrediction = GaussianStateX(this.Poisson.StatePrediction.Means(:,k),this.Poisson.StatePrediction.Covars(:,:,k)) ;
                this.Filter.MeasurementPrediction = GaussianStateX(this.Poisson.MeasurementPrediction.Means(:,k),this.Poisson.MeasurementPrediction.Covars(:,:,k));
                this.Filter.KalmanGain = this.Poisson.MeasurementPrediction.KalmanGains(:,:,k);
                
                posterior = this.Filter.update();
                Covars(:,:,k) = posterior.Covar;
                for j = 1:numMeasurements
                    Means(:,k,j) = posterior.Mean(:,j);
                end
            end
            
            % For each measurement, create a new track as a mixture over
            % all PPP GM components
            for j = 1:numMeasurements
                
                % Association likelihood for new track
                this.AssocLikelihoodMatrix(1,j+1) = C_plus(j);
                
                % Weight for new track
                new_Bernoulli.Weights(j) = C(j)/C_plus(j);
                
                % Reduce mixture to single component
                % ==================================
                % Normalise ck to get mixture weights
                gm = GaussianMixtureX(Means(:,:,j), Covars, Weights(j,:));
                
                % Weighted mean over all components for that measurement
                new_Bernoulli.Means(:,j) = gm.Mean;
                new_Bernoulli.Covars(:,:,j) = gm.Covar;
                new_Bernoulli.Weights(j) = C(j)/C_plus(j);
            end
            
            % Update (i.e., thin) intensity of unknown targets
            poisson.Posterior = copy(this.Poisson.StatePrediction);
            poisson.Posterior.Weights = (1-this.DetectionProbability)*poisson.Posterior.Weights;

            % Not shown in paper--truncate low weight components
            ss = poisson.Posterior.Weights > this.Thresholds.Deletion.Poisson;
            this.Poisson.StatePosterior =  GaussianMixtureStateX(poisson.Posterior.Means(:,ss),...
                                                                 poisson.Posterior.Covars(:,:,ss), ...
                                                                 poisson.Posterior.Weights(ss));
            
            % Perform Data Association
            [a,b] = lbp(this.AssocLikelihoodMatrix(2:end,:),this.AssocLikelihoodMatrix(1,2:end));
            this.MarginalAssociationProbabilities = [[0,b'];a];
            
            bernoulli = this.formTracksTOMB(a,this.Bernoulli,b,new_Bernoulli);
            
            % Truncate tracks with low probability of existence (not shown in algorithm)
            ss = bernoulli.StatePosterior.Weights > this.Thresholds.Deletion.Bernoulli;
            this.Bernoulli.StatePosterior = GaussianMixtureStateX(bernoulli.StatePosterior.Means(:,ss),...
                                                                  bernoulli.StatePosterior.Covars(:,:,ss), ...
                                                                  bernoulli.StatePosterior.Weights(ss));
            this.Bernoulli.StatePosterior.addprop('Trajectories');
            this.Bernoulli.StatePosterior.Trajectories = bernoulli.StatePosterior.Trajectories(ss);
%                                                               
%             this.Bernoulli.StatePosterior = GaussianMixtureStateX(bernoulli.StatePosterior.Means,...
%                                                                   bernoulli.StatePosterior.Covars, ...
%                                                                   bernoulli.StatePosterior.Weights);
%             this.Bernoulli.StatePosterior.addprop('Trajectories');
%             this.Bernoulli.StatePosterior.Trajectories = bernoulli.StatePosterior.Trajectories;
            %[this.Bernoulli, this.Poisson] = this.recycle(this.Bernoulli, this.Poisson, 0.001);
            
        end
        
        function  bernoulli = formTracksTOMB(this,pupd,bernoulli,pnew,bernoulli_new)
                        
            % Infer sizes
            nold = bernoulli.StatePrediction.NumComponents;
            nnew = bernoulli_new.NumComponents;
            
            bernoulli.StatePosterior = GaussianMixtureStateX('empty', this.Model.Transition.NumStateDims, nold+nnew);
            bernoulli.StatePosterior.addprop('Trajectories');
            bernoulli.StatePosterior.Trajectories = cell(1,nold+nnew);

            % Form continuing tracks
            for i = 1:nold
                
                % Existence
                pr = [bernoulli.StatePrediction.Weights(i)*(1-this.DetectionProbability)/this.AssocLikelihoodMatrix(i+1,1)*pupd(i,1), pupd(i,2:end)];
                bernoulli.StatePosterior.Weights(i) = sum(pr);
                assocWeigths = pr/bernoulli.StatePosterior.Weights(i);
                
                % State Update
                this.Filter.StatePrediction = GaussianStateX(bernoulli.StatePrediction.Means(:,i), ...
                                                             bernoulli.StatePrediction.Covars(:,:,i));
                this.Filter.MeasurementPrediction = GaussianStateX(bernoulli.MeasurementPrediction.Means(:,i), ...
                                                                   bernoulli.MeasurementPrediction.Covars(:,:,i));
                this.Filter.KalmanGain = bernoulli.MeasurementPrediction.KalmanGains(:,:,i);
                this.Filter.MeasurementList = this.MeasurementList;
                this.Filter.updatePDA(assocWeigths);
                
                bernoulli.StatePosterior.Means(:,i) = this.Filter.StatePosterior.Mean;
                bernoulli.StatePosterior.Covars(:,:,i) = this.Filter.StatePosterior.Covar;
                
                bernoulli.StatePosterior.Trajectories{i} = bernoulli.StatePrediction.Trajectories{i};
                bernoulli.StatePosterior.Trajectories{i}.StateMean(:,end+1) = bernoulli.StatePosterior.Means(:,i);
                bernoulli.StatePosterior.Trajectories{i}.StateCovar(:,:,end+1) = bernoulli.StatePosterior.Covars(:,:,i);
            end

            % Form new tracks (already single hypothesis)
            for j = 1:nnew
                bernoulli.StatePosterior.Weights(nold+j) = pnew(j) .* bernoulli_new.Weights(j);
                bernoulli.StatePosterior.Means(:,nold+j) = bernoulli_new.Means(:,j);
                bernoulli.StatePosterior.Covars(:,:,nold+j) = bernoulli_new.Covars(:,:,j);
                bernoulli.StatePosterior.Trajectories{nold+j}.StateMean = bernoulli_new.Means(:,j);
                bernoulli.StatePosterior.Trajectories{nold+j}.StateCovar = bernoulli_new.Covars(:,:,j);
            end
        end
        
        function [bernoulli, poisson] = recycle(this, bernoulli, poisson, threshold)
            numBernoulli = bernoulli.StatePosterior.NumComponents;
            numPoisson = poisson.StatePosterior.NumComponents;
            numStateDims = this.Model.Transition.NumStateDims;
            
            recycleInds = find(bernoulli.StatePosterior.Weights< threshold);
            numRecycle = numel(recycleInds);
            
            if(numRecycle>0)
                
                % Copy over old poisson components
                poisson_posterior = GaussianMixtureStateX('empty', numStateDims, numPoisson+numRecycle);
                poisson_posterior.Means(:,1:numPoisson) = poisson.StatePosterior.Means;
                poisson_posterior.Covars(:,:,1:numPoisson) = poisson.StatePosterior.Covars;
                poisson_posterior.Weights(1:numPoisson) = poisson.StatePosterior.Weights;
                
                % Append recycled bernoulli components
                poisson_posterior.Means(:,numPoisson+1:end) = bernoulli.StatePosterior.Means(:,recycleInds);
                poisson_posterior.Covars(:,:,numPoisson+1:end) = bernoulli.StatePosterior.Covars(:,:,recycleInds);
                poisson_posterior.Weights(numPoisson+1:end) = bernoulli.StatePosterior.Weights(recycleInds);
                
                for i=1:numPoisson+numRecycle
                    [S,a]=cholcov(poisson_posterior.Covars(:,:,i));
                    if(a~=0)
                        disp(i);
                    end
                end
                
                % Remove old bernoulli components
                bernoulli.StatePosterior.Means(:,recycleInds) = [];
                bernoulli.StatePosterior.Covars(:,:,recycleInds) = [];
                bernoulli.StatePosterior.Weights(recycleInds) = [];
                bernoulli.StatePosterior.Trajectories(recycleInds) = [];
                
                poisson.StatePosterior = poisson_posterior;

            end
        end
        
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
            measurementPrediction = this.MeasurementPrediction;
        end
        function set.MeasurementPrediction(this, newMeasurementPrediction)
            this.MeasurementPrediction = setMeasurementPrediction(this, newMeasurementPrediction);
        end
        
        % MeasurementLikelihood
        function MeasurementLikelihoods = get.MeasurementLikelihoods(this)
            MeasurementLikelihoods = getMeasurementLikelihoods(this);
        end
        function set.MeasurementLikelihoods(this, newMeasurementLikelihoods)
            this.MeasurementLikelihoods_ = setMeasurementLikelihoods(this, newMeasurementLikelihoods);
        end
        
        % MeasurementLikelihoodsPerComponent
        function MeasurementLikelihoodsPerComponent = get.MeasurementLikelihoodsPerComponent(this)
            MeasurementLikelihoodsPerComponent = getMeasurementLikelihoodsPerComponent(this);
        end
        function set.MeasurementLikelihoodsPerComponent(this, newMeasurementLikelihoodsPerComponent)
            this.MeasurementLikelihoodsPerComponent_ = setMeasurementLikelihoodsPerComponent(this, newMeasurementLikelihoodsPerComponent);
        end
        
        % StatePosterior
        function statePosterior = get.StatePosterior(this)
            statePosterior = this.StatePosterior;
        end
        function set.StatePosterior(this, newStatePosterior)
            this.StatePosterior = setStatePosterior(this, newStatePosterior);
        end
    end
end