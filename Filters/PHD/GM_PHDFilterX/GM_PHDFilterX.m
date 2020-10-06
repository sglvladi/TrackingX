classdef GM_PHDFilterX < FilterX
% GM_PHDFilterX class
%
% Summary of GM_PHDFilterX:
% This is a class implementation of a Gaussian Mixture (GM) Probabilistic
% Hypothesis Density (PHD) Filter.
%
% NOTE: Spawned targets are not currently implemented.
%
% GM_PHDFilterX Properties:
%   + NumComponents - The number of the last filtered GM components employed by the 
%                     PHD Filter
%   + StatePrior - A structure used to store the state prior
%   + StatePrediction - A structure used to store the state prediction
%   + MeasurementPrediction - A structure used to store the measurement prediction
%   + StatePosterior - A structure used to store posterior information  
%   + MeasurementList - A (NumMeasDims x NumMeasurements) matrix used to store
%                       the received measurements  
%   + Filter - A KalmanFilterX (sub)class instance, which is used to
%              predict and update all the GM components.
%   + Gater - A GaterX sub-class, used to perform gating between GM
%             components and measurements. (Optional, default = None)
%   + SurvivalProbability - The probability that a target may cease to exist
%                   between consecutive iterations of the filter.
%   + MaxNumComponents - The maximum allowable number of compunents following 
%                        an update() call 
%   + MergeThreshold - Minimum distance between the means of components. If
%                      the distance between the means of two (or more) components
%                      is less that MergeThreshold, the the Components will
%                      be merged during the next update() call.
%   + PruneThreshold - Minimum weight for components. Any components whose
%                      weight is below PruneThreshold, will be removed on
%                      the next update() call.
%   + Model - An object handle to StateSpaceModelX object
%       + Transition (*)  = Object handle to TransitionModelX SubClass      
%       + Measurement (*)  = Object handle to MeasurementModelX SubClass 
%       + Control (*)  = Object handle to ControlModelX SubClass
%       + Clutter
%       + Birth
%
% GM_PHDFilterX Methods:
%   + GM_PHDFilterX - Constructor method
%   + predict - Performs GM_PHD prediction step
%   + update - Performs GM_PHD update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
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
        MeasWeights = 1;
        
        % Component management
        Filter
        Gater
        MaxNumComponents = 100
        PruneThreshold = 1e-5
        MergeThreshold = 0.5
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
            if (isfield(config,'StatePrior'))
                this.StatePrior = config.StatePrior;
                this.StatePosterior = this.StatePrior;
            end
            if (isfield(config,'SurvivalProbability'))
                this.SurvivalProbability = config.SurvivalProbability;
            end
            if (isfield(config,'Filter'))
                this.Filter = config.Filter;
            end
            if (isfield(config,'Gater'))
                this.Gater = config.Gater;
            end
            if (isfield(config,'MaxNumComponents'))
                this.MaxNumComponents = config.MaxNumComponents;
            end
            if (isfield(config,'PruneThreshold'))
                this.PruneThreshold = config.PruneThreshold;
            end
            if (isfield(config,'MergeThreshold'))
                this.MergeThreshold = config.MergeThreshold;
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
                this.MeasurementLikelihoodsPerComponent_ = this.Model.Measurement.pdf(this.MeasurementList.Vectors,this.StatePrediction.Means,this.StatePrediction.Covars);
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
        function this = GM_PHDFilterX(varargin)
        % GM_PHDFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: struct, optional
        %   A GaussianMixtureStateX subclass object describing the state prior.
        % MaxNumComponents: scalar
        %   The maximum allowed number of mixture components to be employed 
        %   by the GM PHD Filter. 
        %   (default = 10000)
        % Gater: GaterX subclass instance, optional
        %   A pre-configured gater object, which should be an instance of 
        %   a subclass of the GaterX abstract class (e.g. EllipsoidalGaterX).
        %   If supplied, any values passed for ProbOfGating and/or
        %   GateLevel will be ignored.
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        %
        % Usage
        % -----
        % * phd = GM_PHDFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
             
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
        % INITIALISE Initialise the GM PHD Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % StatePrior: struct, optional
        %   A StateX subclass object describing the state prior. If StatePrior 
        %   is not a ParticleStateX instance, then it will be converted into
        %   one by sampling from the provided State's distribution.
        % MaxNumComponents: scalar
        %   The maximum allowed number of mixture components to be employed 
        %   by the GM PHD Filter. 
        %   (default = 10000)
        % Gater: GaterX subclass instance, optional
        %   A pre-configured gater object, which should be an instance of 
        %   a subclass of the GaterX abstract class (e.g. EllipsoidalGaterX).
        %   If supplied, any values passed for ProbOfGating and/or
        %   GateLevel will be ignored.
        % SurvivalProbability: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        %
        % Usage
        % -----
        % * initialise(phd,___,Name,Value,___) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value
        %   pair arguments.
        %
        %  See also predict, update.   
                         
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
         
        function predict(this)
        % PREDICT Perform GM PHD Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted PHD
        %
        % More details
        % ------------
        % * GM_PHDFilterX uses the Model class property, which should 
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
             
            % 1) Birth component extraction 
            % -------------------------- 
            
            prediction = copy(this.StatePosterior);
            numComponents = prediction.NumComponents;
            
            % Perform prediction
            for i = 1:numComponents
                 
                % Use underlying filter to predict components
                this.Filter.StatePosterior = GaussianStateX(prediction.Means(:,i),prediction.Covars(:,:,i));
                statePrediction = this.Filter.predictState();
                 
                % Assing predicted state to components
                prediction.Means(:,i) = statePrediction.Mean;
                prediction.Covars(:,:,i) = statePrediction.Covar;
            end
            prediction.Weights = this.SurvivalProbability.*prediction.Weights;
            
            % Add birth component to prediction
            prediction = prediction + this.Model.Birth.Distribution;
            
            this.StatePrediction = GaussianMixtureStateX(prediction);
                                             
            % 3) Measurement prediction
            % -------------------------
            numMeasDims = this.Model.Measurement.NumMeasDims;
            numStateDims = this.Model.Measurement.NumStateDims;
            numComponents = this.StatePrediction.NumComponents;
            
            % Preallocate memory
            Means = zeros(numMeasDims,numComponents);
            Covars = zeros(numMeasDims,numMeasDims,numComponents);
            Weights = zeros(1,numComponents);
            this.KalmanGains_=zeros(numStateDims,numMeasDims, numComponents);
            
            % Perform measurement prediction
            for i = 1:this.StatePrediction.NumComponents
                     
                % Use underlying filter to predict measurement
                this.Filter.StatePrediction = GaussianStateX(this.StatePrediction.Means(:,i),this.StatePrediction.Covars(:,:,i));
                measurementPrediction = this.Filter.predictMeasurement();
                 
                % Assigned measurement prediction to component
                Means(:,i) = measurementPrediction.Mean;
                Covars(:,:,i) = measurementPrediction.Covar;
                Weights(i) = this.StatePrediction.Weights(i);
                this.KalmanGains_(:,:,i) = this.Filter.KalmanGain;
            end
            
            this.MeasurementPrediction = GaussianMixtureStateX(Means,Covars,Weights);
        end
         
 
        function update(this)
        % UPDATE Perform GM PHD update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
                     
            % 4) Update
            % ---------
            numMeasurements = this.MeasurementList.NumMeasurements;
            numPredComponents = this.StatePrediction.NumComponents;
            numUpdateComponents = (1+numMeasurements)*numPredComponents;
            numStateDims = this.Model.Transition.NumStateDims;
            
            % Allocate memory
            Means = zeros(numStateDims,numUpdateComponents);
            Covars = zeros(numStateDims,numStateDims,numUpdateComponents);
            Weights = zeros(1, numUpdateComponents);
 
            % Create null measurement hypothesis components
            P_D = this.Model.Detection.pdf(this.StatePrediction.Means);
            Means(:,1:numPredComponents) = this.StatePrediction.Means;
            Covars(:,:,1:numPredComponents) = this.StatePrediction.Covars;
            Weights(1:numPredComponents) = (1-P_D).*this.StatePrediction.Weights;
             
            % Compute measurement likelihoods
            g = this.MeasurementLikelihoodsPerComponent;
             
            % Compute normalising constant
            Ck = P_D.*g.*this.StatePrediction.Weights;
            C = sum(Ck,2);
            Ck_plus = C + this.Model.Clutter.pdf(this.MeasurementList.Vectors)';
            
            % Compute the weights for all hypotheses (including null)
            weightsPerHypothesis_ = [Weights(1:numPredComponents);
                                     zeros(numMeasurements, numPredComponents)]; % Null
            weightsPerHypothesis_(2:end,:) = this.MeasWeights'.*Ck./Ck_plus; % True
%             a = sum(weightsPerHypothesis_,1);
            
            % Create true measurement hypothesis components
            for i = 1:numPredComponents
                 
                % Use underlying filter to update components
                this.Filter.StatePrediction = GaussianStateX(this.StatePrediction.Means(:,i),this.StatePrediction.Covars(:,:,i)) ;
                this.Filter.MeasurementPrediction = GaussianStateX(this.MeasurementPrediction.Means(:,i),this.MeasurementPrediction.Covars(:,:,i));
                this.Filter.KalmanGain = this.KalmanGains_(:,:,i);
                
                this.Filter.MeasurementList = this.MeasurementList;
                posterior = this.Filter.update();
                for j = 1:numMeasurements
                    
                    compIndex = numPredComponents + (i-1)*numMeasurements + j;
                    Means(:,compIndex) = posterior.Mean(:,j);
                    Covars(:,:,compIndex) = posterior.Covar;
                end
                
                % Compute corresponding weights
                Weights(numPredComponents + (i-1)*numMeasurements+1 : numPredComponents + i*numMeasurements) = ...
                     weightsPerHypothesis_(2:end,i);
            end  
             
            % 5) Prune
            % --------
            [Means,Covars,Weights] = ...
                gauss_prune(Means,Covars,...
                            Weights, this.PruneThreshold);
             
            % 6) Merge
            % --------
            [Means,Covars,Weights] = ...
               gauss_merge(Means,Covars,...
                           Weights, this.MergeThreshold);
                         
            % 7) Cap
            % ------
            [Means,Covars,Weights] = ...
                gauss_cap(Means,Covars,...
                          Weights, this.MaxNumComponents);
            
            this.StatePosterior = GaussianMixtureStateX(Means,Covars,Weights);
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
     
    methods (Access = protected)
         
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
         
    end
end