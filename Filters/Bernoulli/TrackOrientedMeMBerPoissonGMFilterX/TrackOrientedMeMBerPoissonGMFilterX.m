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

    properties (Access = private, Hidden)
        
    end
    properties
        Filter
        Prediction
        Posterior
        MeasurementList
        BirthModel
        ProbOfBirth
        ProbOfSurvive  
        ProbOfDetection 
        ClutterIntFcn
        AssocLikelihoodMatrix
        MarginalAssociationProbabilities
        Hypothesiser = LoopyBeliefPropagationX();
        Thresholds
        StoreTrajectories = true
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
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
        % ProbOfSurvive: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % ProbOfDetection: scalar
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
           
            this.Prediction.Bernoulli = {};
            this.Posterior.Bernoulli = {};
            this.Prediction.Poisson = {};
            this.Posterior.Poisson = {};
            
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
            
            % reset measLikelihood matrix
            this.pmeasLikelihood_ = [];
            
            % Interpret length of inputs
            numBernoulli = numel(this.Posterior.Bernoulli);
            numPoisson = numel(this.Posterior.Poisson);
            
            % Predict state pdf of existing tracks
            this.Prediction.Bernoulli = [];
            for i = 1:numBernoulli
                
                % Existence
                this.Prediction.Bernoulli{i}.ProbOfExistence = this.ProbOfSurvive*this.Posterior.Bernoulli{i}.ProbOfExistence;
                
                % PDF
                this.Filter.StateMean = this.Posterior.Bernoulli{i}.StateMean;
                this.Filter.StateCovar = this.Posterior.Bernoulli{i}.StateCovar;
                this.Filter.predict();
                this.Prediction.Bernoulli{i}.StateMean = this.Filter.PredStateMean;
                this.Prediction.Bernoulli{i}.StateCovar = this.Filter.PredStateCovar;
                this.Prediction.Bernoulli{i}.MeasMean = this.Filter.PredMeasMean;
                this.Prediction.Bernoulli{i}.InnovErrCovar = this.Filter.InnovErrCovar;
                this.Prediction.Bernoulli{i}.KalmanGain = this.Filter.KalmanGain;
                
                if(this.StoreTrajectories)
                    this.Prediction.Bernoulli{i}.Trajectory = this.Posterior.Bernoulli{i}.Trajectory;
                end
            end
            
            % Predict existing PPP intensity
            this.Prediction.Poisson = [];
            for k = 1:numPoisson
                
                % Intensity
                this.Prediction.Poisson{k}.Weight = this.ProbOfSurvive*this.Posterior.Poisson{k}.Weight;
                
                % PDF
                this.Filter.StateMean = this.Posterior.Poisson{k}.StateMean;
                this.Filter.StateCovar = this.Posterior.Poisson{k}.StateCovar;
                this.Filter.predictState();
                this.Prediction.Poisson{k}.StateMean = this.Filter.PredStateMean;
                this.Prediction.Poisson{k}.StateCovar = this.Filter.PredStateCovar;
            end
                        
            % Incorporate birth intensity into PPP
            [BirthMeans, BirthCovars, BirthWeights] = this.BirthModel.BirthIntFcn();
            numBirths = numel(BirthWeights);
            for k = 1:numBirths
                this.Prediction.Poisson{numPoisson+k}.Weight = BirthWeights(k);
                this.Prediction.Poisson{numPoisson+k}.StateMean = BirthMeans(:,k);
                this.Prediction.Poisson{numPoisson+k}.StateCovar = BirthCovars(:,:,k);
            end

            % Not shown in paper--truncate low weight components
            ss = cellfun(@(c) c.Weight, this.Prediction.Poisson)> this.Thresholds.Deletion.Poisson;
            this.Prediction.Poisson = this.Prediction.Poisson(ss);
            
            % Predict measurement for poisson (surviving and birth)
            for k = 1:numPoisson+numBirths
                this.Filter.PredStateMean = this.Prediction.Poisson{k}.StateMean;
                this.Filter.PredStateCovar = this.Prediction.Poisson{k}.StateCovar;
                this.Filter.predictObs();
                this.Prediction.Poisson{k}.MeasMean = this.Filter.PredMeasMean;
                this.Prediction.Poisson{k}.InnovErrCovar = this.Filter.InnovErrCovar;
                this.Prediction.Poisson{k}.KalmanGain = this.Filter.KalmanGain;
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
        
            [numMeasDims, numMeasurements] = size(this.MeasurementList);
           
            % Extract parameters from model
            H = this.Model.Obs.heval();
            Pd = this.ProbOfDetection;
            lambda_fa = this.ClutterIntFcn(1);
            
            % Interpret sizes from inputs
            numBernoulli = numel(this.Prediction.Bernoulli);
            numPoisson = numel(this.Prediction.Poisson);
            numStateDims = this.Model.Dyn.NumStateDims;
            
            this.AssocLikelihoodMatrix = zeros(numBernoulli+1, numMeasurements+1);
           
            % Prepare for data association
            for i = 1:numBernoulli
                
                % Missed detection hypothesis
                this.AssocLikelihoodMatrix(i+1,1) = 1 - this.Prediction.Bernoulli{i}.ProbOfExistence + this.Prediction.Bernoulli{i}.ProbOfExistence*(1-Pd);

                % True detection hypotheses
                try
                    g = mvnpdf(this.MeasurementList',this.Prediction.Bernoulli{i}.MeasMean',this.Prediction.Bernoulli{i}.InnovErrCovar)';
                catch
                    ads=1;
                end
                for j = 1:numMeasurements
                    this.AssocLikelihoodMatrix(i+1,j+1) = this.Prediction.Bernoulli{i}.ProbOfExistence*Pd*g(j);
                end
            end
            
            % Create a new track for each measurement by updating PPP with measurement
            new_Bernoulli = cell(numMeasurements,0);
            g = zeros(numPoisson,numMeasurements);
            for k = 1:numPoisson
                g(k,:) = mvnpdf(this.MeasurementList',this.Prediction.Poisson{k}.MeasMean',this.Prediction.Poisson{k}.InnovErrCovar)';
            end
            for j = 1:numMeasurements
                ck = zeros(1,numPoisson);
                yk = zeros(numStateDims,numPoisson);
                for k = 1:numPoisson
                    v = this.MeasurementList(:,j) - this.Prediction.Poisson{k}.MeasMean;
                    ck(k) = this.Prediction.Poisson{k}.Weight*Pd*g(k,j);
                    yk(:,k) = this.Prediction.Poisson{k}.StateMean + this.Prediction.Poisson{k}.KalmanGain*v;
                end
                C = sum(ck);
                this.AssocLikelihoodMatrix(1,j+1) = C + lambda_fa;
                new_Bernoulli{j}.ProbOfExistence = C/this.AssocLikelihoodMatrix(1,j+1);
                ck = ck/C;
                new_Bernoulli{j}.StateMean = yk*ck';
                new_Bernoulli{j}.StateCovar = zeros(numStateDims,numStateDims);
                for k = 1:numPoisson
                    v = new_Bernoulli{j}.StateMean - yk(:,k);
                    Pk = this.Prediction.Poisson{k}.StateCovar - this.Prediction.Poisson{k}.KalmanGain*H*this.Prediction.Poisson{k}.StateCovar;
                    new_Bernoulli{j}.StateCovar = new_Bernoulli{j}.StateCovar + ck(k)*(Pk + v*v');
                end
            end
            
            % Update (i.e., thin) intensity of unknown targets
            for k = 1:numPoisson
                this.Posterior.Poisson{k}.StateMean = this.Prediction.Poisson{k}.StateMean;
                this.Posterior.Poisson{k}.StateCovar = this.Prediction.Poisson{k}.StateCovar;
                this.Posterior.Poisson{k}.Weight = this.Prediction.Poisson{k}.Weight*(1-Pd);
            end

            % Not shown in paper--truncate low weight components
            ss = cellfun(@(c) c.Weight, this.Posterior.Poisson)> this.Thresholds.Deletion.Poisson;
            this.Posterior.Poisson = this.Posterior.Poisson(ss);
            
            % Perform Data Association
            [a,b] = lbp(this.AssocLikelihoodMatrix(2:end,:),this.AssocLikelihoodMatrix(1,2:end));
            this.MarginalAssociationProbabilities = [[0,b'];a];
            
            this.Posterior.Bernoulli = this.formTracksTOMB(a,this.Prediction.Bernoulli,b,new_Bernoulli);
            
            % Truncate tracks with low probability of existence (not shown in algorithm)
            ss = cellfun(@(c) c.ProbOfExistence,this.Posterior.Bernoulli)> this.Thresholds.Deletion.Bernoulli;
            this.Posterior.Bernoulli = this.Posterior.Bernoulli(ss);
        end
        
        function  bernoulli = formTracksTOMB(this,pupd,bernoulli,pnew,bernoulli_new)
            
            % Model parameters
            Pd = this.ProbOfDetection;
            
            % Infer sizes
            nold = numel(bernoulli);
            nnew = numel(bernoulli_new);

            % Form continuing tracks
            for i = 1:nold
                
                % Existence
                pr = [bernoulli{i}.ProbOfExistence*(1-Pd)/this.AssocLikelihoodMatrix(i+1,1)*pupd(i,1), pupd(i,2:end)];
                bernoulli{i}.ProbOfExistence = sum(pr);
                assocWeigths = pr/bernoulli{i}.ProbOfExistence;
                
                % State Update
                this.Filter.PredStateMean = bernoulli{i}.StateMean;
                this.Filter.PredStateCovar = bernoulli{i}.StateCovar;
                this.Filter.PredMeasMean = bernoulli{i}.MeasMean;
                this.Filter.InnovErrCovar = bernoulli{i}.InnovErrCovar; 
                this.Filter.KalmanGain = bernoulli{i}.KalmanGain;
                this.Filter.Measurement = this.MeasurementList;
                this.Filter.updatePDA(assocWeigths);
                
                bernoulli{i}.StateMean = this.Filter.StateMean;
                bernoulli{i}.StateCovar = this.Filter.StateCovar;
                
                bernoulli{i}.Trajectory.StateMean(:,end+1) = bernoulli{i}.StateMean;
            end

            % Form new tracks (already single hypothesis)
            for j = 1:nnew
                bernoulli{nold+j}.ProbOfExistence = pnew(j) .* bernoulli_new{j}.ProbOfExistence;
                bernoulli{nold+j}.StateMean = bernoulli_new{j}.StateMean;
                bernoulli{nold+j}.StateCovar = bernoulli_new{j}.StateCovar;
                bernoulli{nold+j}.Trajectory.StateMean = bernoulli{nold+j}.StateMean;
            end
        end
    end
    
    methods (Access = protected)
        
        function initialise_(this, config)
            
            % Apply defaults
            this.Thresholds.Deletion.Bernoulli = 1e-4;
            this.Thresholds.Deletion.Poisson = 1e-4;
            
            if (isfield(config,'Model'))
                this.Model = config.Model;
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
            if (isfield(config,'ClutterIntFcn'))
                this.ClutterIntFcn = config.ClutterIntFcn;
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
        end
        
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