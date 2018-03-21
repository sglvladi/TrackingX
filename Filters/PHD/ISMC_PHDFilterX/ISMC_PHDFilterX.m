classdef ISMC_PHDFilterX < SMC_PHDFilterX
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

    properties
        Detections = [];
        NumBornParticles = 0;
        BornParticles = [];
        BornWeights = [];
        NumParticlesPerTarget = 1000;
        NumParticlesPerMeasurement = 50;
        ExpectedNumBornTargets = 0.2;
    end
    
    methods
        function this = ISMC_PHDFilterX(varargin)
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
            this@SMC_PHDFilterX(varargin{:});

            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'NumBornParticles'))
                        this.NumBornParticles  = config.NumBornParticles;
                    end
                    if (isfield(config,'NumParticlesPerTarget'))
                        this.NumParticlesPerTarget  = config.NumParticlesPerTarget;
                    end
                    if (isfield(config,'NumParticlesPerMeasurement'))
                        this.NumParticlesPerMeasurement  = config.NumParticlesPerMeasurement;
                    end
                    if (isfield(config,'ExpectedNumBornTargets'))
                        this.ExpectedNumBornTargets  = config.ExpectedNumBornTargets;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'NumBornParticles'))
                this.NumBornParticles  = config.NumBornParticles;
            end
            if (isfield(config,'NumParticlesPerTarget'))
                this.NumParticlesPerTarget  = config.NumParticlesPerTarget;
            end
            if (isfield(config,'NumParticlesPerMeasurement'))
                this.NumParticlesPerMeasurement  = config.NumParticlesPerMeasurement;
            end
            if (isfield(config,'ExpectedNumBornTargets'))
                this.ExpectedNumBornTargets  = config.ExpectedNumBornTargets;
            end
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
                        
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'NumBornParticles'))
                        this.NumBornParticles  = config.NumBornParticles;
                    end
                    if (isfield(config,'NumParticlesPerTarget'))
                        this.NumParticlesPerTarget  = config.NumParticlesPerTarget;
                    end
                    if (isfield(config,'NumParticlesPerMeasurement'))
                        this.NumParticlesPerMeasurement  = config.NumParticlesPerMeasurement;
                    end
                    if (isfield(config,'ExpectedNumBornTargets'))
                        this.ExpectedNumBornTargets  = config.ExpectedNumBornTargets;
                    end
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Results;
            if (isfield(config,'NumBornParticles'))
                this.NumBornParticles  = config.NumBornParticles;
            end
            if (isfield(config,'NumParticlesPerTarget'))
                this.NumParticlesPerTarget  = config.NumParticlesPerTarget;
            end
            if (isfield(config,'NumParticlesPerMeasurement'))
                this.NumParticlesPerMeasurement  = config.NumParticlesPerMeasurement;
            end
            if (isfield(config,'ExpectedNumBornTargets'))
                this.ExpectedNumBornTargets  = config.ExpectedNumBornTargets;
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
            
            this.NumMeasurements = size(this.MeasurementList,2);

            % Append born particles from previous iteration to
            % persistant particles (line 5 of Alg.1 in [2])
            this.Particles = [this.Particles, this.BornParticles];
            this.Weights = [this.Weights, this.BornWeights];
            this.NumParticles = this.NumParticles + this.NumBornParticles;

            % Predict persistant particles
            if(this.NumParticles>0)
                this.PredParticles = ...
                    this.Model.Dyn.feval(this.Particles, ...
                        this.Model.Dyn.random(this.NumParticles)); 
                this.PredWeights = ...
                    (1-this.ProbOfDeath)* this.Weights;
            end

            % Generate new born particles
            this.NumBornParticles = this.NumParticlesPerMeasurement*this.NumMeasurements;
            this.BornParticles = zeros(this.Model.Dyn.NumStateDims, this.NumBornParticles);
            this.BornWeights = zeros(1, this.NumBornParticles);
            for i = 1:size(this.MeasurementList,2)
                startInd = 1+this.NumParticlesPerMeasurement*(i-1);
                endInd = startInd + this.NumParticlesPerMeasurement-1;

                % Sample from the birth intensity
                this.BornParticles(:,startInd:endInd) = ...
                    this.BirthIntFcn(this.NumParticlesPerMeasurement,...
                        this.MeasurementList(:,i));
            end
            % Uniform weights across the new born particles
            this.BornWeights = ones(1,this.NumBornParticles)*...
                this.ExpectedNumBornTargets/this.NumBornParticles;
        end
        
        % Update function
        % ----------------
        % Performs the relevant SMC PHD update algo, based on the selected .type
        function update(this)
            
            % Initialise L_k(z) with 1st and 2nd term in Eq. (25) of [2]  
            Lk = zeros(1,this.NumMeasurements) + this.ClutterRate + sum(this.BornWeights);
            
            this.Detections = [];
            % Update persistant targets
            if(this.NumParticles>0)
                % Compute g(z|x) matrix as in [1] 
                g = this.Model.Obs.pdf(this.MeasurementList, this.PredParticles);
                
                % Compute 3rd term of L_k(z) in Eq. (25) of [2]
                for i = 1:this.NumMeasurements   % for all measurements
                    Lk(i) = Lk(i) + sum(this.ProbOfDetection*g(i,:).*this.PredWeights,2);
                end

                % Update weights Eq. (23) of [2]
                this.Weights = ((1-this.ProbOfDetection)+ ...
                    sum(this.ProbOfDetection*g./Lk(ones(1,this.NumParticles),:)',1))...
                    .*this.PredWeights;

                % Calculate w^{n,j} Eq. (21) of [2]
                w_nj = zeros(this.NumMeasurements, this.NumParticles);
                for j = 1:this.NumMeasurements
                    w_nj(j,:) = (this.ProbOfDetection*g(j,:)/Lk(j)).*this.PredWeights;
                end
                
                p_i = sum(w_nj,2)'
                for j = 1:this.NumMeasurements
                    ExistProb = sum(w_nj(j,:),2);
                    if(ExistProb>0.8)
                        warning('Detected target %d! %f',j,ExistProb);
                    end
                end

                % Resample
                weightSum = sum(this.Weights,2);
                this.NumTargets = round(weightSum);
                this.NumParticles = this.NumParticlesPerTarget*this.NumTargets;
                [this.Particles, this.Weights] = ...
                    this.Resampler.resample(this.PredParticles, ...
                        (this.Weights/weightSum),this.NumParticles); % Resample
                this.Weights = this.Weights*weightSum; % Rescale
            end
            
            % Update new born targets 
            % Update weights Eq. (24) of [2]
            this.BornWeights = sum(this.BornWeights(ones(1,this.NumMeasurements),:)...
                ./(Lk(ones(1,this.NumBornParticles),:)'),1);

            % Resample;
            weightSum = sum(this.BornWeights,2);
            [this.BornParticles, this.BornWeights] = ...
                this.Resampler.resample(this.BornParticles, ...
                    (this.BornWeights/weightSum),this.NumBornParticles); 
            %this.BornWeights = this.BornWeights*this.ExpectedNumBornTargets; 
        end
    end
end