classdef GM_PHDFilterX < FilterX
% SMC_PHDFilterX class
%
% Summary of SMC_PHDFilterX:
% This is a class implementation of a Sequential Monte Carlo (SMC) Probabilistic
% Hypothesis Density (PHD) Filter.
%
% SMC_PHDFilterX Properties:
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
%   + ProbOfDetection - The probablity that a target will be detected in
%                       a given measurement scan.
%   + NumTargets - The estimated number of targets following an update step.
%   + Model - An object handle to StateSpaceModelX object
%       + Dyn - Object handle to DynamicModelX SubClass      
%       + Obs - Object handle to ObservationModelX SubClass 
%       + Ctr - Object handle to ControlModelX SubClass 
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

    properties (Access = private, Hidden)
        
    end
    properties
        Components
        NumComponents
        PredComponents
        NumPredComponents
        MaxNumComponents = 100
        PruneThreshold = 1e-5
        MergeThreshold = 4
        Gater
        BirthIntFcn 
        ProbOfDeath  
        ProbOfDetection 
        ClutterRate 
        MeasLikelihood
        MeasWeights = 1;
        WeightsPerMeasurement
        NumTargets = 0       
        MeasurementList
        NumMeasurements
        Filter
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
    end
    
    methods
        function this = GM_PHDFilterX(varargin)
        % GM_PHDFilterX - Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % MaxNumComponents: scalar
        %   The maximum allowed number of mixture components to be employed 
        %   by the GM PHD Filter. 
        %   (default = 10000)
        % PriorComponents: (1 x NumPriorComponents) cell array
        %   A cell array of Component strcuctures, where each Component
        %   structure has the following fields:
        %       - Mean: (NumStateDims x 1) column vector
        %           The mean of the Mixture component
        %       - Covar: (NumStateDims x NumStateDims) matrix
        %           The covariance of the Mixture component
        %       - Weight: scalar
        %           The weight of the Mixture component
        % Gater: GaterX subclass instance, optional
        %   A pre-configured gater object, which should be an instance of 
        %   a subclass of the GaterX abstract class (e.g. EllipsoidalGaterX).
        %   If supplied, any values passed for ProbOfGating and/or
        %   GateLevel will be ignored.
        % GateLevel: scalar, optional
        %   The desired number of standard deviations (away from the mean),
        %   which should be used as a gating threshold. 
        %   Only valid if a Gater is not supplied as a parameter.
        %   If supplied, any value supplied for ProbOfGating will be overriden 
        %   and recomputed based on this value.
        % ProbOfGating: scalar, optional
        %   The normalised ([0,1]) percentage of the predicted measurement 
        %   pdf, that should be incorporated in the validation region 
        %   (i.e. the probability of gating).
        %   Only valid if neither a Gater nor a GateLevel are not supplied 
        %   as parameters.
        % BirthIntFcn: function handle
        %   A function handle, BirthComponents = BirthIntFcn(), 
        %   which when called generates a set of birth components.
        % ProbOfDeath: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % ProbOfDetection: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        %
        % Usage
        % -----
        % * phd = GM_PHDFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * phd = GM_PHDFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * phd = GM_PHDFilterX(ssm,PriorComponents) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the PriorComponents.
        % * phd = GM_PHDFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also predict, update.   
            
            % Call SuperClass method
            this@FilterX(varargin{:});


            if(nargin==0)
                this.Gater = EllipsoidalGaterX('GateLevel',10);
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'PriorComponents'))
                        this.Components  = config.PriorComponents;
                    end
                    if (isfield(config,'Gater'))
                        this.Gater = config.Gater;
                    elseif (isfield(config,'GateLevel'))
                        this.Gater = EllipsoidalGaterX('GateLevel',config.GateLevel);
                    else
                        this.Gater = EllipsoidalGaterX('ProbOfGating',config.ProbOfGating);
                    end
                    if (isfield(config,'BirthIntFcn'))
                        this.BirthIntFcn = config.BirthIntFcn;
                    end
                    if (isfield(config,'ProbOfDeath'))
                        this.ProbOfDeath = config.ProbOfDeath;
                    end
                    if (isfield(config,'ProbOfDetection'))
                        this.ProbOfDetection = config.ProbOfDetection;
                    end
                    if (isfield(config,'ClutterRate'))
                        this.ClutterRate = config.ClutterRate;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            if (isfield(config,'PriorComponents'))
                this.Components  = config.PriorComponents;
            end
            if (isfield(config,'Gater'))
                this.Gater = config.Gater;
            elseif (isfield(config,'GateLevel'))
                this.Gater = EllipsoidalGaterX('GateLevel',config.GateLevel);
            else
                this.Gater = EllipsoidalGaterX('ProbOfGating',config.ProbOfGating);
            end
            if (isfield(config,'BirthIntFcn'))
                this.BirthIntFcn = config.BirthIntFcn;
            end
            if (isfield(config,'ProbOfDeath'))
                this.ProbOfDeath = config.ProbOfDeath;
            end
            if (isfield(config,'ProbOfDetection'))
                this.ProbOfDetection = config.ProbOfDetection;
            end
            if (isfield(config,'ClutterRate'))
                this.ClutterRate = config.ClutterRate;
            end
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the SMC PHD Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % MaxNumComponents: scalar
        %   The maximum allowed number of mixture components to be employed 
        %   by the GM PHD Filter. 
        %   (default = 10000)
        % PriorComponents: (1 x NumPriorComponents) cell array
        %   A cell array of Component strcuctures, where each Component
        %   structure has the following fields:
        %       - Mean: (NumStateDims x 1) column vector
        %           The mean of the Mixture component
        %       - Covar: (NumStateDims x NumStateDims) matrix
        %           The covariance of the Mixture component
        %       - Weight: scalar
        %           The weight of the Mixture component
        % Gater: GaterX subclass instance, optional
        %   A pre-configured gater object, which should be an instance of 
        %   a subclass of the GaterX abstract class (e.g. EllipsoidalGaterX).
        %   If supplied, any values passed for ProbOfGating and/or
        %   GateLevel will be ignored.
        % GateLevel: scalar, optional
        %   The desired number of standard deviations (away from the mean),
        %   which should be used as a gating threshold. 
        %   Only valid if a Gater is not supplied as a parameter.
        %   If supplied, any value supplied for ProbOfGating will be overriden 
        %   and recomputed based on this value.
        % ProbOfGating: scalar, optional
        %   The normalised ([0,1]) percentage of the predicted measurement 
        %   pdf, that should be incorporated in the validation region 
        %   (i.e. the probability of gating).
        %   Only valid if neither a Gater nor a GateLevel are not supplied 
        %   as parameters.
        % BirthIntFcn: function handle
        %   A function handle, BirthComponents = BirthIntFcn(), 
        %   which when called generates a set of birth components.
        % ProbOfDeath: scalar
        %   The probability that a target may cease to exist between consecutive 
        %   iterations of the filter.
        % ProbOfDetection: scalar
        %   The probablity that a target will be detected in a given measurement scan.
        %
        % Usage
        % -----
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
        %  See also predict, update.   
                        
             % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'PriorComponents'))
                        this.Components  = config.PriorComponents;
                    end
                    if (isfield(config,'Gater'))
                        this.Gater = config.Gater;
                    elseif (isfield(config,'GateLevel'))
                        this.Gater = EllipsoidalGaterX('GateLevel',config.GateLevel);
                    else
                        this.Gater = EllipsoidalGaterX('ProbOfGating',config.ProbOfGating);
                    end
                    if (isfield(config,'BirthIntFcn'))
                        this.BirthIntFcn = config.BirthIntFcn;
                    end
                    if (isfield(config,'ProbOfDeath'))
                        this.ProbOfDeath = config.ProbOfDeath;
                    end
                    if (isfield(config,'ProbOfDetection'))
                        this.ProbOfDetection = config.ProbOfDetection;
                    end
                    if (isfield(config,'ClutterRate'))
                        this.ClutterRate = config.ClutterRate;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            if (isfield(config,'PriorComponents'))
                this.Components  = config.PriorComponents;
            end
            if (isfield(config,'Gater'))
                this.Gater = config.Gater;
            elseif (isfield(config,'GateLevel'))
                this.Gater = EllipsoidalGaterX('GateLevel',config.GateLevel);
            else
                this.Gater = EllipsoidalGaterX('ProbOfGating',config.ProbOfGating);
            end
            if (isfield(config,'BirthIntFcn'))
                this.BirthIntFcn = config.BirthIntFcn;
            end
            if (isfield(config,'ProbOfDeath'))
                this.ProbOfDeath = config.ProbOfDeath;
            end
            if (isfield(config,'ProbOfDetection'))
                this.ProbOfDetection = config.ProbOfDetection;
            end
            if (isfield(config,'ClutterRate'))
                this.ClutterRate = config.ClutterRate;
            end
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
            
            % reset measLikelihood matrix
            this.pmeasLikelihood_ = [];
            
            % Prediction for birth components 
            BirthComponents = this.BirthIntFcn();
            
            this.PredComponents.Means = zeros(this.Model.Dyn.NumStateDims,...
                                              this.NumComponents);
            this.PredComponents.Covars = zeros(this.Model.Dyn.NumStateDims,...
                                               this.Model.Dyn.NumStateDims,...
                                               this.NumComponents);
            this.PredComponents.Weights = zeros(1, this.NumComponents);
                                          
            % Prediction for existing components
            for i = 1:this.NumComponents
                
                % Use underlying filter to predict components
                this.Filter.StateMean = this.Components.Means(:,i);
                this.Filter.StateCovar = this.Components.Covars(:,:,i);
                this.Filter.predict();
                
                % Assing predicted state to components
                this.PredComponents.Means(:,i) = this.Filter.PredStateMean;
                this.PredComponents.Covar(:,:,i) = this.Filter.PredStateCovar;
            end
            
            % Compute predicted weights for components
            this.PredComponents.Weights = (1-this.ProbOfDeath)*this.Components.Weights;
            
            % Add birth component to prediction
            this.PredComponents.Means = [this.PredComponents.Means, BirthComponents.Means];
            this.PredComponents.Covars = [this.PredComponents.Covars, BirthComponents.Covars];
            this.PredComponents.Weights = [this.PredComponents.Weights, BirthComponents.Weights];
            
            this.NumPredComponents = this.NumComponents + size(BirthComponents.Means,2);
        end
        

        function update(this)
        % UPDATE Perform SMC PHD update step
        %
        % Usage
        % -----
        % * update(this) calculates the corrected intensity.
        %
        % See also UnscentedParticleFilterX, predict, smooth.
        
            this.NumMeasurements = size(this.MeasurementList,2);
            
            numUpdateComponents = (1+this.NumMeasurements)*this.NumPredComponents;
            updateComponents.Means = zeros(this.Model.Dyn.NumStateDims,...
                                           numUpdateComponents);
            updateComponents.Covars = zeros(this.Model.Dyn.NumStateDims,...
                                            this.Model.Dyn.NumStateDims,...
                                            numUpdateComponents);
            updateComponents.Weights = zeros(1, numUpdateComponents);
            updateComponents.PredMeasMeans = zeros(this.Model.Obs.NumObsDims,...
                                                   numUpdateComponents);
            updateComponents.InnovErrCovars = zeros(this.Model.Obs.NumObsDims,...
                                                    this.Model.Obs.NumObsDims,...
                                                    numUpdateComponents);
            updateComponents.CrossCovars = zeros(this.Model.Dyn.NumStateDims,...
                                                 this.Model.Obs.NumObsDims,...
                                                 numUpdateComponents);
            
            % Create null measurement hypothesis components
            updateComponents.Means(:,1:this.NumPredComponents) = this.PredComponents.Means;
            updateComponents.Covar(:,:,1:this.NumPredComponents) = this.PredComponents.Covars;
            updateComponents.Weights(1:this.NumPredComponents) = (1-this.ProbOfDetection)*this.PredComponents.Weights;
            
            likelihoods = zeros(1,numUpdateComponents);
            % Create true measurement hypothesis components
            for i = 1:this.NumPredComponents
                for j = 1:this.NumMeasurements
                    
                    % Use underlying filter to update components
                    this.Filter.resetStateEstimates();
                    this.Filter.Measurement = this.MeasurementList(:,j);
                    this.Filter.PredStateMean = this.PredComponents.Means(:,i);
                    this.Filter.PredStateCovar = this.PredComponents.Covars(:,:,i);
                    this.Filter.update();

                    % Assing predicted state to components
                    compIndex = this.NumPredComponents + (i-1)*this.NumMeasurements + j;
                    updateComponents.Means(:,compIndex) = this.Filter.PredStateMean;
                    updateComponents.Covars(:,:,compIndex) = this.Filter.PredStateCovar;
                    updateComponents.PredMeasMeans(:,compIndex) = this.Filter.PredMeasMean;
                    updateComponents.InnovErrCovars(:,:,compIndex) = this.Filter.InnovErrCovar;
                    updateComponents.CrossCovars(:,:,compIndex) = this.Filter.CrossCovar;
                end
                offsetInd = this.NumPredComponents + (i-1)*this.NumMeasurements;
                likelihoods(offsetInd:offsetInd+this.NumMeasurements) = ...
                    gauss_pdf(this.MeasurementList,...
                              updateComponents.PredMeasMeans(:,compIndex),...
                              this.Filter.InnovErrCovar);
            end
             
            updateComponents.Weights(this.NumPredComponents+1:end) = ...
                this.ProbOfDetection*this.PredComponents.Weights(i)...
                *gauss_pdf(this.MeasurementList(:,j),this.Filter.PredMeas,this.Filter.InnovErrCovar);
            denom = denom + updateComponents{compIndex}.Weight;
            
            % Normalise weights
            for i = 1:numel(updateComponents)
                updateComponents{i}.Weight = updateComponents{i}.Weight/denom;
            end
            
            % Prune
            for i = 1:numel(updateComponents)
                if(updateComponents{i}.Weight<=this.PruneThreshold)
                    updateComponents{i} = [];
                end
                
            end
            updateComponents = updateComponents(~cellfun('isempty',updateComponents));
            
            % Merge 
            
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