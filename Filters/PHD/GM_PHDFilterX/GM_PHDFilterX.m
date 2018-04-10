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
%   + Components - A structure used to store the last filtered component with the 
%                  following fields (updated on every this.update() call):
%                   + Means: (NumStateDims x NumComponents) matrix 
%                       Used to store the last computed/set filtered component means
%                   + Covars: (NumStateDims x NumStateDims x NumComponents) matrix
%                       Used to store the last computed/set filtered component covariances.
%                   + Weights: (1 x NumComponents) row vector
%                       Used to store the last computed/set filtered component weights.
%   + NumPredComponents - The number of the last predicted GM components employed by the 
%                         PHD Filter
%   + PredComponents - A structure used to store the last predicted components:
%                      with the following fields (updated on every this.predict() call):
%                   + Means: (NumStateDims x NumPredComponents) matrix 
%                       Used to store the last computed/set predicted component means
%                   + Covars: (NumStateDims x NumStateDims x NumPredComponents) matrix
%                       Used to store the last computed/set predicted component covariances.
%                   + Weights: (1 x NumPredComponents) row vector
%                   + PredMeasMeans: (NumObsDims x NumPredComponents) matrix
%                       Used to store the last computed/set predicted component 
%                       measurement means.
%                   + InnovErrCovars: (NumObsDims x NumObsDims x NumPredComponents) matrix
%                       Used to store the last computed/set predicted component 
%                       measurement (innovation) covariances.
%                   + CrossCovars: (NumStateDims x NumObsDims x NumPredComponents) matrix
%                       Used to store the last computed/set predicted component 
%                       state-measurement cross covariances.
%   + MeasurementList - A (NumObsDims x NumObs) matrix used to store the received 
%                       measurements
%   + BirthIntFcn - A function handle, BirthComponents = BirthIntFcn(), 
%                   which when called generates a set of birth components,
%                   in the form of a BirthComponents structure with the
%                   following fields:
%                       + Means: (NumStateDims x NumBirthComponents) matrix 
%                           Used to store the birth component means
%                       + Covars: (NumStateDims x NumStateDims x NumBirthComponents) matrix
%                           Used to store the birth component covariances.
%                       + Weights: (1 x NumBirthComponents) row vector 
%                           Used to store the birth component weights.
%   + Model - An object handle to StateSpaceModelX object
%       + Dyn - Object handle to DynamicModelX SubClass      
%       + Obs - Object handle to ObservationModelX SubClass 
%       + Ctr - Object handle to ControlModelX SubClass     
%   + Filter - A KalmanFilterX (sub)class instance, which is used to
%              predict and update all the GM components.
%   + Gater - A GaterX sub-class, used to perform gating between GM
%             components and measurements. (Optional, default = None)
%   + ProbOfDeath - The probability that a target may cease to exist
%                   between consecutive iterations of the filter.
%   + ProbOfDetection - The probablity that a target will be detected in
%                       a given measurement scan.
%   + ClutterRate - The mean number of clutter measurements per unit volume
%                   of the surveillance region, under the assumption of uniformly
%                   distributed clutter, with Poisson distributed number of clutter
%                   measurtements
%   + MaxNumComponents - The maximum allowable number of compunents following 
%                        an update() call 
%   + MergeThreshold - Minimum distance between the means of components. If
%                      the distance between the means of two (or more) components
%                      is less that MergeThreshold, the the Components will
%                      be merged during the next update() call.
%   + PruneThreshold - Minimum weight for components. Any components whose
%                      weight is below PruneThreshold, will be removed on
%                      the next update() call.
%   + NumTargets - The estimated number of targets following an update step.
 
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
 
    properties (Access = private, Hidden)
         
    end
    properties (Dependent)
        NumComponents
        NumMeasurements
        NumTargets
    end
    properties
        Components
        PredComponents
        NumPredComponents
        MaxNumComponents = 100
        PruneThreshold = 1e-5
        MergeThreshold = 0.5
        Gater
        BirthIntFcn 
        ProbOfDeath  
        ProbOfDetection 
        ClutterRate 
        MeasLikelihood
        MeasWeights = 1;
        WeightsPerMeasurement      
        MeasurementList
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
                        this.Gater = EllipsoidalGaterX('ProbOfGating',config.ProbOfGating,'NumObsDims',2);
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
                    if (isfield(config,'Filter'))
                        this.Filter = config.Filter;
                    else
                        this.Filter = KalmanFilterX(this.Model);
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
            if (isfield(config,'Filter'))
                this.Filter = config.Filter;
            else
                this.Filter = KalmanFilterX(this.Model);
            end
        end
         
        function initialise(this,varargin)
        % INITIALISE Initialise the GM PHD Filter with a certain 
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
        % * initialise(phd,ssm) initialises the GM_PHDFilterX object 
        %   phd with the provided StateSpaceModelX object ssm.
        % * initialise(phd,priorParticles,priorWeights)initialises the 
        %   GM_PHDFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(phd,ssm,priorDistFcn) initialises the GM_PHDFilterX
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
                    if (isfield(config,'Filter'))
                        this.Filter = config.Filter;
                    else
                        this.Filter = KalmanFilterX(this.Model);
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
            if (isfield(config,'Filter'))
                this.Filter = config.Filter;
            else
                this.Filter = KalmanFilterX(this.Model);
            end
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
             
            % reset measLikelihood matrix
            this.pmeasLikelihood_ = [];
             
            % 1) Birth component extraction 
            % --------------------------
            BirthComponents = this.BirthIntFcn();
                                                     
            % 2) State Prediction for existing components
            % ----------------------------------------
            % Preallocate memory
            this.PredComponents.Means = zeros(this.Model.Dyn.NumStateDims,...
                                              this.NumComponents);
            this.PredComponents.Covars = zeros(this.Model.Dyn.NumStateDims,...
                                               this.Model.Dyn.NumStateDims,...
                                               this.NumComponents);
            this.PredComponents.Weights = zeros(1, this.NumComponents);
             
            % Perform prediction
            for i = 1:this.NumComponents
                 
                % Use underlying filter to predict components
                this.Filter.StateMean = this.Components.Means(:,i);
                this.Filter.StateCovar = this.Components.Covars(:,:,i);
                this.Filter.predictState();
                 
                % Assing predicted state to components
                this.PredComponents.Means(:,i) = this.Filter.PredStateMean;
                this.PredComponents.Covars(:,:,i) = this.Filter.PredStateCovar;
            end
             
            % Compute predicted weights for components
            this.PredComponents.Weights = (1-this.ProbOfDeath)*this.Components.Weights;
             
            % Add birth component to prediction
            this.PredComponents.Means = [this.PredComponents.Means, BirthComponents.Means];
            this.PredComponents.Covars = cat(3,this.PredComponents.Covars, BirthComponents.Covars);
            this.PredComponents.Weights = [this.PredComponents.Weights, BirthComponents.Weights];
             
            % Compute number of prediction components
            this.NumPredComponents = this.NumComponents + size(BirthComponents.Means,2);
                                             
            % 3) Measurement prediction
            % ----------------------
            % Preallocate memory
            this.PredComponents.PredMeasMeans = zeros(this.Model.Obs.NumObsDims,...
                                                      this.NumPredComponents);
            this.PredComponents.InnovErrCovars = zeros(this.Model.Obs.NumObsDims,...
                                                       this.Model.Obs.NumObsDims,...
                                                       this.NumPredComponents);
            this.PredComponents.CrossCovars = zeros(this.Model.Dyn.NumStateDims,...
                                                    this.Model.Obs.NumObsDims,...
                                                    this.NumPredComponents);
            % Perform measurement prediction
            for i = 1:this.NumPredComponents
                     
                % Use underlying filter to predict measurement
                this.Filter.resetStateEstimates();
                this.Filter.PredStateMean = this.PredComponents.Means(:,i);
                this.Filter.PredStateCovar = this.PredComponents.Covars(:,:,i);
                this.Filter.predictObs();
                 
                % Assigned measurement prediction to component
                this.PredComponents.PredMeasMeans(:,i) = this.Filter.PredMeasMean;
                this.PredComponents.InnovErrCovars(:,:,i) = this.Filter.InnovErrCovar;
                this.PredComponents.CrossCovars(:,:,i) = this.Filter.CrossCovar;
            end
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
            numUpdateComponents = (1+this.NumMeasurements)*this.NumPredComponents;
            % Preallocate memory
            updateComponents.Means = zeros(this.Model.Dyn.NumStateDims,...
                                           numUpdateComponents);
            updateComponents.Covars = zeros(this.Model.Dyn.NumStateDims,...
                                            this.Model.Dyn.NumStateDims,...
                                            numUpdateComponents);
            updateComponents.Weights = zeros(1, numUpdateComponents);
             
            % Create null measurement hypothesis components
            updateComponents.Means(:,1:this.NumPredComponents) = this.PredComponents.Means;
            updateComponents.Covars(:,:,1:this.NumPredComponents) = this.PredComponents.Covars;
            updateComponents.Weights(1:this.NumPredComponents) = (1-this.ProbOfDetection)*this.PredComponents.Weights;
             
            % Compute measurement likelihoods
            likelihoods = this.MeasLikelihood;
             
            % Compute normalising constant
            Ck = zeros(1,this.NumMeasurements);
            for i = 1:this.NumMeasurements   % for all measurements
                Ck(i) = sum(this.ProbOfDetection*likelihoods(i,:).*this.PredComponents.Weights,2);
            end
             
            % Create true measurement hypothesis components
            for i = 1:this.NumPredComponents
                 
                % Use underlying filter to update components
                this.Filter.resetStateEstimates();
                this.Filter.Measurement = this.MeasurementList;
                this.Filter.PredStateMean = this.PredComponents.Means(:,i);
                this.Filter.PredStateCovar = this.PredComponents.Covars(:,:,i);
                this.Filter.PredMeasMean = this.PredComponents.PredMeasMeans(:,i);
                this.Filter.InnovErrCovar = this.PredComponents.InnovErrCovars(:,:,i);
                this.Filter.CrossCovar = this.PredComponents.CrossCovars(:,:,i);
                this.Filter.update();
                     
                for j = 1:this.NumMeasurements
                    % Assing predicted state to components
                    compIndex = this.NumPredComponents + (i-1)*this.NumMeasurements + j;
                    updateComponents.Means(:,compIndex) = this.Filter.StateMean(:,j);
                    updateComponents.Covars(:,:,compIndex) = this.Filter.StateCovar;
                    % Compute corresponding weights
                    updateComponents.Weights(compIndex) = ...
                        this.ProbOfDetection*this.PredComponents.Weights(i)...
                        *likelihoods(j,i)/(this.ClutterRate+Ck(j));
                end
            end  
             
            % 5) Prune
            % --------
            [updateComponents.Means,updateComponents.Covars,updateComponents.Weights] = ...
                gauss_prune(updateComponents.Means,updateComponents.Covars,...
                            updateComponents.Weights, this.PruneThreshold);
             
            % 6) Merge
            % --------
            [updateComponents.Means,updateComponents.Covars,updateComponents.Weights] = ...
               gauss_merge(updateComponents.Means,updateComponents.Covars,...
                           updateComponents.Weights, this.MergeThreshold);
                         
            % 7) Cap
            % ------
            [this.Components.Means,this.Components.Covars,this.Components.Weights] = ...
                gauss_cap(updateComponents.Means,updateComponents.Covars,...
                            updateComponents.Weights, this.MaxNumComponents);
        end
     
        function MeasLikelihood = get.MeasLikelihood(this)
            MeasLikelihood = getMeasLikelihood(this);
        end
         
        function NumComponents = get.NumComponents(this)
            NumComponents = length(this.Components.Weights);
        end
         
        function NumTargets = get.NumTargets(this)
            NumTargets = sum(this.Components.Weights);
        end
         
        function NumMeasurements = get.NumMeasurements(this)
            NumMeasurements = size(this.MeasurementList,2);
        end
    end
     
    methods (Access = protected)
         
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
         
        function MeasLikelihood = getMeasLikelihood(this)
            if(isempty(this.pmeasLikelihood_))
                this.pmeasLikelihood_ = this.Model.Obs.pdf(this.MeasurementList, this.PredComponents.Means, this.PredComponents.Covars);               
            end
            MeasLikelihood = this.pmeasLikelihood_;
        end
    end
end