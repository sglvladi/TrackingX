classdef ParticleFilterX < FilterX
% ParticleFilterX class
%
% Summary of ParticleFilterX:
% This is a class implementation of a SIR Particle Filter. (Alg. 4 of [1])
%
% ParticleFilterX Properties: (*)
%   + NumParticles      The number of particles employed by the Particle Filter
%   + Particles         A (NumStateDims x NumParticles) matrix used to store 
%                       the last computed/set filtered particles  
%   + Weights           A (1 x NumParticles) vector used to store the weights
%                       of the last computed/set filtered particles
%   + PredParticles     A (NumStateDims x NumParticles) matrix used to store 
%                       the last computed/set predicted particles  
%   + PredWeights       A (1 x NumParticles) vector used to store the weights
%                       of the last computed/set predicted particles
%   + Measurement       A (NumObsDims x 1) matrix used to store the received measurement
%   + ControlInput      A (NumCtrDims x 1) matrix used to store the last received 
%                       control input
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
%   ¬ StateMean         A (NumStateDims x 1) vector used to store the last 
%                       computed filtered state mean.  
%   ¬ StateCovar        A (NumStateDims x NumStateDims) matrix used to store
%                       the last computed filtered state covariance
%   ¬ PredStateMean     A (NumStateDims x 1) vector used to store the last 
%                       computed prediicted state mean  
%   ¬ PredStateCovar    A (NumStateDims x NumStateDims) matrix used to store
%                       the last computed/set predicted state covariance
%   ¬ PredMeasMean      A (NumObsDims x 1) vector used to store the last 
%                       computed predicted measurement mean
%   ¬ InnovErrCovar     A (NumObsDims x NumObsDims) matrix used to store the
%                       last computed innovation error covariance
%   ¬ CrossCovar        A (NumStateDims x NumObsDims) matrix used to store 
%                       the last computed cross-covariance Cov(X,Y)  
%   + Model             An object handle to StateSpaceModelX object
%       + Dyn = Object handle to DynamicModelX SubClass      
%       + Obs = Object handle to ObservationModelX SubClass 
%       + Ctr = Object handle to ControlModelX SubClass    
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
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
    properties
        NumParticles = 1000
        Particles
        Weights
        PredParticles
        PredWeights
        ControlInput
        ResamplingScheme = 'Systematic'
        ResamplingPolicy = {'TimeInterval',1}
        Resampler = SystematicResamplerX();
        PriorDistFcn
    end
    
    properties (Access=protected)
        pmeasLikelihood_ = []
    end
    
    properties (Dependent)
        StateMean
        StateCovar
        PredStateMean
        PredStateCovar
        PredMeasMean
        InnovErrCovar
        CrossCovar
        MeasLikelihood
    end
    
    methods
        function this = ParticleFilterX(varargin)
        % PARTICLEFILTERX Constructor method
        %   
        % DESCRIPTION: 
        % * pf = ParticleFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * pf = ParticleFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * pf = ParticleFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * pf = ParticleFilterX(ssm,priorDistFcn) returns an object handle, 
        %   preconfigured with the provided StateSpaceModel object handle ssm
        %   and the prior information about the state, provided in the form 
        %   of the priorDistFcn function.
        % * pf = ParticleFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * NumParticles        (Scalar) The number of particles to be employed by the  
        %                       Particle Filter. [Default = 1000]
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
        %                       return 2 outputs. If a
        %                       PriorDistFcn is specified, then any values provided for the
        %                       PriorParticles and PriorWeights arguments are ignored.
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
        %
        %  See also predict, update, smooth. 
                 
            % Call SuperClass method
            this@FilterX(varargin{:});
            
            if(nargin==0)
                this.Resampler = SystematicResamplerX();
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'NumParticles'))
                        this.NumParticles  = config.NumParticles;
                    end
                    if (isfield(config,'PriorDistFcn'))
                        [this.Particles,this.Weights] = ...
                            config.PriorDistFcn(this.NumParticles);
                    elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                         this.Particles = config.PriorParticles;
                         this.Weights = config.PriotWeights;
                    end
                     if (isfield(config,'Resampler'))
                         this.Resampler = config.Resampler;
                     elseif (isfield(config,'ResamplingScheme'))
                         if(strcmp(config.ResamplingScheme,'Systematic'))
                             this.Resampler = SystematicResamplerX();
                         elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                             this.Resampler = MultinomialResamplerX();
                         end
                     else
                         this.Resampler = SystematicResamplerX();
                     end
                     if (isfield(config,'ResamplingPolicy'))
                         this.ResamplingPolicy = config.ResamplingPolicy;
                     end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            if (isfield(config,'NumParticles'))
                this.NumParticles  = config.NumParticles;
            end
            if (isfield(config,'PriorDistFcn'))
                [this.Particles,this.Weights] = ...
                    config.PriorDistFcn(this.NumParticles);
            elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                 this.Particles = config.PriorParticles;
                 this.Weights = config.PriotWeights;
            end
            if (isfield(config,'Resampler'))
                 this.Resampler = config.Resampler;
            elseif (isfield(config,'ResamplingScheme'))
                if(strcmp(config.ResamplingScheme,'Systematic'))
                     this.Resampler = SystematicResamplerX();
                 elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                     this.Resampler = MultinomialResamplerX();
                end
            else
                 this.Resampler = SystematicResamplerX();
            end
            if (isfield(config,'ResamplingPolicy'))
                 this.ResamplingPolicy = config.ResamplingPolicy;
            end 
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the Particle Filter with a certain set of
        % parameters.  
        %   
        % DESCRIPTION: 
        % * initialise(pf,ssm) initialises the ParticleFilterX object kf
        %   with the provided StateSpaceModelX object ssm.
        % * initialise(pf,priorParticles,priorWeights)initialises the 
        %   ParticleFilterX object pf with the provided StateSpaceModel object    
        %   ssm and the prior information about the state, provided in the form  
        %   of the priorParticles and priorWeights variables.
        % * initialise(pf,ssm,priorDistFcn) initialises the ParticleFilterX
        %   object pf with the provided StateSpaceModel object handle ssm
        %   and the prior information about the state, provided in the form 
        %   of the priorDistFcn function.
        % * initialise(pf,___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * NumParticles        (Scalar) The number of particles to be employed by the  
        %                       Particle Filter. [Default = 1000]
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
        %                       return 2 outputs. If a
        %                       PriorDistFcn is specified, then any values provided for the
        %                       PriorParticles and PriorWeights arguments are ignored.
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
        %
        %  See also predict, update, smooth. 
                    
            initialise@FilterX(this);
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    if (isfield(config,'Model'))
                        this.Model = config.Model;
                    end
                    if (isfield(config,'NumParticles'))
                        this.NumParticles  = config.NumParticles;
                    end
                    if (isfield(config,'PriorDistFcn'))
                        [this.Particles,this.Weights] = ...
                            config.PriorDistFcn(this.NumParticles);
                    elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                         this.Particles = config.PriorParticles;
                         this.Weights = config.PriotWeights;
                    end
                     if (isfield(config,'Resampler'))
                         this.Resampler = config.Resampler;
                     elseif (isfield(config,'ResamplingScheme'))
                         if(strcmp(config.ResamplingScheme,'Systematic'))
                             this.Resampler = SystematicResamplerX();
                         elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                             this.Resampler = MultinomialResamplerX();
                         end
                     else
                         this.Resampler = SystematicResamplerX();
                     end
                     if (isfield(config,'ResamplingPolicy'))
                         this.ResamplingPolicy = config.ResamplingPolicy;
                     end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            if (isfield(config,'Model'))
                this.Model = config.Model;
            end
            if (isfield(config,'NumParticles'))
                this.NumParticles  = config.NumParticles;
            end
            if (isfield(config,'PriorDistFcn'))
                [this.Particles,this.Weights] = ...
                    config.PriorDistFcn(this.NumParticles);
            elseif ((isfield(config,'priorParticles'))&&(isfield(config,'priorParticles')))
                 this.Particles = config.PriorParticles;
                 this.Weights = config.PriotWeights;
            end
            if (isfield(config,'Resampler'))
                 this.Resampler = config.Resampler;
            elseif (isfield(config,'ResamplingScheme'))
                if(strcmp(config.ResamplingScheme,'Systematic'))
                     this.Resampler = SystematicResamplerX();
                 elseif(strcmp(config.ResamplingScheme,'Multinomial'))
                     this.Resampler = MultinomialResamplerX();
                end
            else
                 this.Resampler = SystematicResamplerX();
            end
            if (isfield(config,'ResamplingPolicy'))
                 this.ResamplingPolicy = config.ResamplingPolicy;
            end 
        end
        
        function predict(this)
        % PREDICT Perform SIR Particle Filter prediction step
        %   
        % DESCRIPTION: 
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % MORE DETAILS:
        % * ParticleFilterX uses the Model class property, which should be an
        %   instance of the TrackingX.Models.StateSpaceModel class, in order
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
            
            % Extract model parameters
            f  = @(x,wk) this.Model.Dyn.feval(x, wk); % Transition function
            wk = this.Model.Dyn.random(this.NumParticles); % Process noise
        
            % Propagate particles through the dynamic model
            this.PredParticles = ParticleFilterX_Predict(f,this.Particles,wk);
            this.PredWeights = this.Weights;
            
            predict@FilterX(this);
        end
        
        function update(this)
        % UPDATE Perform SIR Particle Filter update step
        %   
        % DESCRIPTION: 
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        %   See also ParticleFilterX, predict, smooth.
        
            if(size(this.Measurement,2)>1)
                error('[PF] More than one measurement have been provided for update. Use ParticleFilterX.UpdateMulti() function instead!');
            elseif size(this.Measurement,2)==0
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
            end
            
            % Perform update
            this.Weights = ...
                ParticleFilterX_UpdateWeights(@(y,x) this.Model.Obs.pdf(y,x),...
                        this.Measurement,this.PredParticles,this.PredWeights);
                    
            if(strcmp(this.ResamplingPolicy(1),"TimeInterval"))
                [this.Particles,this.Weights] = ...
                    this.Resampler.resample(this.PredParticles,this.Weights);
            end
            
            update@FilterX(this);
        end
        
        function updatePDA(this, assocWeights, MeasLikelihood)
        % UPDATEPDA - Performs PF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % DESCRIPTION:
        %  * updatePDA(assocWeights) Performs KF-PDA update step for multiple 
        %    measurements based on the provided (1-by-Nm+1) association weights 
        %    matrix assocWeights.
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        %   See also KalmanFilterX, Predict, Iterate, Smooth, resample.
        
            NumMeasurements = size(this.Measurement,2);  
            if(~NumMeasurements)
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
                this.Particles = this.PredParticles;
                this.Weights = this.PredWeights;
                return;
            end
            
            if(~exist('assocWeights','var'))
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            if(nargin<3)
                MeasLikelihood = this.MeasLikelihood;  
            end
            
            % Perform update
            this.Weights = ...
                ParticleFilterX_UpdatePDA(@(y,x) this.Model.Obs.heval(y,x),this.Measurement,this.PredParticles,...
                                          this.PredWeights, assocWeights, MeasLikelihood);
            
            if(strcmp(this.ResamplingPolicy(1),"TimeInterval"))
                [this.Particles,this.Weights] = ...
                    this.Resampler.resample(this.PredParticles,this.Weights);
            end
        end
                
        function smoothed_estimates = Smooth(this, filtered_estimates)
        % Smooth - Performs FBS smoothing on a provided set of estimates
        %           (Based on [1])
        %   
        %   Inputs:
        %       filtered_estimates: a (1 x N) cell array, where N is the total filter iterations and each cell is a copy of this.Params after each iteration
        %   
        %   Outputs:
        %       smoothed_estimates: a copy of the input (1 x N) cell array filtered_estimates, where the .x and .P fields have been replaced with the smoothed estimates   
        %
        %   (Virtual inputs at each iteration)        
        %           -> filtered_estimates{k}.particles          : Filtered state mean estimate at timestep k
        %           -> filtered_estimates{k}.P          : Filtered state covariance estimate at each timestep
        %           -> filtered_estimates{k+1}.x_pred   : Predicted state at timestep k+1
        %           -> filtered_estimates{k+1}.P_pred   : Predicted covariance at timestep k+1
        %           -> smoothed_estimates{k+1}.x        : Smoothed state mean estimate at timestep k+1
        %           -> smoothed_estimates{k+1}.P        : Smoothed state covariance estimate at timestep k+1 
        %       where, smoothed_estimates{N} = filtered_estimates{N} on initialisation
        %
        %   (NOTE: The filtered_estimates array can be accumulated by running "filtered_estimates{k} = ukf.Params" after each iteration of the filter recursion) 
        %   
        %   Usage:
        %       ukf.Smooth(filtered_estimates);
        %
        %   [1] Mike Klaas, Mark Briers, Nando de Freitas, Arnaud Doucet, Simon Maskell, and Dustin Lang. 2006. Fast particle smoothing: if I had a million particles. In Proceedings of the 23rd international conference on Machine learning (ICML '06). ACM, New York, NY, USA, 481-488.
        %
        %   See also ParticleFilterX, Predict, Update, Iterate.
        
            % Allocate memory
            N                           = length(filtered_estimates);
            smoothed_estimates          = cell(1,N);
            smoothed_estimates{N}       = filtered_estimates{N}; 
            
            % Perform Rauch–Tung–Striebel Backward Recursion
            for k = N-1:-1:1
                lik = this.DynModel.eval(filtered_estimates{k}.k, filtered_estimates{k+1}.particles, filtered_estimates{k}.particles);
                denom = sum(filtered_estimates{k}.w(ones(this.Params.Np,1),:).*lik,2)'; % denom(1,j)
                smoothed_estimates{k}.w = filtered_estimates{k}.w(1,:) .* sum(smoothed_estimates{k+1}.w(ones(this.Params.Np,1),:).*lik'./denom(ones(this.Params.Np,1),:),2)';
                smoothed_estimates{k}.particles =  filtered_estimates{k}.particles;
                smoothed_estimates{k}.x = sum(smoothed_estimates{k}.w.*smoothed_estimates{k}.particles,2);
            end
        end        
        
        % ===============================>
        % ACCESS METHODS
        % ===============================>
        
        function StateMean = get.StateMean(this)
            StateMean = getStateMean(this);
        end
        
        function StateCovar = get.StateCovar(this)
            StateCovar = getStateCovar(this);
        end
        
        function PredStateMean = get.PredStateMean(this)
            PredStateMean = getPredStateMean(this);
        end
        
        function PredStateCovar = get.PredStateCovar(this)
            PredStateCovar = getPredStateCovar(this);
        end
        
        function PredMeasMean = get.PredMeasMean(this)
            PredMeasMean = getPredMeasMean(this);
        end
        
        function InnovErrCovar = get.InnovErrCovar(this)
            InnovErrCovar = getInnovErrCovar(this);
        end
        
        function MeasLikelihood = get.MeasLikelihood(this)
            MeasLikelihood = getMeasLikelihood(this);
        end

%         function set.Measurement(this,newMeasurement)
%             setMeasurement(this,newMeasurement);
%         end
    end
    
    methods (Access = protected)
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        function StateMean = getStateMean(this)
            StateMean = sum(this.Weights.*this.Particles,2);
        end
        
        function StateCovar = getStateCovar(this)
            StateCovar = weightedcov(this.Particles,this.Weights);
        end
        
        function PredStateMean = getPredStateMean(this)
            PredStateMean = sum(this.PredWeights.*this.PredParticles,2);
        end
        
        function PredStateCovar = getPredStateCovar(this)
            PredStateCovar = weightedcov(this.PredParticles,this.PredWeights);
        end
        
        function PredMeasMean = getPredMeasMean(this)
            PredMeasMean = this.Model.Obs.heval(this.PredStateMean);
        end
        
        function InnovErrCovar = getInnovErrCovar(this)
            transParticles = this.Model.Obs.heval(this.PredParticles,true);
            InnovErrCovar = weightedcov(transParticles,this.PredWeights);
        end
        
        function MeasLikelihood = getMeasLikelihood(this)
            if(isempty(this.pmeasLikelihood_))
                this.pmeasLikelihood_ = this.Model.Obs.pdf(this.Measurement,this.PredParticles);               
            end
            MeasLikelihood = this.pmeasLikelihood_;
        end
        
%         function setMeasurement(this,newMeasurement)
%             this.Measurement = newMeasurement;
%         end
    end
end