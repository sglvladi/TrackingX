classdef UnscentedParticleFilterX < ParticleFilterX
% UnscentedParticleFilterX class
%
% Summary of UnscentedParticleFilterX:
% This is a class implementation of a SIR Particle Filter. (Alg. 4 of [1])
%
% UnscentedParticleFilterX Properties: (*)
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
%   + Alpha             ||
%   + Kappa             || UKF scaling parameters, as described in [1]
%   + Beta              || 
%
%   (*) NumStateDims, NumObsDims and NumCtrDims denote the dimentionality of 
%       the state, measurement and control vectors respectively.
%
% UnscentedParticleFilterX Methods:
%   + UnscentedParticleFilterX  - Constructor method
%   + predict        - Performs UKF prediction step
%   + update         - Performs UKF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% [1] E. A. Wan and R. Van Der Merwe, "The unscented Kalman filter for nonlinear estimation," 
%     Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and 
%     Control Symposium (Cat. No.00EX373), Lake Louise, Alta., 2000, pp. 153-158.
%
% See also ParticleFilterX, ExtendedParticleFilterX

    properties %(Access = private)
        ukf   % Instance of an UnscentedKalmanFilterX
    end
    
    methods
        function this = UnscentedParticleFilterX(varargin)
        % UNSCENTEDPARTICLEFILTERX Constructor method
        %   
        % DESCRIPTION: 
        % * epf = UnscentedParticleFilterX() returns an unconfigured object 
        %   handle. Note that the object will need to be configured at a 
        %   later instance before any call is made to it's methods.
        % * epf = UnscentedParticleFilterX(ssm) returns an object handle,
        %   preconfigured with the provided StateSpaceModelX object handle ssm.
        % * epf = UnscentedParticleFilterX(ssm,priorParticles,priorWeights) 
        %   returns an object handle, preconfigured with the provided  
        %   StateSpaceModel object handle ssm and the prior information   
        %   about the state, provided in the form of the priorParticles 
        %   and priorWeights variables.
        % * epf = UnscentedParticleFilterX(ssm,priorDistFcn) returns an object
        %   handle, preconfigured with the provided StateSpaceModel object 
        %   handle ssm and the prior information about the state, provided  
        %   in the form of the priorDistFcn function.
        % * epf = UnscentedParticleFilterX(___,Name,Value,___) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        % INPUT ARGUMENTS:
        % * NumParticles        (Scalar) The number of particles to be employed by the  
        %                       Unscented Particle Filter. [Default = 1000]
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
        % * Alpha,Kappa,Beta    (Scalar) values used to parameterise the
        %                       Unscented Kalman Filter. Alpha, Kappa and
        %                       Beta should be passed individually.
        %
        %  See also predict, update, smooth. 
        
           % Call SuperClass method
           this@ParticleFilterX(varargin{:});
           
           % Instantiate EKF
           this.ukf = UnscentedKalmanFilterX(varargin{:});
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the Unscented Particle Filter with a certain 
        % set of parameters.  
        %   
        % DESCRIPTION: 
        % * initialise(epf,ssm) initialises the UnscentedParticleFilterX object 
        %   epf with the provided StateSpaceModelX object ssm.
        % * initialise(epf,priorParticles,priorWeights)initialises the 
        %   UnscentedParticleFilterX object pf with the provided StateSpaceModel     
        %   object ssm and the prior information about the state, provided in  
        %   the form  of the priorParticles and priorWeights variables.
        % * initialise(epf,ssm,priorDistFcn) initialises the UnscentedParticleFilterX
        %   object pf with the provided StateSpaceModel object handle ssm
        %   and the prior information about the state, provided in the form 
        %   of the priorDistFcn function.
        % * initialise(epf,___,Name,Value,___) instantiates an object handle, 
        %   configured with the options specified by one or more Name,Value
        %   pair arguments.
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
        % * Alpha,Kappa,Beta    (Scalar) values used to parameterise the Unscented Kalman Filter.
        %                       Alpha, Kappa and Beta should be passed individually.
        %
        %  See also predict, update, smooth. 
                    
           % Call SuperClass method
           initialise@ParticleFilterX(this,varargin{:});
           
           % Initialise Unscented Kalman Filter
           this.ukf.initialise(varargin{:});
        end
        
        function predict(this)
        % PREDICT Perform Unscented Particle Filter prediction step
        %   
        % DESCRIPTION: 
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % MORE DETAILS:
        % * UnscentedParticleFilterX uses the Model class property, which should 
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
                                
            % Compute UKF prior mean and covariance
            this.ukf.StateMean  = this.StateMean;
            this.ukf.StateCovar = this.StateCovar;
            
            % Iterate UKF to obtain Optimal Proposal
            this.ukf.predict();
            
            predict@FilterX(this);
        end
        
        function update(this)
        % UPDATE Perform Unscented Particle Filter update step
        %   
        % DESCRIPTION: 
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        %   See also UnscentedParticleFilterX, predict, smooth.
            
            % Update UKF measurement
            this.ukf.Measurement = this.Measurement;
            
            % Perform UKF update to obtain Optimal Proposal
            this.ukf.update();
            
            % Sample from UKF proposal
            this.PredParticles = mvnrnd(this.ukf.StateMean', this.ukf.StateCovar,this.NumParticles)'; 
            this.PredWeights = this.Weights;
            
            % Call SuperClass method
            update@ParticleFilterX(this);             
        end
        
        function UpdatePDA(this, assocWeights, LikelihoodMatrix)
        % UpdatePDA - Performs UPF update step, for multiple measurements
        %   
        %   Inputs:
        %       assoc_weights: a (1 x Nm+1) association weights matrix. The first index corresponds to the dummy measurement and
        %                       indices (2:Nm+1) correspond to
        %                       measurements. Default = [0, ones(1,ObsNum)/ObsNum];
        %       LikelihoodMatrix: a (Nm x Np) likelihood matrix, where Nm is the number of measurements and Np is the number of particles.
        %
        %   (NOTE: The measurement "this.Params.y" needs to be updated, when necessary, before calling this method) 
        %   
        %   Usage:
        %       (pf.Params.y = y_new; % y_new is the new measurement)
        %       upf.UpdateMulti(assocWeights, LikelihoodMatrix); 
        %
        %   See also UParticleFilterX, Predict, Iterate, Smooth, resample.
        
            % Update EKF measurement
            this.ukf.Params.y = this.Params.y;
            
            % Perform EKF update to obtain Optimal Proposal
            this.ukf.UpdatePDA(assocWeights);
            
            % Sample from EKF proposal
            this.Params.particles = mvnrnd(this.ukf.Params.x', this.ukf.Params.P,this.Params.Np)'; % Sample from optimal proposal
            
            if(exist('LikelihoodMatrix','var'))
                UpdatePDA@ParticleFilterX(this, assocWeights, LikelihoodMatrix); 
            else
                UpdatePDA@ParticleFilterX(this, assocWeights);
            end
        end
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
            PredStateMean = this.ukf.StateMean;
        end
        
        function PredStateCovar = getPredStateCovar(this)
            PredStateCovar = this.ukf.StateCovar;
        end
        
        function PredMeasMean = getPredMeasMean(this)
            PredMeasMean = this.Model.Obs.heval(this.PredStateMean);
        end
        
        function InnovErrCovar = getInnovErrCovar(this)
            InnovErrCovar = this.ukf.InnovErrCovar;
        end
        
        function setModel(this,newModel)
            this.Model = newModel;
            this.ukf.Model = newModel;
        end
        
        function setMeasurement(this,newMeasurement)
            this.Measurement = newMeasurement;
            this.ukf.Measurement = newMeasurement;
        end
    end
end