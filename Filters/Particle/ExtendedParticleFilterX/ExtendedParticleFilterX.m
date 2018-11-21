classdef ExtendedParticleFilterX < ParticleFilterX
% ExtendedParticleFilterX class
%
% Summary of ExtendedParticleFilterX:
% This is a class implementation of a Particle Filter that utilises an Extended
% Kalman Filter to generate its importance/proposal density.
%
% ExtendedParticleFilterX Properties: (*)
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
%       + Dyn - Object handle to TransitionModelX SubClass      
%       + Obs - Object handle to MeasurementModelX SubClass 
%       + Ctr - Object handle to ControlModelX SubClass    
%
%   (*) NumStateDims, NumObsDims and NumCtrDims denote the dimentionality of 
%       the state, measurement and control vectors respectively.
%
% ExtendedParticleFilterX Methods:
%   + ExtendedParticleFilterX  - Constructor method
%   + predict        - Performs UKF prediction step
%   + update         - Performs UKF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes

    properties %(Access = private)
        ekf   % Instance of an ExtendedKalmanFilterX
    end
    
    methods
        function this = ExtendedParticleFilterX(varargin)
        % EXTENDEDPARTICLEFILTERX Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % Prior: struct, optional
        %   A (NumStateDims x 1) column vector, representing the prior
        %   state mean, which is copied over to Posterior.
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
           this@ParticleFilterX(varargin{:});
           
           % Instantiate EKF
           this.ekf = ExtendedKalmanFilterX(varargin{:});
        end
        
        function initialise(this,varargin)
        % INITIALISE Initialise the Extended Particle Filter with a certain 
        % set of parameters.  
        %   
        % Parameters
        % ----------       
        % Model: StateSpaceModelX
        %   An object handle to StateSpaceModelX object.
        % Prior: struct, optional
        %   A (NumStateDims x 1) column vector, representing the prior
        %   state mean, which is copied over to Posterior.
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
                    
           % Call SuperClass method
           initialise@ParticleFilterX(this,varargin{:});
           
           % Initialise Extended Kalman Filter
           this.ekf.initialise(varargin{:});
        end
        
        function predict(this)
        % PREDICT Perform Extended Particle Filter prediction step
        %   
        % Usage
        % -----
        % * predict(this) calculates the predicted system state and measurement,
        %   as well as their associated uncertainty covariances.
        %
        % More details
        % ------------
        % * ExtendedParticleFilterX uses the Model class property, which should 
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
        
            % Compute EKF prior mean and covariance
            this.ekf.StatePosterior = this.StatePosterior;
            
            % Iterate EKF to obtain Optimal Proposal
            this.StatePrediction = this.ekf.predict();
   
        end
        
       function update(this)
        % UPDATE Perform Extended Particle Filter update step
        %   
        % Usage
        % -----
        % * update(this) calculates the corrected sytem state and the 
        %   associated uncertainty covariance.
        %
        % See also ExtendedParticleFilterX, predict, smooth.
                        
            % Perform EKF update to obtain Optimal Proposal
            this.ekf.update();
            
            % Sample from EKF proposal
            this.StatePrediction = ParticleStateX(this.ekf.StatePosterior.Mean, this.ekf.StatePosterior.Covar, this.NumParticles);
            
            % Call SuperClass method
            update@ParticleFilterX(this);             
        end
        
        function updatePDA(this, assocWeights, MeasLikelihood)
        % UPDATEPDA - Performs EPF update step, for multiple measurements
        %             Update is performed according to the generic (J)PDAF equations [1] 
        % 
        % Usage
        % -----
        %  * updatePDA(assocWeights) Performs EPF-PDA update step for multiple 
        %    measurements based on the provided (1-by-Nm+1) association weights 
        %    matrix assocWeights.
        %
        %   [1] Y. Bar-Shalom, F. Daum and J. Huang, "The probabilistic data association filter," in IEEE Control Models, vol. 29, no. 6, pp. 82-100, Dec. 2009.
        %
        % See also KalmanFilterX, Predict, Iterate, Smooth, resample.
            
            NumMeasurements = size(this.Measurement,2);  
            if(~NumMeasurements)
                warning('[PF] No measurements have been supplied to update track! Skipping Update step...');
                this.Particles = this.Model.Dyn.feval(this.Particles,true);
                %this.Weights = this.PredWeights;
                return;
            end
            
            if(~exist('assocWeights','var'))
                assocWeights = [0, ones(1,ObsNum)/ObsNum]; % (1 x Nm+1)
            end
            
            % Perform EKF update to obtain Optimal Proposal
            this.ekf.updatePDA(assocWeights);
            
            % Sample from EKF proposal
            this.PredParticles = mvnrnd(this.ekf.StateMean', this.ekf.StateCovar,this.NumParticles)'; 
            this.PredWeights = this.Weights;
            
            if(exist('MeasLikelihood','var'))
                updatePDA@ParticleFilterX(this, assocWeights, MeasLikelihood); 
            else
                updatePDA@ParticleFilterX(this, assocWeights);
            end
        end
    end
    
    methods (Access = protected)
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        function Model = setModel(this,newModel)
            Model = newModel;
            this.ekf.Model = newModel;
        end
        
        function MeasurementList = setMeasurementList(this,newMeasurementList)
            MeasurementList = newMeasurementList;
            this.ekf.MeasurementList = newMeasurementList;
        end
    end
end