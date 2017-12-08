classdef GenericObservationModelX < ObservationModelX
    % PositionalObsModel class
    %
    % Summary of PositionalObsModel
    % This is a class implementation of a time & state invariant Linear-Gaussian Observation Model.
    %
    % ObservationModelX Interface:
    %   Properties:
    %       - Params = structure with the following compulsory fields:
    %           .xDim = Number of model/state dimensions 
    %           .yDim = Number of observation dimensions 
    %           .smartArgs       = Set (true|false) to (enable|disable) 'name',value pair input parsing.
    %
    %   Methods:
    %       obs                - State-space to measurement-space state projection function h(~,~,v_k)
    %       obs_cov            - Observation noise covariance R_k 
    %       obs_noise          - Observation noise sample generator
    %       sample             - Sample from measurement model, given a mean and number of samples
    %       eval_likelihood    - Evaluates the likelihood p(y_k|x_k)= N(y_k; x_k, R) of a set of measurements, given a set of (particle) states   

    properties
    end
    
    methods
        function this = GenericObservationModelX(Init)
        % PositionalObsModelX - Constructor method
        %   
        %   Inputs:
        %       Required
        %       ========        
        %               
        %       Optional
        %       ========
        %       xDim      = State dimensionality, 3 possible settings: 2-4-6 => (x,v_x)-(x,y,v_x,v_y)-(x,y,z,v_x,v_y,v_z), default 4)
        %       yDim      = Observation dimensionality, 3 possible settings: 2-4-6 => (x,v_x)-(x,y,v_x,v_y)-(x,y,z,v_x,v_y,v_z), default 2)
        %       r         = Process noise diffusion coefficient (default 0.01 m/s^2)           
        %       smartArgs = Set (true|false) to (enable|disable) 'name',value pair input parsing. (Default = false)
        %                   The class methods support ('name',value) pair input arguments, however this can have a great impact on speed. 
        %                   To ensure optimum speed, such input parsing can be turned OFF by setting smartArgs to false.
        %
        %   Usage:
        %       po = PositionalObsModelX(xDim,yDim,r,smartArgs)
        %       po = ConstantVelocityModelX('xDim',xDim,'yDim',yDim,'r',r,'smartArgs',smartArgs)
        %
        %   See also obs, obs_cov, obs_noise, eval.   
            
            % Define input parser
            
            % Add .xDim
            Params.xDim = Init.xDim;
            Params.yDim = Init.yDim;
            
            % Validate .r
            if isfield(Init,'r')
                Params.r = Init.r;
            else
                Params.r = 0.1;
                fprintf('[CVModel] Process noise diffusion coefficient missing... Applying default setting "Params.r = 0.1"..\n');
            end
            
            % Define .h
            switch(Params.yDim)
                case(1)
                    Params.H = @(~) [1 zeros(1,Params.xDim-Params.yDim)]; 
                case(2)
                    Params.H = @(~) [[0.5 0; 0 0.5],zeros(2,Params.xDim-Params.yDim)]; 
                case(3)
                    Params.H = @(~) [[1 0 0; 0 1 0; 0 0 1], zeros(3,Params.xDim-Params.yDim)]; 
            end
            
            % Define .R
            if ~isfield(Params,'R')
                switch(Params.yDim)
                    case(1)
                        Params.R = @(~) Params.r^2; 
                    case(2)
                        Params.R = @(~) eye(2)*(Params.r/2)^2;
                    case(3)
                        Params.R = @(~) eye(3)*Params.r^2;
                end
            end
            
            this@ObservationModelX(Params);
      
        end
        
        function yk = obs(this,k,xk,vk)
        % obs - State-space to measurement-space state projection function h(k,xk,vk)
        %
        %   Inputs:
        %       k      : Time variable (!NOT USED HERE!)
        %       xk     : a (.xDim x Ns) matrix of Ns state vectors from time k (Optional)
        %       vk     : a (.yDim x Ns) matrix Ns of observation noise vectors, corresponding to the state vectors xk. (Optional)  
        %
        %   The above arguments should be expected STRICTLY in the above given order, UNLESS arguments are passed as name,value pairs using the MATLAB inputParser scheme.
        %   In any case, variables should be provided in preceding order, e.g. if vk is supplied, this means both k and xk must also be supplied by some means, even if the values are unnecessary/arbitrary
        %
        %   Outputs:
        %       yk: a (.yDim x Ns) matrix of Ns measurement vectors, which have been produced by projecting the predicted state vectors into measurement space
        %        
        %   Usage:
        %     - yk = ObsModel.obs(k,xk,vk);
        %     or equivalently:
        %     - yk = ObsModel.obs('k', k, 'xk', xk, 'vk', vk); % if inputParser is used;
        %
        %     With a properly defined ObsModel.obs_noise(k,Ns) function, you can also use it to 
        %     generate vk internally and project a number of samples, with noise, as so:
        %     - yk = ObsModel.obs(k,xk,ObsModel.obs_noise(k,size(xk,2)));
        %
        %   See also obs_cov, obs_noise, sample, eval_likelihood.
            switch(nargin)
                case 1
                    k  = 0;
                    xk = 1;
                    vk = 0;
                case 2
                    xk = 1;
                    vk = 0;
                case 3
                    vk = 0;
            end
            yk = this.Params.H()*xk + vk ;           
        end
        
        function Rk = obs_cov(this,k)
        % obs_cov - Returns measurement covariance Rk. Only applies to Gaussian models.
        %           REQUIRED for your model to work with KalmanFilterX, EKalmanFilterX, UKalmanFilterX and all other Gaussian filters.
        %           Currently ProjectX filters do not support state-dependent measurement covariances. Extensions to state dependent covariance should be possible.
        %           Contact Lyudmil @ sglvladi@gmail.com, if you require such a functionality.
        %
        %   Inputs:
        %       k: Time variable (NOT USED HERE!) 
        %
        %   Outputs:
        %       Rk: a (.yDim x .yDim) process covariance matrix.  
        %
        %   Usage:
        %   Rk = ObsModel.obs_cov(); % For time invariant process covariance
        %   Rk = ObsModel.obs_cov(k); % For time variant process covariance
        %   Rk = ObsModel.obs_cov('k',k); % For time variant process covariance
        %
        %   See also obs, obs_noise.
                    
            % Return process covariance
            Rk = this.Params.R();    % Time invariant
        end
        
        function vk = obs_noise(this,k,Ns)
        % obs_noise - Process noise sample generator 
        %
        %   Inputs:
        %       k  : Time variable (!NOT USED HERE!)
        %       Ns : The number samples to be generated (Optional, default is 1)  
        %
        %   Outputs:
        %       vk: a (.yDim x Ns) matrix of Ns process noise samples, where nx is the dimensionality of the state   
        %
        %   Usage:
        %   vk = ObsModel.obs_noise()
        %   vk = ObsModel.obs_noise(k) % k can have any scalar numeric value, or be [] 
        %   vk = ObsModel.obs_noise(k, Ns) % k can have any scalar numeric value, or be [] 
        %
        %   See also obs, obs_cov.
        
            % Generate noise samples
            vk = mvnrnd(zeros(this.Params.yDim,Ns)',this.Params.R())';
        end
        
         function LikelihoodMatrix = eval(this,k,yk,xk)
        % eval - Evaluates the probability p(y_k|x_k) of a set of new states, given a set of (particle) state vectors  
        % 
        %   Inputs:
        %       k    : Time variable
        %       yk   : a (.yDim x Np) matrix of Np new state vectors, where nx is the dimensionality of the state
        %       xk   : a (.xDim x Ns) matrix of Ns old state vectors, where nx is the dimensionality of the state
        %   All inputs are REQUIRED. Set k=[] for time invariant systems
        %
        %   Outputs:
        %       LikelihoodMatrix: a (Np x Ns) matrix of probabilities p(y_k|x_k)    
        %
        %   Usage:
        %   LikelihoodMatrix = eval_likelihood(this, Dt, yk, xk) Evaluates and produces a (Np x Ns) likelihood matrix.
        %
        %   See also obs, obs_cov, obs_noise, sample.
        
            LikelihoodMatrix = zeros(size(yk,2), size(xk,2));
            if(size(xk,2)>size(yk,2))
                % If R is symmetric and the number of state vectors is higher than the number of measurements,
                %   then evaluate N(x_k; y_k, R) = N(y_k; x_k, R) to increase speed
                for i=1:size(yk,2)
                    LikelihoodMatrix(i,:) = mvnpdf(yk(:,i)', this.obs(k,xk)', this.Params.R())';
                end
            else
                for i=1:size(xk,2)
                    LikelihoodMatrix(:,i) = mvnpdf(yk', this.obs(k,xk(:,i))', this.Params.R())';  
                end
             end
                        
        end
        
        function samples = sample(this, ~, mu, Ns)
        % sample - Sample from measurement model, given a mean and number of samples  
        %
        %   Inputs:
        %       mu: a (ny x 1) mean vector, where ny is the dimensionality of the measurement
        %       Ns: number of samples (Optional, default = 1)
        %       ~Dt: time index/interval (!NOT USED HERE!)
        %
        %   Outputs:
        %       samples: a (ny x Ns) matrix of samples, drawn using the observation model   
        %
        %   Usage:
        %   samples = sample(this, mu, Ns) generates a (ny x Ns) set of Ns samples using the measurement model, given a mean mu.
        %
        %   See also obs, obs_cov, obs_noise, eval_likelihood.
            if(~exist('Ns', 'var'))
                Ns = 1;
            end
            samples = mvnrnd(mu', this.Params.R(), Ns)';
        end
    end
end