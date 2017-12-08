classdef PositionalObsModelX < ObservationModelX
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
        function this = PositionalObsModelX(Init)
        % PSOTIONALOBSMODELX Constructor method
        %   
        % INPUTS:   Init   - Structure with the following fields: 
        %           .xDim  - State dimensionality, 3 possible settings: 
        %                    2-4-6 => (x,v_x)-(x,y,v_x,v_y)-(x,y,z,v_x,v_y,v_z)
        %                    (Optional, Default = 4)
        %           .yDim  - Observation dimensionality, 3 possible settings:
        %                    1-2-3 => (x)-(x,y)-(x,y,z)
        %                    (Optional, Default = 2)
        %           .r     - Process noise diffusion coefficient 
        %                    (Optional, Default = 0.01 m/s^2)                   
        %
        % USAGE:    po = PositionalObsModelX(Init)
        %
        % See also OBS, OBS_COV, OBS_NOISE, EVAL.   
            
            % Validate .xDim
            if isfield(Init,'xDim')
                Params.xDim = Init.xDim;
            else
                Params.xDim = 4;
            end
            
            % Validate .yDim
            if isfield(Init,'yDim')
                Params.yDim = Init.yDim;
            else
                Params.yDim = 2;
            end
            
            % Validate .r
            if isfield(Init, 'r')
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
                    Params.H = @(~) [[1 0; 0 1],zeros(2,Params.xDim-Params.yDim)]; 
                case(3)
                    Params.H = @(~) [[1 0 0; 0 1 0; 0 0 1], zeros(3,Params.xDim-Params.yDim)]; 
            end
            
            % Define .h_inv
            Params.H_inv = @(~) 1; 
            
            % Validate .r
            if ~isfield(Params,'r')
                Params.r = 0.1;
                fprintf('[POModel] Observation noise standard deviation missing... Applying default setting "r = 0.1"..\n');
            end
            
            % Define .R
            if ~isfield(Params,'R')
                switch(Params.yDim)
                    case(1)
                        Params.R = @(~) Params.r^2; 
                    case(2)
                        Params.R = @(~) eye(2)*Params.r^2;
                    case(3)
                        Params.R = @(~) eye(3)*Params.r^2;
                end
            end
            
            this@ObservationModelX(Params);
      
        end
        
        function yk = obs(this, k, xk, vk)
        % OBS State-space to measurement-space state projection function h(k,xk,vk)
        %
        % INPUTS:   k 	- Time variable (!NOT USED HERE!)
        %           xk  - A (xDim x Ns) matrix of Ns state vectors from time k 
        %                 (Optional)
        %           vk  - A (yDim x Ns) matrix Ns of observation noise vectors, 
        %                 corresponding to the state vectors xk. 
        %                 (Optional)  
        %
        % OUTPUTS:  yk  - A (yDim x Ns) matrix of Ns measurement vectors, 
        %                 which have been produced by projecting the predicted
        %                 state vectors into measurement space
        %        
        % USAGE:    yk = ObsModel.obs(k,xk,vk);
        %
        %           With a properly defined ObsModel.obs_noise(k,Ns) function,
        %           you can also use it to generate vk internally and project 
        %           a number of samples, with noise, as so:
        %
        %           yk = ObsModel.obs(k,xk,ObsModel.obs_noise(k,size(xk,2)));
        %
        %   See also OBS_COV, OBS_NOISE, SAMPLE, EVAL.
            
            switch(nargin)
                case 1
                    k = 1;
                    xk = 1;
                    vk = 0;
                case 2
                    xk = 1;
                    vk = 0;
                case 3
                    vk = 0;
            end
            
            % Compute predicted measurement
            yk = this.Params.H(k)*xk + vk ;           
        end
        
        function xk = obs_inv(this, k, yk, vk)
        % OBS_INV Measurement-space to state-space state projection function h_inv(k,xk,vk)
        %
        % INPUTS:   k 	- Time variable (!NOT USED HERE!)
        %           yk  - A (yDim x Ns) matrix of Ns measurement vectors from time k 
        %                 (Optional)
        %           vk  - A (yDim x Ns) matrix Ns of observation noise vectors, 
        %                 corresponding to the state vectors xk. 
        %                 (Optional)  
        %
        % OUTPUTS:  xk  - A (xDim x Ns) matrix of Ns state vectors, 
        %                 which have been produced by projecting the 
        %                 measurement vectors into state space
        %        
        % USAGE:    xk = ObsModel.obs_inv(k,yk,vk);
        %
        %           With a properly defined ObsModel.obs_noise(k,Ns) function,
        %           you can also use it to generate vk internally and project 
        %           a number of samples, with noise, as so:
        %
        %           xk = ObsModel.obs_inv(k,yk,ObsModel.obs_noise(k,Ns));
        %
        %   See also OBS_COV, OBS_NOISE, SAMPLE, EVAL.
            
            switch(nargin)
                case 1
                    k = 1;
                    yk = 1;
                    vk = 0;
                case 2
                    yk = 1;
                    vk = 0;
                case 3
                    vk = 0;
            end
            
            % Compute predicted measurement
            xk = this.Params.H_inv(k)*yk - vk ;           
        end
        
        function Rk = obs_cov(this, k)
        % OBS_COV Returns measurement covariance Rk. Only applies to Gaussian models..
        %
        % INPUTS:   k   - Time variable 
        %
        % OUTPUTS:  Rk  - A (yDim x yDim) process covariance matrix.  
        %
        % USAGE:    Rk = ObsModel.obs_cov(); % For time invariant process covariance
        %           Rk = ObsModel.obs_cov(k); % For time variant process covariance
        %
        % See also OBS, OBS_NOISE, SAMPLE, EVAL.
                    
            % Return process covariance
            Rk = this.Params.R();    % Time invariant
        end
        
        function vk = obs_noise(this, k, Ns)
        % OBS_NOISE Process noise sample generator 
        %
        % INPUTS:   k  - Time variable 
        %           Ns - The number samples to be generated 
        %                (Optional, Default = 1)  
        %
        % OUTPUTS:  vk - A (yDim x Ns) matrix of Ns process noise samples   
        %
        % USAGE:    vk = ObsModel.obs_noise()
        %           vk = ObsModel.obs_noise(k) % k can have any scalar numeric value, or be [] 
        %           vk = ObsModel.obs_noise(k, Ns) % k can have any scalar numeric value, or be [] 
        %
        % See also OBS, OBS_COV, SAMPLE, EVAL.
            
            if nargin~=3
                Ns = 1;
            end
        
            % Generate noise samples
            vk = mvnrnd(zeros(this.Params.yDim,Ns)',this.Params.R())';
        end
        
         function LikelihoodMatrix = eval(this, k, yk, xk)
        % EVAL Evaluates the probability p(y_k|x_k) of a set of new states, 
        % given a set of (particle) state vectors  
        % 
        % INPUTS:   k   - Time variable
        %           yk  - A (yDim x Np) matrix of Np new state vectors
        %           xk  - A (xDim x Ns) matrix of Ns old state vectors
        %   All inputs are REQUIRED. Set k=[] for time invariant systems
        %
        % OUTPUTS:  LikelihoodMatrix - A (Np x Ns) matrix of probabilities 
        %                              p(y_k|x_k)    
        %
        % USAGE:    LikelihoodMatrix = ObsModel.eval(k, yk, xk)
        %
        % See also OBS, OBS_COV, OBS_NOISE, SAMPLE.
        
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
        
        function samples = sample(this, k, mu, Ns)
        % SAMPLE Sample from measurement model, given a mean and number of 
        % samples  
        %
        % INPUTS:   k   - Time variable (NOT USED HERE)
        %           mu  - A (yDim x 1) mean vector
        %           Ns  - Number of samples (Optional, Default = 1)
        %
        % OUTPUTS:  samples - A (yDim x Ns) matrix of samples, drawn using 
        %                     the observation model   
        %
        % USAGE:    samples = sample(this, mu, Ns)
        %
        % See also OBS, OBS_COV, OBS_NOISE, EVAL.
            if(~exist('Ns', 'var'))
                Ns = 1;
            end
            samples = mvnrnd(mu', this.Params.R(), Ns)';
        end
    end
end