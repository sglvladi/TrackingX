classdef ConstantHeadingModelX <  DynamicModelX % Handle class with copy functionality
    % ConstantHeadingModelX class
    %
    % Summary of ConstantHeadingModelX
    % This is a class implementation of a 2D nonlinear-Gaussian Constant Heading Dynamic Model [1].
    % The model is described by the following SDEs:
    %
    %   dx = v*cos(h)*dt                      | Position on X-axis (m)
    %   dy = v*sin(h)*dt                      | Position on Y axis (m)
    %   dv = q_vel*dW_t,  W_t~N(0,q_vel^2)    | Absolute velocity  (m/s)
    %   dh = q_head*dB_t, B_t~N(0,q_head^2)   | Heading            (rad)
    %
    % [1] P. A. Kountouriotis and S. Maskell, "Maneuvering target tracking using an unbiased nearly constant heading model," 2012 15th International Conference on Information Fusion, Singapore, 2012, pp. 2249-2255.
    %
    % ConstantHeadingModelX Properties:
    %    - Params   = structure with fields:
    %       .xDim            = 4 (x,y,v,h), should not be changed!
    %    	.q_vel           = Velocity noise standard deviation. Default = 0.01 m/s^2
    %    	.q_head          = Heading noise standard deviation. Default = 0.16 rad/s
    %    	.f(k, xk)        = State-time variant process transition function handle f(x_k|x_{k-1}), returns matrix (Defined internally, but can be overloaded)
    %    	.Q(k)            = Time variant process noise covariance function handle, returns (.xDim x .xDim) matrix (Defined internally, but can be overloaded)
    %
    % ConstantHeadingModelX Methods:
    %    sys        - State vector process transition function f(x_{k-1},w_k) 
    %    sys_cov    - State covariance process Q_k (Only REQUIRED for Gaussian Filters)
    %    sys_noise  - Process noise sample generator 
    %    eval       - Evaluates the probability p(x_k|x_{k-1}) = N(x_k; x_{k-1}, Q) of a set of new states, given a set of (particle) state vectors  
  
    properties
    end
    
    methods
        function this = ConstantHeadingModelX(Init)
        % CONSTANTHEADINGMODELX Constructor method
        %   
        % INPUTS:    Init   - Structure with the following fields:     
        %            .q_vel  - Velocity noise standard deviation. Default = 0.01 m/s^2
        %            .q_head - Heading noise standard deviation. Default = 0.16 rad/s           
        %
        % OUTPUTS:   this   - Object instance   
        %
        % USAGE:     Init.q_vel  = 0.01;
        %            Init.q_head = 0.13;
        %            cv = ConstantHeadingModelX(Init)
        % 
        % See also sys, sys_cov, sys_noise, eval.      
                        
            % Validate .q_vel
            if isfield(Init, 'q_vel')
                Params.q_vel = Init.q_vel;
            else
                Params.q_vel = 0.01;
                fprintf('[CHModel] Velocity noise standard deviation missing... Applying default setting "q_vel = 0.01"..\n');
            end
            
            % Validate .q_head
            if isfield(Init, 'q_head')
                Params.q_head = Init.q_head;
            else
                Params.q_head = 0.01;
                fprintf('[CHModel] Velocity noise standard deviation missing... Applying default setting "q_vel = 0.01"..\n');
            end
            
            % Define .xDim
            Params.xDim = 4;
            
            % Define .f
            Params.f = @(k, xkm1) [xkm1(1,:)+k*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+k*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:); xkm1(4,:)];
            
            % Define .Q
            Params.Q = @(k) k*diag([Params.q_vel^4, Params.q_vel^4, Params.q_vel^2, Params.q_head^2]);
                
            this@DynamicModelX(Params);
      
        end
        
        function xk = sys(this, k, xkm1, wk)
        % SYS State vector process transition function f(x_{k-1},w_k) 
        %
        % INPUTS:   k      - Time variable 
        %           xkm1   - A (xDim x Ns) matrix of Ns state vectors from
        %                    time k-1.
        %           wk     - A (xDim x Ns) matrix Ns of process noise 
        %                    vectors, corresponding to the state vectors xkm1. 
        %                    (Optional, Default = 0)  
        %
        % OUTPUTS:  xk     - A (xDim x Ns) matrix of Ns state vectors, which
        %                    have been propagated through the dynamic model
        %        
        % USAGE:
        %     - xk = DynModel.sys(k)
        %     - xk = DynModel.sys(k,xkm1)
        %     - xk = DynModel.sys(k,xkm1,wk)
        %     Also use xk = DynModel.sys(k,xkm1,DynModel.sys_noise(k,Ns)) to propagate Ns samples with process noise
        %
        % See also CONSTANTVELOCITYMODELX, SYS_COV, SYS_NOISE, EVAL.
        
            switch(nargin)
                case 3
                    wk   = 0;
            end
            
            % Compute result
            xk = this.Params.f(k,xkm1) + wk;
        end
        
        function Qk = sys_cov(this, k)
        % SYS_COV Returns process covariance Q_k. Only applies to Gaussian 
        %         models.
        %
        % INPUTS:   k   - Time variable
        %
        % OUTPUS:   Qk  - The (xDim x xDim) process covariance matrix.  
        %
        % USAGE:    Qk = DynModel.sys_cov(k); 
        %           Qk = DynModel.sys_cov('k',k);
        %
        % See also CONSTANTHEADINGMODELX, SYS, SYS_NOISE, EVAL.
            
            % Return process covariance
            Qk = this.Params.Q(k); % Time variant
        end
        
        function wk = sys_noise(this, k, Ns)
        % SYS_NOISE Process noise sample generator 
        %
        % INPUTS:   k   - Time variable
        %           Ns  - The number samples to be generated 
        %                 (Optional, default is 1)  
        %
        % OUTPUTS:  wk  - A (xDim x Ns) matrix of Ns process noise samples   
        %
        % USAGE:    wk = DynModel.sys_noise(k) 
        %           wk = DynModel.sys_noise(k, Ns) 
        %
        % See also CONSTANTHEADINGMODELX, SYS, SYS_COV, EVAL.
       
            switch(nargin)
                case 2
                    Ns = 1;
            end
            wk = mvnrnd(zeros(this.Params.xDim,Ns)',this.Params.Q(k))';
        end
        
        function ProbabilityMatrix = eval(this, k, xk, xkm1)
        % EVAL Evaluates the probability p(x_k|x_{k-1}) of a set of new 
        %      states, given a set of (particle) state vectors  
        % 
        % INPUTS:   k     - Time variable
        %           xk    - A (xDim x Np) matrix of Np new state vectors.
        %           xkm1  - A (xDim x Ns) matrix of Ns old state vectors.
        %
        % OUTPUTS:  ProbabilityMatrix - A (Np x Ns) matrix of probabilities 
        %                               p(x_k|x_{k-1})    
        %
        % USAGE:    ProbabilityMatrix = DynModel.eval(this, Dt, x_k, x_km1) 
        %
        % See also CONSTANTHEADINGMODELX, SYS, SYS_COV, SYS_NOISE.
        
            xk_km1 = this.sys(k,xkm1);
            ProbabilityMatrix = gauss_pdf(xk_km1, xk, this.Params.Q(k));
        end
    end
end