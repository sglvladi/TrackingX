classdef ConstantHeadingModel2X <  DynamicModelX % Handle class with copy functionality
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
    %       .smartArgs       = Set (true|false) to (enable|disable) 'name',value pair input parsing.
    %
    % ConstantHeadingModelX Methods:
    %    sys        - State vector process transition function f(x_{k-1},w_k) 
    %    sys_cov    - State covariance process Q_k (Only REQUIRED for Gaussian Filters)
    %    sys_noise  - Process noise sample generator 
    %    eval       - Evaluates the probability p(x_k|x_{k-1}) = N(x_k; x_{k-1}, Q) of a set of new states, given a set of (particle) state vectors  
  
    properties
    end
    
    methods
        function this = ConstantHeadingModel2X(varargin)
        % ConstantHeadingModelX - Constructor method
        %   
        %   Inputs:
        %       Required
        %       ========
        %       None      
        %               
        %       Optional
        %       ========
        %       q_vel     = Velocity noise standard deviation. Default = 0.01 m/s^2
        %    	q_head    = Heading noise standard deviation. Default = 0.16 rad/s           
        %       smartArgs = Set (true|false) to (enable|disable) 'name',value pair input parsing. (Default = false)
        %                   The class methods support ('name',value) pair input arguments, however this can have a great impact on speed. 
        %                   To ensure optimum speed, such input parsing can be turned OFF by setting smartArgs to false.
        %
        %   Usage:
        %       cv = ConstantVelocityModelX(xDim,q,smartArgs)
        %       cv = ConstantVelocityModelX('xDim',xDim,'q',q,'smartArgs',smartArgs)
        %
        %   See also sys, sys_cov, sys_noise, eval.      
            
            % Define input parser
            p = inputParser;
            
            % Add q_vel as a parameter
            default_q_vel = [];
            validate_q_vel = @(x) isnumeric(x)&&isscalar(x);
            addOptional(p, 'q_vel', default_q_vel, validate_q_vel);
            
            % Add q_head as a parameter
            default_q_head = [];
            validate_q_head = @(x) isnumeric(x)&&isscalar(x);
            addOptional(p, 'q_head', default_q_head, validate_q_head);
            
            % Add smartArgs as a parameter
            default_sa = [];
            validate_sa = @(x) islogical(x);
            addOptional(p, 'smartArgs', default_sa, validate_sa)
            
            parse(p, varargin{:});
            
            % Validate .q_vel
            if ~isempty(p.Results.q_vel)
                Params.q_vel = p.Results.q_vel;
            else
                Params.q_vel = 0.01;
                fprintf('[CHModel] Velocity noise standard deviation missing... Applying default setting "q_vel = 0.01"..\n');
            end
            
            % Validate .q_head
            if ~isempty(p.Results.q_head)
                Params.q_head = p.Results.q_head;
            else
                Params.q_head = 0.01;
                fprintf('[CHModel] Velocity noise standard deviation missing... Applying default setting "q_vel = 0.01"..\n');
            end
            
            % Define .xDim
            Params.xDim = 4;
            
            % Validate .smartArgs
            if ~isempty(p.Results.smartArgs)
                Params.smartArgs = p.Results.smartArgs;
            else
                Params.smartArgs = true;
                warning('[CHModel] Smart (name,value) arguments are turned ON by default! Disable to improve speed.');
            end
            
            % Define .f
            Params.f = @(k, xkm1) [xkm1(1,:)+k*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+k*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:); xkm1(4,:)];
            
            % Define .Q
            Params.Q = @(k) k*diag([Params.q_vel^4, Params.q_vel^4, Params.q_vel^2, Params.q_head^2]);
                
            this@DynamicModelX(Params);
      
        end
        
        function xk = sys(this, k, varargin)
         % sys - State vector process transition function f(x_{k-1},w_k) 
        %
        %   Inputs:
        %       k      : (REQUIRED!) Time variable 
        %       xkm1   : (REQUIRED!) a (.xDim x Ns) matrix of Ns state vectors from time k-1, where nx is the dimensionality of the state (Optional)
        %       wk     : a (.xDim x Ns) matrix Ns of process noise vectors, corresponding to the state vectors xkm1. (Optional)  
        %
        %   The above arguments should be expected STRICTLY in the above given order, UNLESS arguments are passed as name,value pairs using the MATLAB inputParser scheme.
        %   In any case, variables should be provided in preceding order, e.g. if wk is supplied, this means both k and xkm1 must also be supplied by some means, even if the values are unnecessary/arbitrary
        %
        %   Outputs:
        %       xk: a (.xDim x Ns) matrix of Ns state vectors, which have been propagated through the dynamic model
        %        
        %   Usage:
        %     - xk = DynModel.sys(k) or xk = DynModel.sys('k', k) return xk = F(k)
        %     - xk = DynModel.sys(k,xkm1) or xk = DynModel.sys('k', k,'xkm1', xkm1) return xk = F(k)*xkm1
        %     - xk = DynModel.sys(k,xkm1,wk); or xk = DynModel.sys('k', k,'xkm1', xkm1, 'wk', wk); return xk = F(k)*xkm1 + wk
        %     Also use xk = DynModel.sys(k,xkm1,DynModel.sys_noise(k,Ns)) to propagate Ns samples with process noise
        %   See also sys_cov, sys_noise.
        
            if(this.Params.smartArgs)
                % Define input parser
                p = inputParser;

                % Define dafault values for xkm1 and wk
                default_xkm1 = 0;
                default_wk = 0;

                % Add k as an optional parameter
                validate_k = @(x) (isscalar(x) && isnumeric(x)); % Only valid either if k is numeric scalar
                addRequired(p, 'k', validate_k);
                % Add xkm1 as an optional parameter
                validate_xkm1 = @(x) (isempty(x) || (isnumeric(x)&&(size(x,1) == this.Params.xDim)));
                addOptional(p, 'xkm1', default_xkm1, validate_xkm1);
                % Add wk as an optional parameter
                validate_wk = @(x) (isempty(x) || (isnumeric(x)&&(size(x,1) == this.Params.xDim)));
                addOptional(p, 'wk', default_wk, validate_wk);
                parse(p,k,varargin{:});
                k    = p.Results.k;
                xkm1 = p.Results.xkm1;
                wk   = p.Results.wk;
            else
                switch(length(varargin))
                    case 0
                        xkm1 = 1;
                        wk   = 0;
                    case 1
                        xkm1 = varargin{1};
                        wk   = 0;
                    case 2
                        xkm1 = varargin{1};
                        wk   = varargin{2};
                    otherwise
                        error('Invalid number of arguments passed!');
                end
            end
            
            % Compute result
            xk = this.Params.f(k,xkm1) + wk;
        end
        
        function Qk = sys_cov(this, varargin)
        % sys_cov - Returns process covariance Q_k. Only applies to Gaussian models.
        %           Currently ProjectX filters do not support state-dependent process covariances. Extensions to state dependent covariance should be possible.
        %           Contact Lyudmil @ sglvladi@gmail.com, if you require such a functionality.
        %
        %   Inputs:
        %        k      : (REQUIRED!) Time variable 
        %
        %   Outputs:
        %       Qk: a (.xDim x .xDim) process covariance matrix.  
        %
        %   Usage:
        %   Qk = DynModel.sys_cov(k); 
        %   Qk = DynModel.sys_cov('k',k);
        %
        %   See also sys, sys_noise.
            
            if(this.Params.smartArgs)
                % Define input parser
                p = inputParser;

                % Add k as a REQUIRED argument
                validate_k = @(x) (~isempty(x) && isscalar(x) && isnumeric(x)); % Only valid if k is scallar and numeric
                addRequired(p, 'k', validate_k);
                parse(p,varargin{:});
                k = p.Results.k;
            else
                k = varargin{1};
            end
            
            % Return process covariance
            Qk = this.Params.Q(k); % Time variant
        end
        
        function w_k = sys_noise(this, varargin)
        % sys_noise - Process noise sample generator 
        %
        %   Inputs:
        %       k  : Time variable
        %       Ns : The number samples to be generated (Optional, default is 1)  
        %
        %   Outputs:
        %       w_k: a (.xDim x Ns) matrix of Ns process noise samples, where nx is the dimensionality of the state   
        %
        %   Usage:
        %   wk = DynModel.sys_noise() should ONLY be allowed for time invariant noise
        %   wk = DynModel.sys_noise(Dt) 
        %   wk = DynModel.sys_noise(Dt, Ns) 
        %
        %   See also sys, sys_cov.
        
            if(this.Params.smartArgs)
                % Define input parser
                p = inputParser;

                % Define default Ns
                default_Ns = 1;

                % Add k as a REQUIRED argument
                validate_k = @(x) (isscalar(x) && isnumeric(x)); % Only valid if k is scallar and numeric
                addRequired(p, 'k', validate_k);
                validate_Ns = @(x) (isscalar(x) && isnumeric(x) && x>0); % Only valid if k is scallar and numeric
                addOptional(p, 'Ns', default_Ns, validate_Ns);
                parse(p, varargin{:});
                k  = p.Results.k;
                Ns = p.Result.Ns;
            else
                switch(length(varargin))
                    case 1
                        k  = varargin{1};
                        Ns = 1;
                    case 3
                        k  = varargin{1};
                        Ns = varargin{2};
                        x_k = varargin{3};
                    otherwise
                        error("Invalid number of arguments passed!");
                end
            end
            w_k = mvnrnd(zeros(this.Params.xDim,Ns)',this.Params.Q(k))';
            w_k = [w_k(3,:)*k.*cos(w_k(3,:))/2-x_k(3,:).*w_k(4,:).*sin(x_k(4,:)); w_k(3,:)*k.*sin(w_k(3,:))/2+x_k(3,:).*w_k(4,:).*cos(x_k(4,:));w_k(3,:); w_k(4,:)];
        end
        
        function ProbabilityMatrix = eval(this, varargin)
        % eval - Evaluates the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors  
        % 
        %   Inputs:
        %       k     : Time variable
        %       x_k   : a (.xDim x Np) matrix of Np new state vectors, where nx is the dimensionality of the state
        %       x_km1 : a (.xDim x Ns) matrix of Ns old state vectors, where nx is the dimensionality of the state
        %   All inputs are REQUIRED. Set k=[] for time invariant systems
        %
        %   Outputs:
        %       ProbabilityMatrix: a (Np x Ns) matrix of probabilities p(x_k|x_{k-1})    
        %
        %   Usage:
        %   ProbabilityMatrix = eval(this, Dt, x_k, x_km1) Evaluates and produces a (Np x Ns) probability matrix.
        %
        %   See also obs, obs_cov, obs_noise, sample.
            
            if(this.Params.smartArgs)
                % Define input parser
                p = inputParser;

                % Add k as a required parameter
                validate_k = @(x) (isscalar(x) && isnumeric(x)); % Only valid either if k is empty (time invariant models), or if k is 
                addRequired(p, 'k', validate_k);
                % Add xkm1 as a required parameter
                validate_xk = @(x) (isnumeric(x)&&(size(x,1) == this.Params.xDim));
                addRequired(p, 'xk', validate_xk);
                % Add wk as a required parameter
                validate_xkm1 = @(x) (isnumeric(x)&&(size(x,1) == this.Params.xDim));
                addRequired(p, 'xkm1', validate_xkm1);
                parse(p, varargin{:});
                
                k    = p.Results.k;
                xk   = p.Results.xk;
                xkm1 = p.Results.xkm1;
            else
                k    = varargin{1};
                xk   = varargin{2};
                xkm1 = varargin{3};
            end
            
            ProbabilityMatrix = zeros(size(xk,2), size(xkm1,2));
            if(issymmetric(this.Params.Q(k)) && size(xkm1,2)>size(xk,2))
                % If R is symmetric and the number of state vectors is higher than the number of measurements,
                %   then evaluate N(x_km1; x_k, R) = N(x_k; x_km1, R) to increase speed
                for i=1:size(xk,2)
                    ProbabilityMatrix(i,:) = mvnpdf(xkm1', xk(:,i)', this.Params.Q(k))';
                end
            else
                for i=1:size(p.Results.xkm1,2)
                    ProbabilityMatrix(:,i) = mvnpdf(xk', xkm1(:,i)', this.Params.Q(k))';  
                end
             end
        end
    end
end