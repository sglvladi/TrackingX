classdef ConstantVelocityModelX <  DynamicModelX % Handle class with copy functionality
    % ConstantVelocityModelX class
    %
    % Summary of ConstantVelocityModel
    % This is a class implementation of a linear-Gaussian Constant Velocity Dynamic Model. 
    %
    % The model is described by the following SDEs:
    %
    %   dx   = v_x*dt                    | Position on X-axis (m)
    %   dy   = v_y*dt                    | Position on Y axis (m)
    %   dv_x = q*dW_t,  W_t~N(0,q^2)     | Speed on X-axis    (m/s)
    %   dv_y = q*dB_t,  B_t~N(0,q^2)     | Speed on Y-axis    (m/s)
    %
    % ConstantVelocityModelX Properties:
    %    - Params   = structure with fields:
    %       .xDim            = dimensionality, 3 possible settings: 2-4-6 => (x,v_x)-(x,y,v_x,v_y)-(x,y,z,v_x,v_y,v_z), default 4)
    %       .q               = Process noise diffusion coefficient (default 0.01 m/s^2)
    %       .f(Dt)           = Linear time-variant process transition function handle f(x_k|x_{k-1}), returns matrix (Defined internally, but can be overloaded)
    %    	.Q(Dt)           = Time-variant process noise covariance function handle, returns (nx x nx) matrix (Defined internally, but can be overloaded)
    %       .smartArgs       = Set (true|false) to (enable|disable) 'name',value pair input parsing.
    %
    % ConstantVelocityModelX Methods:
    %    sys        - State vector process transition function f(x_{k-1},w_k) 
    %    sys_cov    - State covariance process Q_k (Only REQUIRED for Gaussian Filters)
    %    sys_noise  - Process noise sample generator 
    %    eval       - Evaluates the probability p(x_k|x_{k-1}) = N(x_k; x_{k-1}, Q) of a set of new states, given a set of (particle) state vectors  
  
    properties
    end
    
    methods
        function this = ConstantVelocityModelX(Init)
        % CONSTANTVELOCITYMODELX Constructor method
        %   
        % INPUTS:    Init   - Structure with the following fields:
        %            .xDim   - State dimensionality, 3 possible settings: 
        %                      2-4-6 => (x,v_x)-(x,y,v_x,v_y)-(x,y,z,v_x,v_y,v_z)       
        %            .q      - Process noise diffusion coefficient 
        %                      (Optional, Default = 0.01 m/s^2)
        %
        % OUTPUTS:   this   - Object instance   
        %
        % USAGE:     Init.xDim = 4;
        %            Init.q    = 0.01;
        %            cv = ConstantVelocityModelX(Init)
        %
        %  See also sys, sys_cov, sys_noise, eval.   
                   
            % Add .xDim
            Params.xDim = Init.xDim;
            
            % Validate .q
            if isfield(Init, 'q')
                Params.q = Init.q;
            else
                Params.q = 0.01;
                fprintf('[CVModel] Process noise diffusion coefficient missing... Applying default setting "Params.q = 0.01"..\n');
            end
            
            % Define .F
            switch(Params.xDim)
                case(2)
                    Params.F = @(Dt,~) [1 Dt;
                                        0 1]; 

                case(4)
                    Params.F = @(Dt,~) [1 0 Dt 0;
                                        0 1 0 Dt;
                                        0 0 1 0;
                                        0 0 0 1];
                case(6)
                    Params.F = @(Dt,~) [1 0 0 Dt 0 0;
                                        0 1 0 0 Dt 0;
                                        0 0 1 0 0 Dt;
                                        0 0 0 1 0 0;
                                        0 0 0 0 1 0;
                                        0 0 0 0 0 1];
                otherwise
                    error('[CVModel] Invalid model dimensionality (Params.xDim)... Valid options are 2,4 and 6 ..\n');  
            end
            
            
            % Validate .Q
            switch(Params.xDim)
                case(2)
                    Params.Q = @(Dt) [Dt^3/3 Dt^2/2;
                                      Dt^2/2 Dt]; 

                case(4)
                    Params.Q = @(Dt) [Dt^3/3, 0, Dt^2/2, 0;
                                      0, Dt^3/3, 0, Dt^2/2; 
                                      Dt^2/2, 0, Dt, 0;
                                      0, Dt^2/2, 0, Dt]*Params.q;
                case(6)
                    Params.Q = @(Dt) [Dt^3/3 0 0 Dt^2/2 0 0;
                                      0 Dt^3/3 0 0 Dt^2/2 0;
                                      0 0 Dt^3/3 0 0 Dt^2/2;
                                      Dt^2/2 0 0 Dt 0 0;
                                      0 Dt^2/2 0 0 Dt 0;
                                      0 0 Dt^2/2 0 0 Dt]*Params.q;
            end
            
            % Call SuperClass method
            this@DynamicModelX(Params);
        end
        
        function xk = sys(this, k, xkm1, wk)
        % SYS State vector process transition function f(x_{k-1},w_k) 
        %
        % INPUTS:   k      - Time variable 
        %           xkm1   - a (xDim x Ns) matrix of Ns state vectors from
        %                    time k-1. (Optional, Default = 1)
        %           wk     - a (xDim x Ns) matrix Ns of process noise 
        %                    vectors, corresponding to the state vectors xkm1. 
        %                    (Optional, Default = 0)  
        %
        % OUTPUTS:  xk     - a (xDim x Ns) matrix of Ns state vectors, which
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
                case 1 
                    k = 1;
                    xkm1 = 1;
                    wk   = 0;
                case 2
                    xkm1 = 1;
                    wk   = 0;
                case 3
                    wk   = 0;
            end
            
            % Compute result
            xk = this.Params.F(k)*xkm1 + wk;
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
        % See also CONSTANTVELOCITYMODELX, SYS, SYS_NOISE, EVAL.
            
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
        % See also CONSTANTVELOCITYMODELX, SYS, SYS_COV, EVAL.
       
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
        % See also CONSTANTVELOCITYMODELX, SYS, SYS_COV, SYS_NOISE.
        
            xk_km1 = this.sys(k,xkm1);
            ProbabilityMatrix = zeros(size(xk,2), size(xkm1,2));
            if(size(xkm1,2)>size(xk,2))
                for i=1:size(xk,2)
                    ProbabilityMatrix(i,:) = gauss_pdf(xk(:,i), xk_km1, this.Params.Q(k));
                end
            else
                for i=1:size(xkm1,2)
                    ProbabilityMatrix(:,i) = gauss_pdf(xk, xk_km1(:,i), this.Params.Q(k))';  
                end
             end
                        
        end
        
    end
end