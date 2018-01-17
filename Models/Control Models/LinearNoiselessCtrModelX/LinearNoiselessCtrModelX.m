classdef LinearNoiselessCtrModelX < ControlModelX
    % LinearNoiselessCtrModelX class base class/template
    %
    % Summary of ControlModelX
    %  This is the class implementation of a noiseless control model, with 
    %  linear gain matrix B. 
    %
    % LinearNoiselessCtrModelX Interface:
    %   Base properties:
    %        Params     - structure with the following compulsory fields:
    %        .xDim  - Number of state vector dimensions
    %        .uDim  - Number of control input vector dimensions
    %        .B     - Linear control gain matrix
    %
    %   Base methods:
    %        ctr         - Control input vector process transition function b(u_k,w_k) 
    %        ctr_noise   - Control input noise sample generator function  
    %
    %
    % For examples see also the source of LINEARNOISELESSCTRMODELX
    
    properties
    end
    
    methods
        function this = LinearNoiselessCtrModelX(Init)
        % LINEARNOISELESSCTRMODELX Constructor method
        %   
        % INPUTS:   Virtually anything you want! 
        %   
        % USAGE:    CtrModel = LinearNoiselessCtrModelX(Init); 
        %
        %   You are free to do ANYTHING you want in here.
        %   Since no current TrackingX components auto-initialise any models, the number and type of constructor arguments can be set to anything that works for you. 
        %
        %   See also LINEARNOISELESSCTRMODELX.
        
            % Add .xDim
            Params.xDim = Init.xDim;
            
            % Add .uDim
            Params.uDim = Init.uDim;
            
            % Define .B
            Params.B = Init.B;
                        
            % Call SuperClass method
            this@ControlModelX(Params);     
        end
        
        function uk_x = ctr(this, k, uk, wk)
        % CTR Control input transition function b(u_{k},w_k) 
        %
        % INPUTS:   k    - Time variable
        %           uk   - A (uDim x 1) control input vector from time k
        %           wk   - A (uDim x 1) control input noise vector  
        %
        %           The above arguments should be expected STRICTLY in the 
        %           above given order, to ensure compatibility with all TrackingX
        %           components.
        %
        % OUTPUTS:  uk_x  - A (xDim x 1) state vector, produced by propagating
        %                   the control input through the control model. 
        %           
            switch(nargin)
                case 1 
                    k = 1;
                    uk = 1;
                case 2
                    uk = 1;
            end
            uk_x = this.Params.B()*uk;        
        end
        
        function Uk = ctr_cov(this, k)
        % CTR_COV Returns control model covariance U_k. Only applies to Gaussian models.
        %         * This method is REQUIRED if you want your model to be compatible
        %           with Gaussian TrackingX Filters (i.e. KalmanFilterX, etc).
        %         * Currently TrackingX filters do not support state-dependent
        %           process covariances. Extensions to state dependent 
        %           covariance should be possible.
        %           Contact Lyudmil @ sglvladi@gmail.com, if you require such a functionality.
        %
        % INPUTS:   k  - Time variable 
        %
        % OUTPUTS:  Uk - A (xDim x xDim) process covariance matrix.  
        %
        % USAGE:    Uk = CtrModel.ctr_cov(); % For time invariant process covariance
        %           Uk = CtrModel.ctr_cov(k); % For time variant process covariance
        %
        %   See also CTR, CTR_NOISE.
            Uk = 0;
        end
        
        function wk = ctr_noise(this, k, Ns)
        % CTR_NOISE Control input noise sample generator 
        %
        % INPUT:    k   - Time variable
        %           Ns  - The number samples to be generated 
        %                 (Optional, default is 1)  
        %
        % OUTPUTS:  wk - A (uDim x Ns) matrix of Ns process noise samples   
        %
        % USAGE:    wk = CtrModel.ctr_noise() should ONLY be allowed for time invariant noise
        %           wk = CtrModel.ctr_noise(Dt) 
        %           wk = CtrModel.ctr_noise(Dt, Ns) 
        %
        % See also CTR, CTR_COV.
            wk = zeros(this.Params.uDim,Ns);
        end
        
%         function ProbabilityMatrix = eval(this, k, xk, xkm1)
%         % EVAL Evaluates the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors  
%         % 
%         % INPUTS:   k    - Time variable
%         %           xk   - A (xDim x Np) matrix of Np new state vectors
%         %           xkm1 - A (xDim x Ns) matrix of Ns old state vectors
%         %   All inputs are REQUIRED. Set k=[] for time invariant systems
%         %
%         % OUTPUTS:  ProbabilityMatrix - A (Np x Ns) matrix of probabilities
%         %                               p(x_k|x_{k-1})    
%         %
%         % USAGE:    ProbabilityMatrix = eval(this, k, x_k, x_km1) 
%         %
%         % See also SYS, SYS_COV, SYS_NOISE.
%                               
%         end
        
    end
end