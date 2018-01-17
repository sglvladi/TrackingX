classdef ControlModelX
    % ControlModelX class base class/template
    %
    % Summary of ControlModelX
    %  This is the base TrackingX control model class, which should be used as
    %  the superclass of any custom defined dynamic models, compatible with
    %  the TrackingX toolkit. 
    %
    % ControlModelX Interface:
    %   Base properties:
    %        Params     - structure with the following compulsory fields:
    %        .xDim  - Number of state vector dimensions
    %        .uDim  - Number of control input vector dimensions
    %
    %   Base methods:
    %        ctr         - Control input vector process transition function b(u_k,w_k) 
    %        ctr_noise   - Control input noise sample generator function  
    %
    % The above parameters and methods are accessed by the majority of 
    % existing TrackingX library components and are COMPULSORY to guarantee 
    % compatibility with TrackingX.
    %   
    % You can also define any custom properties or methods that your model 
    % might require, or might assist you in defining your model, as long as
    % they do not collide with any of the above-defined properties and methods.
    %       
    % To keep things tidy, you are advised to place any custom properties 
    % within the "Params" structure property, e.g.:
    %       Params.b(Dt,~)      = Time-variant control input transition function 
    %                             handle b(k,wk)
    %       Params.U(Dt)        = Time-variant control input noise covariance 
    %                             function handle, returns (uDim x uDim) matrix
    %
    % For examples see also the source of LINEARNOISELESSCTRMODELX
    
    properties
        Params
    end
    
    methods
        function this = ControlModelX(Init)
        % CONTROLMODELX Constructor method
        %   
        % INPUTS:   Virtually anything you want! 
        %   
        % USAGE:    CtrModel = ControlModelX(Init); 
        %
        %   You are free to do ANYTHING you want in here.
        %   Since no current TrackingX components auto-initialise any models, the number and type of constructor arguments can be set to anything that works for you. 
        %
        %   See also LINEARNOISELESSCTRMODELX.
        
            % Validate .xDim 
%             if ~isfield(Params,'xDim')
%                 error('[CtrModel] Model dimensionality has not been specified!');
%             end
            
%             % Define an example transition function handle
%             Params.f = @() [5 1; 0 4]; % Linear invariant
%             %this.Params.f = @(k,xkm1,wk) [xkm1(2)*sin(k*xkm1(1))*wk(1); atan(xkm(2),xkm1(1)) + wk(2)]; % Non-Linear, state-time variant
%             
%             % Define an example process covariance function handle
%             Params.Q = @() [1 0; 0 1]; % Time invariant
%             %this.Params.Q = @(k) [k^2/2 k; k 1]; %Time variant
            
            % Initiate thisect properties
            this.Params = Init;      
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
        % USAGE:    uk_x = CtrModel.sys(k,uk,wk);
        %
        %           With a properly defined CtrModel.sys_noise(k) function,
        %           you can also use it to generate wk internally and propagate
        %           a number of samples as so:
        %           
        %           uk_x = CtrModel.sys(k,uk,CtrModel.sys_noise(k));        
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
%             
%             % Return control model covariance
%             Uk = this.Params.U();    % Time invariant
%             Uk = this.Params.U(k); % Time variant
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