classdef DynamicModelX
    % DynamicModelX class base class/template
    %
    % Summary of DynamicModelX
    %  This is a base ProjectX dynamic model class, which should be used as
    %  the superclass of any custom defined dynamic models, compatible with
    %  the ProjectX toolkit. 
    %
    % DynamicModelX Interface:
    %   Base properties:
    %        Params     - structure with the following compulsory fields:
    %        .xDim  - Number of model/state dimensions
    %
    %   Base methods:
    %        sys         - State vector process transition function f(x_{k-1},w_k) 
    %        sys_noise   - Process noise sample generator function
    %        eval        - Function to evaluate the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors  
    %
    % The above parameters and methods are accessed by the majority of 
    % existing ProjectX library components and are COMPULSORY to guarantee 
    % compatibility with ProjectX.
    %   
    % You can also define any custom properties or methods that your model 
    % might require, or might assist you in defining your model, as long as
    % they do not collide with any of the above-defined properties and methods.
    %       
    % To keep things tidy, you are advised to place any custom properties 
    % within the "Params" structure property, e.g.:
    %       Params.q            = Process noise diffusion coefficient
    %       Params.f(Dt,~)      = Time-variant process transition function 
    %                             handle f(x_k|x_{k-1}), returns matrix.
    %       Params.Q(Dt)        = Time-variant process noise covariance 
    %                             function handle, returns (nx x nx) matrix
    %
    % For examples see also the source of CONSTANTVELOCITYMODELX and CONSTANDHEADINGMODELX
    
    properties
        Params
    end
    
    methods
        function this = DynamicModelX(Init)
        % DYNAMICMODELX Constructor method
        %   
        % INPUTS:   Virtually anything you want! 
        %   
        % USAGE:    DynModel = DynamicModelX(Init); 
        %
        %   You are free to do ANYTHING you want in here.
        %   Since no current ProjecX components auto-initialise any models, the number and type of constructor arguments can be set to anything that works for you. 
        %
        %   See also CONSTANTVELOCITYMODELX, CONSTANTHEADINGMODELX.
        
            % Validate .xDim 
%             if ~isfield(Params,'xDim')
%                 error('[DynModel] Model dimensionality has not been specified!');
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
        
        function xk = sys(this, k, xkm1, wk)
        % SYS State vector process transition function f(x_{k-1},w_k) 
        %
        % INPUTS:   k    - Time variable
        %           xkm1 - A (xDim x Ns) matrix of Ns state vectors from 
        %                  time k-1
        %           wk   - A (xDim x Ns) matrix Ns of process noise vectors, 
        %                  corresponding to the state vectors xkm1.  
        %
        %           The above arguments should be expected STRICTLY in the 
        %           above given order, to ensure compatibility with all ProjectX
        %           components.
        %
        % OUTPUTS:  xk  - A (xDim x Ns) matrix of Ns state vectors, which 
        %                 have been propagated through the dynamic model
        %        
        % USAGE:    xk = DynModel.sys(k,xkm1,wk);
        %
        %           With a properly defined DynModel.sys_noise(k,Ns) function,
        %           you can also use it to generate wk internally and propagate
        %           a number of samples as so:
        %           
        %           xk = DynModel.sys(k,xkm1,DynModel.sys_noise(k,size(xkm1,2)));
        %
        %   1. IMPORTANT NOTES:
        %       * In general, k is a REQUIRED argument. ALL ProjectX filters
        %         pass their own time variable as a function argument and 
        %         thus you should generally expect to handle/ignore it.
        %         k will ONLY work as an OPTIONAL input for the small subclass 
        %         of linear invariant models (i.e. both linear and state,time invariant),
        %         where both the time and state variables can be skipped. 
        %         (See Section 2 for more details)
        %       * xkm1 is an OPTIONAL input for linear models. For such models
        %         sys(k,x_km1) is equivalent to sys(k)¬xkm1. For example, 
        %         KalmanFilterX calls sys(k) to evaluate and extract the 
        %         transition matrix F(k) from the model.
        %       * xkm1 is ALWAYS passed as an argument by ALL non-linear 
        %         ProjectX filters (i.e. EKalmanFilterX, UKalmanFilterX, ParticleFilterX etc.)  
        %       * wk is an OPTIONAL input for all models and filters.
        %       * In any case, the first 3 argument spaces should ALWAYS be 
        %         reserved for k,xkm1,wk in the given order.
        %       * If you desire to add more input/output arguments, you can
        %         start from varargin{4} for inputs and varargout{2} for outputs 
        %
        %   2. Dynamic model notation conventions used herein:
        %       * (Non)Linear (state-time) invariant    : xk = sys(~,~,w(~))     = F(~) ¬ ~ ¬ wk(~) = f(~,~,wk(~)) (#)
        %       * Linear, state variant                 : xk = sys(~,xkm1,w(~))  = F(~) ¬ xkm1 ¬ wk(~) 
        %       * Linear, time variant                  : xk = sys(k,~,wk(~))    = F(k) ¬ ~ ¬ wk(k)
        %       * Linear, state-time variant            : xk = sys(k,xkm1,wk(k)) = F(k) ¬ xkm1 ¬ wk(k)
        %       * Non-Linear, time variant              : xk = sys(k,~,wk(k))    = f(k,~,wk(k))
        %       * Non-Linear, state variant             : xk = sys(~,xkm1,wk(~)) = f(~,xkm1,wk(~)) 
        %       * Non-Linear, state-time variant        : xk = sys(k,xkm1,wk(k)) = f(k,xkm1,wk(k))
        %
        %      (#)Non-Linear invariant models reduce to Linear invariant models. 
        %       ¬ denotes any linear operation
        %       ~ denotes any unnecessary arguments
        %
        %   3. Proposed sys() pseudo code, guaranteed to work with all ProjectX Filters:
        %       Define:
        %           anyValue(x): Returns true if argument x exists with any kind of value
        %           validValue(x): Returns true if argument x exist and has a valid value
        %           noArgument(x): Returns true if argument x does not exists
        %           noArguments : Returns true if no arguments are passed
        %       Switch(ModelType):
        %       * (Non)Linear (state-time) invariant:
        %           If (noArguments):
        %               return F;
        %           Elif (anyValue(k,xkm1) AND validValue(wk))
        %               return F ¬ wk;
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, state variant :
        %           If (validValue(xkm1,wk) AND anyValue(k))
        %               return F ¬ xkm1 ¬ wk;
        %           Elif (validValue(xkm1) AND anyValue(k) AND noArgument(wk))
        %               return F ¬ xkm1;
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, time variant :
        %           If (validValue(k,wk) AND anyValue(xkm1))
        %               return F(k) ¬ wk(k);
        %           Elif (validValue(k) AND (anyValue(xkm1) OR noArgument(xkm1)) AND noArgument(wk))
        %               return F(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, state-time variant :
        %           If (validValue(k,xkm1,wk))
        %               return F(k) ¬ xkm1 ¬ wk(k);
        %           Elif (validValue(k,xkm1) AND noArgument(wk))
        %               return F(k) ¬ xkm1;
        %           Elif (validValue(k)) AND noArgument(xkm1, wk))
        %               return F(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, time variant              
        %           If (validValue(k,wk) AND anyValue(xkm1))
        %               return f(k,~,wk(k));
        %           Elif (validValue(k) AND (anyValue(xkm1) OR noArgument(xkm1)) AND noArgument(wk))
        %               return f(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, state variant 
        %           If (validValue(xkm1,wk) AND anyValue(k))
        %               return f(~,xkm1,wk(~));
        %           Elif (validValue(xkm1) AND (anyValue(k)) AND noArgument(wk))
        %               return f(~,xkm1,~);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, state-time variant
        %           If (validValue(k,xkm1,wk))
        %               return f(k,xkm1,wk(k));
        %           Elif (validValue(k,xkm1) AND noArgument(wk))
        %               return f(k,xkm1,~);
        %           Else
        %               ERROR;
        %           EndIf;
        %       EndSwitch; 
        %
        %   [README!] It is worth noting that, if you ensure that your code abides to the above-define sudo code rules,
        %             then you could skip input validation completely, to remove the additional computation overhead.
        %
        %   See also SYS_COV, SYS_NOISE.
        
        end
        
        function Qk = sys_cov(this, k)
        % SYS_COV Returns process covariance Q_k. Only applies to Gaussian models.
        %         * This method is REQUIRED if you want your model to be compatible
        %           with Gaussian ProjectX Filters (i.e. KalmanFilterX, etc).
        %         * Currently ProjectX filters do not support state-dependent
        %           process covariances. Extensions to state dependent 
        %           covariance should be possible.
        %           Contact Lyudmil @ sglvladi@gmail.com, if you require such a functionality.
        %
        % INPUTS:   k  - Time variable 
        %
        % OUTPUTS:  Qk - A (xDim x xDim) process covariance matrix.  
        %
        % USAGE:    Qk = DynModel.sys_cov(); % For time invariant process covariance
        %           Qk = DynModel.sys_cov(k); % For time variant process covariance
        %
        %   See also SYS, SYS_NOISE.
%             
%             % Return process covariance
%             Qk = this.Params.Q();    % Time invariant
%             Qk = this.Params.Q(p.Results.k); % Time variant
        end
        
        function wk = sys_noise(this, k, Ns)
        % SYS_NOISE Process noise sample generator 
        %
        % INPUT:    k   - Time variable
        %           Ns  - The number samples to be generated 
        %                 (Optional, default is 1)  
        %
        % OUTPUTS:  wk - A (xDim x Ns) matrix of Ns process noise samples   
        %
        % USAGE:    wk = DynModel.sys_noise() should ONLY be allowed for time invariant noise
        %           wk = DynModel.sys_noise(Dt) 
        %           wk = DynModel.sys_noise(Dt, Ns) 
        %
        % See also SYS, SYS_COV.
        %
        %     % Generate noise samples
        %     wk = mvnrnd(zeros(onj.Params.,Ns)',this.config.Q(Dt))';
        end
        
        function ProbabilityMatrix = eval(this, k, xk, xkm1)
        % EVAL Evaluates the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors  
        % 
        % INPUTS:   k    - Time variable
        %           xk   - A (xDim x Np) matrix of Np new state vectors
        %           xkm1 - A (xDim x Ns) matrix of Ns old state vectors
        %   All inputs are REQUIRED. Set k=[] for time invariant systems
        %
        % OUTPUTS:  ProbabilityMatrix - A (Np x Ns) matrix of probabilities
        %                               p(x_k|x_{k-1})    
        %
        % USAGE:    ProbabilityMatrix = eval(this, k, x_k, x_km1) 
        %
        % See also SYS, SYS_COV, SYS_NOISE.
                              
        end
        
    end
end