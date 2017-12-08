classdef ObservationModelX
    % ObservationModelX class template
    %
    % Summary of ObservationModelX
    % This is a class template, which can be used as a guidline when defining custom observation model subclasses, compatible with the ProjectX toolkit.
    % Any custom defined Observation Model should be derived from this ObservationModelX base class.
    %
    % ObservationModelX Interface:
    %   Properties:
    %       - Params = structure with the following compulsory fields:
    %           .xDim = Number of model/state dimensions 
    %           .yDim = Number of observation dimensions 
    %
    %   Methods:
    %       obs                - State-space to measurement-space state projection function h(k,x_k,v_k)
    %       obs_cov            - Observation noise covariance R_k (Only REQUIRED for Gaussian Filters)
    %       obs_noise          - Observation noise sample generator
    %       sample             - Sample from measurement model, given a mean and number of samples
    %       eval_likelihood    - Evaluates the likelihood p(y_k|x_k) of a set of measurements, given a set of (particle) states  
    %
    % The above parameters and methods are accessed by the majority of existing ProjectX library components and are COMPULSORY to guarantee compatibility with ProjectX.
    %   
    % You can also define any custom properties or methods that your model might require, or might assist you in defining your model,
    %  as long as they do not colide with any of the above-defined properties and methods.
    % Such properties and methods can be used for your own purposes and directly accessed from outside the class, (unless you explicitly choose not to).
    %       
    % To keep things tidy, you are advised to place any custom properties within the Params structure, e.g.:
    %       Params.r            = Process noise diffusion coefficient
    %       Params.h(Dt,~)      = Time-variant observation model function handle h(y_k|x_k), returns (.yDim x N) matrix.
    %       Params.R(Dt)        = Time-variant observation noise covariance function handle, returns (.yDim x .yDim) matrix
    %
    % For examples see also the source of PositionalObsModel
    
    properties
        Params
    end
    
    methods
        function obj = ObservationModelX(Params)
        % ObservationModelX - Constructor method
        %   
        %   Inputs:
        %       Virtually anything you want! For this example we assume a "Params" structure, but you are free to do anything. 
        %   
        %   Usage:
        %       ObsModel = ObservationModelX(Params); 
        %
        %   You are free to do ANYTHING you want in here, as long as the obj.Params.xDim and obj.Params.yDim properties are correctly initialised by the end of the method.
        %   Since no current ProjecX components auto-initialise any models, the number and type of constructor arguments can be set to anything that works for you. 
        %
        %   See also PositionalObsModelX.
        
            % Validate .xDim 
            if ~isfield(Params,'xDim')
                error('[ObsModel] State dimensionality has not been specified!');
            end
            if ~isfield(Params,'yDim')
                error('[ObsModel] Observation dimensionality has not been specified!');
            end
            
%             % Define an example transition function handle
%             Params.h = @() [5 1; 0 4]; % Linear invariant
%             %obj.Params.h = @(k,xk,vk) [xk(2)*sin(k*xk(1))*vk(1); atan(xk(2),xk(1)) + vk(2)]; % Non-Linear, state-time variant
%             
%             % Define an example process covariance function handle
%             Params.R = @() [1 0; 0 1]; % Time invariant
%             %obj.Params.R = @(k) [k^2/2 k; k 1]; %Time variant
            
            % Initiate object properties
            obj.Params = Params;      
        end
        
        function xk = obs(obj, varargin)
        % obs - State-space to measurement-space state projection function h(k,xk,vk)
        %
        %   Inputs:
        %       k      : Time variable
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
        %   1. IMPORTANT NOTES:
        %       * In general, k is a REQUIRED argument. ALL ProjectX filters pass their own time variable as a function argument and thus you should generally expect to handle/ignore it.
        %           k will ONLY work as an OPTIONAL input for the small subclass of linear invariant models (i.e. both linear and state,time invariant), where both the time and state variables can be skipped. (See Section 2 for more details)
        %       * xk is an OPTIONAL input for linear models. For such models obs(k,xk) is equivalent to obs(k)¬xk(¬ denotes any linear operation). For example, KalmanFilterX calls obs(k) to evaluate and extract the transition matrix H(k) from the model.
        %           xk is ALWAYS passed as an argument by ALL non-linear ProjectX filters (i.e. EKalmanFilterX, UKalmanFilterX, ParticleFilterX etc.)  
        %       * vk is an OPTIONAL input for all models and filters.
        %       * In any case, the first 3 argument spaces should ALWAYS be reserved for k,xk,vk in the given order.
        %       * If you desire to add more input/output arguments, you can start from varargin{4} for inputs and varargout{2} for outputs
        %       * An input parser example is available below, which allows inputs to be passed in name,value format. 
        %
        %   2. Observation model notation conventions used herein:
        %       * (Non)Linear (state-time) invariant    : yk = obs(~,~,v(~))     = H(~)¬ ~ ¬ vk(~) = h(~,~,vk(~)) (#)
        %       * Linear, state variant                 : yk = obs(~,xk,v(~))    = H(~) ¬ xk ¬ vk(~) 
        %       * Linear, time variant                  : yk = obs(k,~,vk(~))    = H(k) ¬ ~ ¬ vk(k)
        %       * Linear, state-time variant            : yk = obs(k,xk,vk(k))   = H(k) ¬ xk ¬ vk(k)
        %       * Non-Linear, time variant              : yk = obs(k,~,vk(k))    = h(k,~,vk(k))
        %       * Non-Linear, state variant             : yk = obs(~,xk,vk(~))   = h(~,xk,vk(~)) 
        %       * Non-Linear, state-time variant        : yk = obs(k,xk,vk(k))   = h(k,xk,vk(k))
        %
        %      (#)Non-Linear invariant models reduce to Linear in variant models. 
        %       ¬ denotes any linear operation
        %       ~ denotes any unnecessary arguments
        %
        %   3. Proposed obs() pseudo code, guaranteed to work with all ProjectX Filters:
        %       Define:
        %           anyValue(x): Returns true if argument x exists with any kind of value
        %           validValue(x): Returns true if argument x exist and has a valid value
        %           noArgument(x): Returns true if argument x does not exists
        %           noArguments : Returns true if no arguments are passed
        %       Switch(ModelType):
        %       * (Non)Linear (state-time) invariant:
        %           If (noArguments):
        %               return F;
        %           Elif (anyValue(k,xk) AND validValue(vk))
        %               return H ¬ vk;
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, state variant :
        %           If (validValue(xk,vk) AND anyValue(k))
        %               return H ¬ xk ¬ vk;
        %           Elif (validValue(xkm1) AND anyValue(k) AND noArgument(vk))
        %               return H ¬ xk;
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, time variant :
        %           If (validValue(k,vk) AND anyValue(xk))
        %               return H(k) ¬ vk(k);
        %           Elif (validValue(k) AND (anyValue(xk) OR noArgument(xk)) AND noArgument(vk))
        %               return H(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Linear, state-time variant :
        %           If (validValue(k,xk,vk))
        %               return H(k) ¬ xk ¬ vk(k);
        %           Elif (validValue(k,xk) AND noArgument(vk))
        %               return H(k) ¬ xk;
        %           Elif (validValue(k)) AND noArgument(xk, vk))
        %               return H(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, time variant              
        %           If (validValue(k,vk) AND anyValue(xk))
        %               return h(k,~,vk(k));
        %           Elif (validValue(k) AND (anyValue(xk) OR noArgument(xk)) AND noArgument(vk))
        %               return h(k);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, state variant 
        %           If (validValue(xk,vk) AND anyValue(k))
        %               return h(~,xk,vk(~));
        %           Elif (validValue(xk) AND (anyValue(k)) AND noArgument(vk))
        %               return h(~,xk,~);
        %           Else
        %               ERROR;
        %           EndIf;
        %       * Non-Linear, state-time variant
        %           If (validValue(k,xk,vk))
        %               return h(k,xk,vk(k));
        %           Elif (validValue(k,xk) AND noArgument(vk))
        %               return h(k,xk,~);
        %           Else
        %               ERROR;
        %           EndIf;
        %       EndSwitch; 
        %
        %   [README!] It is worth noting that, if you ensure that your code abides to the above-define sudo code rules,
        %             then you could skip input validation completely, to remove the additional computation overhead.
        %
        %   See also obs_cov, obs_noise.
        
            % Example for linear (state-time) invariant models
            % =================================================>
%             % Define input parser
%             p = inputParser;
%             
%             % Define dafault values to suit your needs
%             default_k = [];
%             default_xk = [];
%             default_vk = [];
%             
%             % Add k as an optional parameter
%             validate_k = @(x) (isempty(x) || (isscalar(x) && isnumeric(x))); % Only valid either if k is empty, or if k is 
%             addOptional(p, 'k', default_k, validate_k);
%             % Add xk as an optional parameter
%             validate_xk = @(x) (isempty(x) || (isnumeric(x)&&(size(x,1) == obj.Params.xDim)));
%             addOptional(p, 'xk', default_xk, validate_xk);
%             % Add vk as an optional parameter
%             validate_vk = @(x) (isempty(x) || (isnumeric(x)&&(size(x,1) == obj.Params.xDim)));
%             addOptional(p, 'vk', default_vk, validate_vk);
%             parse(p,varargin{:});
%             
%             if(isempty(p.Results.vk))
%                 yk = obj.Params.h();
%             else
%                 yk = obj.Params.h() + p.Results.vk;
%             end
        
        end
        
        function Rk = obs_cov(obj, varargin)
        % obs_cov - Returns measurement covariance R_k. Only applies to Gaussian models.
        %           REQUIRED for your model to work with KalmanFilterX, EKalmanFilterX, UKalmanFilterX and all other Gaussian filters.
        %           Currently ProjectX filters do not support state-dependent measurement covariances. Extensions to state dependent covariance should be possible.
        %           Contact Lyudmil @ sglvladi@gmail.com, if you require such a functionality.
        %
        %   Inputs:
        %       k: Time variable 
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
            
%             % Define input parser
%             p = inputParser;
%             
%             % Define dafault values to suit your needs
%             default_k = [];
%             validate_k = @(x) (isempty(x) || (isscalar(x) && isnumeric(x))); % Only valid if k is empty, or if k is scallar and numeric
%             addOptional(p, 'k', default_k, validate_k);
%             parse(p,varargin{:});
%             
%             % Return process covariance
%             Rk = obj.Params.R();    % Time invariant
%             % Rk = obj.Params.R(p.Results.k); % Time variant
        end
        
        function vk = obs_noise(obj, varargin)
        % obs_noise - Process noise sample generator 
        %
        %   Inputs:
        %       k  : Time variable
        %       Ns : The number samples to be generated (Optional, default is 1)  
        %
        %   Outputs:
        %       vk: a (.yDim x Ns) matrix of Ns process noise samples, where nx is the dimensionality of the state   
        %
        %   Usage:
        %   vk = ObsModel.obs_noise() should ONLY be allowed for time invariant noise
        %   vk = ObsModel.obs_noise(k) 
        %   vk = ObsModel.obs_noise(k, Ns) 
        %
        %   See also obs, obs_cov.
%             % Define input parser
%             p = inputParser;
%             
%             % Define default Ns
%             default_Ns = 1;
%             
%             % Add k as a REQUIRED argument
%             validate_k = @(x) (isscalar(x) && isnumeric(x)); % Only valid if k is scallar and numeric
%             addRequired(p, 'k', validate_k);
%             validate_Ns = @(x) (isscalar(x) && isnumeric(x) && x>0); % Only valid if k is scallar and numeric
%             addOptional(p, 'Ns', default_Ns, validate_Ns);
%             parse(p, varargin{:});
%         
%             % Generate noise samples
%             vk = mvnrnd(zeros(obj.Params.yDim,Ns)',obj.config.R(p.Results.k))';
        end
        
        function LikelihoodMatrix = eval(obj, varargin)
        % eval - Evaluates the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors  
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
        %   LikelihoodMatrix = eval_likelihood(obj, Dt, x_k, x_km1) Evaluates and produces a (Np x Ns) likelihood matrix.
        %
        %   See also obs, obs_cov, obs_noise, sample.
        
%             % Define input parser
%             p = inputParser;
%             
%             % Add k as a required parameter
%             validate_k = @(x) (isscalar(x) && isnumeric(x)); % Only valid either if k is empty (time invariant models), or if k is 
%             addRequired(p, 'k', validate_k);
%             % Add yk as a required parameter
%             validate_yk = @(x) (isnumeric(x)&&(size(x,1) == obj.Params.yDim));
%             addRequired(p, 'yk', validate_yk);
%             % Add xk as a required parameter
%             validate_xk = @(x) (isnumeric(x)&&(size(x,1) == obj.Params.xDim));
%             addOptional(p, 'xk', validate_xk);
%             parse(p, varargin{:});
%         
%             LikelihoodMatrix = zeros(size(p.Results.yk,2), size(p.Results.xk,2));
%             if(issymmetric(obj.config.R(p.Results.k)) && size(p.Results.xk,2)>size(p.Results.yk,2))
%                 % If R is symmetric and the number of state vectors is higher than the number of measurements,
%                 %   then evaluate N(x_k; y_k, R) = N(y_k; x_k, R) to increase speed
%                 for i=1:size(yk,2)
%                     LikelihoodMatrix(i,:) = mvnpdf(p.Results.xk', p.Results.yk(:,i)', obj.config.R(p.Results.k))';
%                 end
%             else
%                 for i=1:size(xk,2)
%                     LikelihoodMatrix(:,i) = mvnpdf(p.Results.yk', p.Results.xk(:,i)', obj.config.R(p.Results.k))';  
%                 end
%              end
%                         
        end
        
    end
end