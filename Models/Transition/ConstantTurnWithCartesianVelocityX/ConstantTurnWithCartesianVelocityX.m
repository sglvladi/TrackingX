classdef ConstantTurnWithCartesianVelocityX <  TransitionModelX 
% ConstantTurnWithCartesianVelocityX class
%
% Summary of ConstantTurnWithCartesianVelocityX
% This is a class implementation of a time-varying 2D Linear-Gaussian 
% Constant TurnRate Transition Model. [1]
%
% The model is described by the following SDEs:
%   TODO: Fill in the SDEs
%
% Or equivalently:
%   x_t = f_ct(x_{t-1},dt) + G_ct*w_t,  w_t ~ N(0,Q_ct)
%
% where: 
%   x = [x; x_vel; y; y_vel; omega]
%
%   f_ct(x,dt) = [ x + x_vel*sin(omega*dt)/omega - y_vel*(1-cos(omega*dt))/omega; 
%                   x_vel*cos(omega*dt) - y_vel*sin(omega*dt); 
%                   y + x_vel*(1-cos(omega*dt))/omega + y_vel*sin(omega*dt)/omega; 
%                   x_vel*sin(omega*dt) + y_vel*cos(omega*dt);
%                   omega ];
%
%   w = [ dt^2/2        0               0;     
%         dt            0               0;
%         0             dt^2/2          0;  *  [w_x, w_y, w_w]';
%         0             dt              0;
%         0             0               dt] 
%
% [1] M. Roth, G. Hendeby and F. Gustafsson, "EKF/UKF maneuvering target tracking 
%     using coordinated turn models with polar/Cartesian velocity," 17th International 
%     Conference on Information Fusion (FUSION), Salamanca, 2014, pp. 1-8.
%
% ConstantVelocityModelX_2D Properties:
%   - VelocityErrVariance     The process noise (co)variance
%   - TurnRateErrVariance     The turn rate error variance
%   - TimestepDuration        Value of the time variant dt
%
% ConstantVelocityModelX_2D Methods:
%   - feval(~)         Equivalent to applying the model transition equations 
%   - random(~)        Process noise sample generator function
%   - pdf(~)           Function to evaluate the probability p(x_k|x_{k-1}) of 
%                       a set of new states, given a set of (particle) state vectors
%                       e.g. eval = @(xk,xkm1) mvnpdf(xk,xkm1,Q);
%   - covar(~)         Returns the state covariance process Q_k 
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool
    
    properties
        VelocityErrVariance = [1 1]
        TurnRateErrVariance = pi/180
        TimestepDuration = duration(0,0,1);
    end
    
    properties (Dependent)
      NumStateDims
    end
    
    properties (Access = private)
        f_ct = @(x,dt) [x(1,:)+x(2,:).*sin(x(5,:)*dt)./x(5,:)-x(4,:).*(1-cos(x(5,:)*dt))./x(5,:); 
                         x(2,:)*cos(dt*x(5,:)) - x(4,:)*sin(dt*x(5,:)); 
                         x(3,:)+x(2,:).*(1-cos(x(5,:)*dt))./x(5,:)+x(4,:).*sin(x(5,:)*dt)./x(5,:); 
                         x(2,:).*sin(dt*x(5,:)) + x(4,:).*cos(dt*x(5,:));
                         x(5,:)];
                      
        G_ct = @(dt) [dt^2/2,    0,      0;
                       dt,        0,      0;
                       0          dt^2/2, 0;
                       0          dt,     0;
                       0          0 ,     dt];
        
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if (isfield(config,'TimestepDuration'))
                this.TimestepDuration = config.TimestepDuration;
            end
            if (isfield(config,'VelocityErrVariance'))
                this.VelocityErrVariance = config.VelocityErrVariance;
                if(length(this.VelocityErrVariance)==1)
                    this.VelocityErrVariance = this.VelocityErrVariance(ones(1,2),1);
                end
            end
            if (isfield(config,'TurnRateErrVariance'))
                this.TurnRateErrVariance = config.TurnRateErrVariance;
            end
        end
        
        function Q_k =  Q_ct(this)
           Q_k =  blkdiag(diag(this.VelocityErrVariance),this.TurnRateErrVariance);
        end
    end
    
    methods
        function this = ConstantTurnWithCartesianVelocityX(varargin)
        % ConstantTurnWithCartesianVelocityX Constructor method
        %   
        % Parameters
        % ----------
        % VelocityErrVariance: scalar or (1 x 2) vector
        %   Describes the variance of the zero-mean (Cartesian) velocity 
        %   process noise in m/s^2. If VelocityErrVariance is scalar, then 
        %   the noise is modeled as iid on each Cartesian dimension. Otherwise, 
        %   VelocityErrVariance can be a (1 x 2) vector, whose elements 
        %   specify the variance for each Cartesian dimension.
        % TurnRateErrVariance: scalar
        %   The variance of the zero-mean turn-rate diffusion process,
        %   specified in rad/s^2
        % TimestepDuration: duration
        %   The timestep duration. Specifying this on construction allows
        %   to avoid supplying it in every predict/update cycle.
        %
        % Usage
        % ----- 
        % * ConstantTurnWithCartesianVelocityX(___,Name,Value) instantiates an object 
        %   handle, configured with parameters specified by one or more Name,Value 
        %   pair arguments. 
        % * ConstantTurnWithCartesianVelocityX(config) instantiates an object  
        %   handle configured with the parameters specified inside the 'config'
        %   structure, whose fieldnames correspond to the parameter names as
        %   given in the Parameters section above.
        %
        %  See also feval, random, pdf, covariance.   
            
            % Call SuperClass method
            this@TransitionModelX;
            
            % Return quickly if no arguments are passed
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config)
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config)
        end
        
        function xk = feval(this, xkm1, wk, dt)
        % FEVAL Propagate a given state through the dynamic model
        %
        % Parameters
        % ----------
        % xkm1: (NumStateDims x Ns) matrix, optional
        %   - A matrix, whose columns correspond to individual state vectors.
        %   - If not provided, then the state transition matrix (F) will be 
        %     returned.
        % wk: (NumStateDims x Ns) matrix or boolean, optional
        %   - If wk is a boolean and set True, then the default model noise 
        %     generation function (this.random()) will be used to generate the 
        %     noise. 
        %   - Otherwise, wk should be a matrix of equal dimensions to xk, 
        %     where each column corresponds to the random noise which will
        %     be added to each state vector. 
        %   - If not provided, then no noise will be added to the state vectors.
        % dt: duration, optional
        %   A time variant. (default: duration(0,0,1))
        %
        % Returns
        % -------
        % xk: (NumStateDims x Ns) matrix or function handle
        %   - If no parameters are passed to the function, then xk will be the
        %     state transition function f.
        %   - Else, xk will be a (NumStateDims x Ns) matrix, whose columns 
        %     correspond to the columns of xkm1, each propagated through the
        %     dynamic model
        %   
        % Usage
        % -----
        % * xk = FEVAL(this,xkm1,wk,dt) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time dt and adding random noise wk. xk, xkm1 and wk must have  
        %   the same dimensions, i.e (NumStateDim x Ns), where Ns is the 
        %   number of states/columns in xkm1.
        % * xk = FEVAL(this,xkm1,wk) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time this.TimestepDuration and adding random noise wk.
        %   (i.e. Default dt = this.TimestepDuration)
        % * xk = FEVAL(this,xkm1) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time this.TimestepDuration without the addition of noise.
        %   (i.e. Default wk = 0)
        % * xk = FEVAL(this) returns the model's transition matrix resulting
        %   through the application of this.TimestepDuration. 
        %   (i.e. Default xkm1 = 1)
        %
        % See also PDF, RANDOM, COVARIANCE.
        
            switch(nargin)
                case 2
                    wk   = 0;
                    dt = this.TimestepDuration;
                case 3
                    dt = this.TimestepDuration;
                    if(islogical(wk) && wk)
                        wk = this.random(size(xkm1,2),dt);
                    end
                case 4
                    if(islogical(wk) && wk)
                        wk = this.random(size(xkm1,2),dt);
                    end
            end
            
            % Compute result
            xk = this.f_ct(xkm1,seconds(dt)) + wk;
        end
        
        function cov = covar(this,dt)
        % COVARIANCE Returns process covariance matrix.
        %
        % Parameters
        % ----------
        % dt: duration, optional
        %   A time variant. (default: duration(0,0,1))
        %
        % Returns
        % -------
        % cov: (NumStateDims x NumStateDims) matrix
        %   The process noise covariance matrix
        %
        % Usage
        % -----
        % * Qk = covariance(this,dt) returns process covariance matrix, upon 
        %   application of dt.
        % * Qk = covariance(this) returns process covariance matrix, upon 
        %   application of this.TimestepDuration. 
        % (i.e. Default dt = this.TimestepDuration)
        %
        % See also FEVAL, RANDOM, PDF.
            switch(nargin)
                case 1 
                    dt = this.TimestepDuration;
            end
            
            % Return process covariance
            cov = this.G_ct(seconds(dt))*this.Q_ct()*this.G_ct(seconds(dt))'; % Time variant
        end
        
        function wk = random(this, Ns, dt)
        % RANDOM Generates random samples from the dynamic model's noise
        % distribution.
        %
        % Parameters
        % ----------
        % Ns: scalar, optional
        %   The number of samples to be generated.
        %   (default = 1)
        % dt: duration, optional
        %   A time variant. 
        %   (default=this.TimestepDuration)
        %
        % Returns
        % -------
        % wk: (NumStateDims x Ns)
        %   A matrix, whose columns correspond to indivual/independent
        %   noise samples.
        %
        % Usage
        % -----
        % * wk = random(this,Ns,dt) generates and returns a set of Ns samples
        %   generated from the noise distribution of the dynamic model, i.e.
        %   wk ~ N(0,Q(dt)), where Q is the noise covariance upon application 
        %   of the time variant dt.
        % * wk = random(this,Ns) generates and returns a set of Ns samples
        %   generated from the noise distribution of the dynamic model, i.e.
        %   wk ~ N(0,Q), where Q is the noise covariance upon application 
        %   of the time variant this.TimestepDuration.
        % (i.e. Default dt = this.TimestepDuration)
        %
        % See also FEVAL, PDF, COVARIANCE
       
            switch(nargin)
                case 1
                    Ns = 1;
                    dt = this.TimestepDuration;
                case 2
                    dt = this.TimestepDuration;
            end
              
            wk = this.G_ct(seconds(dt))*gauss_rnd(zeros(3,1),this.Q_ct(),Ns);
        end
        
        function prob = pdf(this, xk, xkm1, dt)
        % PDF Evaluates the probability/likelihood p(x_k|x_{k-1}) of a 
        % (set of) new state vector(s), given a (set of) old state vector(s)  
        % 
        % Parameters
        % ----------
        % xk: (NumStateDims x Ns) matrix
        %   A matrix, whose columns correspond to individual new state vectors.
        % xkm1: (NumStateDims x Np) matrix
        %   A matrix, whose columns correspond to individual old state vectors
        % dt: duration, optional
        %   A time variant. 
        %   (default=this.TimestepDuration)
        %
        % Returns
        % -------
        % prob: (Np x Ns) matrix
        %   A matrix, where each element (i,j) corresponds to the evaluated
        %   probability p(xk(:,j)|xkm1(:,i))
        %
        % Usage
        % -----
        % * prob = pdf(x_k, x_km1, dt) evaluates and returns a (a x b)
        %   probability/likelihood matrix given the (NumStatesDim x a) x_k
        %   and (NumStatesDim x b) x_km1 state matrices. dt is an optional 
        %   argument (Default = this.TimestepDuration) time variable, which is 
        %   used for computing the model's covariance.
        %
        % See also FEVAL, RANDOM, COVARIANCE
            
             if(nargin<4)
                dt = this.TimestepDuration;
            end
            
            xk_km1 = this.feval(xkm1,false,dt);
            prob = zeros(size(xk,2), size(xkm1,2));
            if(size(xkm1,2)>size(xk,2))
                for i=1:size(xk,2)
                    prob(i,:) = gauss_pdf(xk(:,i), xk_km1, this.covar(dt));
                end
            else
                for i=1:size(xkm1,2)
                    prob(:,i) = gauss_pdf(xk, xk_km1(:,i), this.covar(dt))';  
                end
            end
                        
        end
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        function numStateDims = get.NumStateDims(this)
            numStateDims = 5;
        end
    end
end