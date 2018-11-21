classdef ConstantVelocityX < TransitionModelX
% ConstantVelocityX class
%
% Summary of ConstantVelocityX
% This is a class implementation of a time-varying 2D Linear-Gaussian 
% Constant Velocity Dynamic Model.
%
% The model is described by the following SDEs:
%   dx  = v*dt                    | Position (m)
%   dv  = q*dW(t),  W(t)~N(0,q^2)  | Speed    (m/s)
%
% Or equivalently:
%   x(t) = F(dt)*x(t-1) + w(t),  w(t)~ N(0,Q)
%
% where:
%   x = [x;
%        v];
%   F = [1  dt; 
%        0  1];
%   Q = [dt^3/3, dt^2/2, 
%        dt^2/2, dt];
%
% ConstantVelocityX Properties:
%   + VelocityErrVariance  Value of the noise diffusion coefficient q
%   + TimestepDuration         Value of the time variant dt
%
% ConstantVelocityX Methods:
%   + feval - Equivalent to applying the model transition equations 
%   + random - Process noise sample generator function
%   + pdf - Function to evaluate the probability p(x_k|x_{k-1}) of 
%           a set of new states, given a set of (particle) state vectors
%   + covariance - Returns the state covariance process Q_k 
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool

    properties (Access = private)
        F_cv = @(dt) [1 dt; 
                      0 1];
        Q_cv = @(dt,q) [dt^3/3 dt^2/2; 
                        dt^2/2, dt]*q;
        F
    end
    
    properties
        VelocityErrVariance = 1
        NumDims = 1
        TimestepDuration = duration(0,0,1);
    end
    
    properties (Dependent)
      NumStateDims
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if (isfield(config,'NumDims'))
                this.NumDims = config.NumDims;
            end
            if (isfield(config,'TimestepDuration'))
                this.TimestepDuration = config.TimestepDuration;
            end
            if (isfield(config,'VelocityErrVariance'))
                this.VelocityErrVariance = config.VelocityErrVariance;
                if(this.NumDims>length(this.VelocityErrVariance))
                    this.VelocityErrVariance = this.VelocityErrVariance(ones(1,this.NumDims),1);
                end
            end
                     
            this.F = @(dt) kron(eye(this.NumDims),this.F_cv(dt));
        end
        
        function Qk = Q(this,dt)
            Qk = [];
            for i = 1:this.NumDims
               Qk = blkdiag(Qk,this.Q_cv(dt,this.VelocityErrVariance(i)));
            end
        end
    end
    
    methods
        function this = ConstantVelocityX(varargin)
        % ConstantVelocityX Constructor method
        %   
        % Parameters
        % ----------
        % VelocityErrVariance: scalar or vector
        %   Describes the process noise diffusion coefficient. If NumDims
        %   is higher that 1 and VelocityErrVariance is scalar, then the
        %   same noise coefficient will be applied to all dimensions.
        %   Otherwise, VelocityErrVariance can be a vector of length
        %   NumDims, where the i-th element specifies the coefficient for
        %   the i-th dimension.
        % NumDims: scalar
        %   The number of dimensions. (default = 1)
        % TimestepDuration: duration
        %   The timestep duration. Specifying this on construction allows
        %   to avoid supplying it in every predict/update cycle.
        %
        % Usage
        % ----- 
        % * ConstantVelocityX(___,Name,Value) instantiates an object 
        %   handle, configured the parameters specified by one or
        %   more Name,Value pair arguments. 
        % * ConstantVelocityX(config) instantiates an object handle configured  
        %   with the parameters specified inside the 'config' structure, whose 
        %   fieldnames correspond to the parameter names as given in the 
        %   Parameters section above.
        %
        %  See also feval, random, pdf, covariance.   
            
            % Call SuperClass method
            this@TransitionModelX;
                        
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
            config = struct_concat(parser.Results, parser.Unmatched);
            this.initialise_(config);
            
        end
        
        function mat = matrix(this,dt)
        % matrix Return the model transition matrix
        %
        % Parameters
        % ----------
        % dt: duration, optional
        %   - A time duration variable (default: 1 sec)
        % 
        % Returns
        % -------
        % (NumStateDims x NumStateDims) matrix
        %   - The model transition matrix
            if(nargin==1)
                dt = this.TimestepDuration;
            end
            mat = this.F(seconds(dt));
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
        % (NumStateDims x Ns) or (NumStateDims x NumStateDims)  matrix
        %   - If no parameters are passed to the function, then xk will be the
        %     (NumStateDims x NumStateDims) state transition matrix.
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
                case 1 
                    xkm1 = 1;
                    wk   = 0;
                    dt = this.TimestepDuration;
                case 2
                    wk   = 0;
                    dt = this.TimestepDuration;
                case 3
                    dt = this.TimestepDuration;
                    if(islogical(wk) && wk)
                        wk = this.random(size(xkm1,2),dt);
                    end
            end
            
            % Compute result
            xk = this.matrix(dt)*xkm1 + wk;
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
        % (NumStateDims x NumStateDims) matrix
        %   The process noise covariance matrix
        %
        % Usage
        % -----
        % * Qk = covariance(this,dt) returns process covariance matrix, upon 
        %   application of dt.
        %
        % See also FEVAL, RANDOM, PDF.
            switch(nargin)
                case 1 
                    dt = this.TimestepDuration;
            end
            
            % Return process covariance
            cov = this.Q(seconds(dt)); % Time variant
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
              
            wk = gauss_rnd(zeros(this.NumStateDims,1),this.covar(dt),Ns);
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
        % See also FEVAL, RANDOM, COVAR
            
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
            numStateDims = 2*this.NumDims;
        end
        
    end
end