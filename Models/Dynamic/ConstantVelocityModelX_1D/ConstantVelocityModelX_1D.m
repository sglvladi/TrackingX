classdef ConstantVelocityModelX_1D < DynamicModelX
% ConstantVelocityModelX_1D class
%
% Summary of ConstantVelocityModelX_1D
% This is a class implementation of a time-varying 1D Linear-Gaussian 
% Constant Velocity Dynamic Model.
%
% The model is described by the following SDEs:
%   dX  = V*dt                     | Position on X-axis (m)
%   dV  = q*dW(t),  W(t)~N(0,q^2)  | Speed on X-axis    (m/s)
%
% Or equivalently:
%   x(t) = F(Dt)*x(t-1) + w(t),  w(t)~ N(0,Q)
%
% where:
%   x = [X; V]
%   F = [1 Dt; 0 1]
%   Q = [Dt^3/3, Dt^2/2; Dt^2/2, Dt]*q; 
%
% ConstantVelocityModelX_1D Properties:
%   + VelocityErrVariance  Value of the noise diffusion coefficient q
%   + TimeVariant         Value of the time variant Dt
%
% ConstantVelocityModelX_1D Methods:
%   + feval - Equivalent to applying the model transition equations 
%   + random - Process noise sample generator function
%   + pdf - Function to evaluate the probability p(x_k|x_{k-1}) of 
%           a set of new states, given a set of (particle) state vectors
%   + covariance - Returns the state covariance process Q_k 
%
    properties (Access = private)
        F = @(Dt) [1 Dt; 0 1];
        Q = @(Dt) [Dt^3/3 Dt^2/2; Dt^2/2, Dt];
    end
    properties
        VelocityErrVariance
        TimeVariant
    end
    
    methods
        function this = ConstantVelocityModelX_1D(varargin)
        % CONSTANTVELOCITYMODEL_1D Constructor method
        % 
        % Parameters
        % ----------
        % VelocityErrVariance: scalar
        %   A value describing the process noise diffusion coefficient.
        %
        % Usage
        % -----
        % * ConstantVelocityModelX_1D(q) instantiates an object handle 
        %   configured with the provided process noise diffusion
        %   coefficient q (scalar). 
        % * ConstantVelocityModelX_1D(config) instantiates an object handle 
        %   configured with the provided process noise diffusion
        %   coefficient config.VelocityErrVariance (scalar).
        % * ConstantVelocityModelX_1D(___,Name,Value) instantiates an object 
        %   handle, configured with additional options specified by one or
        %   more Name,Value pair arguments. 
        %
        %  See also feval, random, pdf, covariance.   
            
        
            % Call SuperClass method
            this@DynamicModelX();
            
            % Return quickly if no arguments are passed
            if(nargin==0)
                this.TimeVariant = 1;
                this.VelocityErrVariance = 1;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    this.VelocityErrVariance = varargin{1}.VelocityErrVariance;
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addOptional('VelocityErrVariance',[]);
            parser.parse(varargin{:});
            
            this.VelocityErrVariance = parser.Results.VelocityErrVariance;
            
            this.NumStateDims = 2;
            
            this.TimeVariant = 1;
            
            this.Q = @(Dt) [Dt^3/3 Dt^2/2; Dt^2/2, Dt]*this.VelocityErrVariance;
            
        end
        
        function xk = feval(this, xkm1, wk, Dt)
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
        % Dt: scalar, optional
        %   A time variant. (default=this.TimeVariant)
        %
        % Returns
        % -------
        % xk: (NumStateDims x Ns) or (NumStateDims x NumStateDims)  matrix
        %   - If no parameters are passed to the function, then xk will be the
        %     (NumStateDims x NumStateDims) state transition matrix.
        %   - Else, xk will be a (NumStateDims x Ns) matrix, whose columns 
        %     correspond to the columns of xkm1, each propagated through the
        %     dynamic model
        %   
        % Usage
        % -----
        % * xk = FEVAL(this,xkm1,wk,Dt) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time Dt and adding random noise wk. xk, xkm1 and wk must have  
        %   the same dimensions, i.e (NumStateDim x Ns), where Ns is the 
        %   number of states/columns in xkm1.
        % * xk = FEVAL(this,xkm1,wk) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time this.TimeVariant and adding random noise wk.
        % * xk = FEVAL(this,xkm1) returns the new state xk, produced by
        %   propagating the given state xkm1 through the dynamic model for
        %   time this.TimeVariant without the addition of noise.
        %   (i.e. Default wk = 0)
        % * xk = FEVAL(this) returns the model's transition matrix resulting
        %   through the application of this.TimeVariant. 
        %   (i.e. Default xkm1 = 1)
        %
        % See also PDF, RANDOM, COVARIANCE.
        
            switch(nargin)
                case 1 
                    xkm1 = 1;
                    wk   = 0;
                    Dt = this.TimeVariant;
                case 2
                    wk   = 0;
                    Dt = this.TimeVariant;
                case 3
                    Dt = this.TimeVariant;
                    if(islogical(wk) && wk)
                        wk = this.random(size(xkm1,2),Dt);
                    end
            end
            
            % Compute result
            xk = this.F(Dt)*xkm1 + wk;
        end
        
        function cov = covariance(this,Dt)
        % COVARIANCE Returns process covariance matrix.
        %
        % Parameters
        % ----------
        % Dt: scalar, optional
        %   A time variant. (default=this.TimeVariant)
        %
        % Returns
        % -------
        % cov: (NumStateDims x NumStateDims) matrix
        %   The process noise covariance matrix
        %
        % Usage
        % -----
        % * Qk = covariance(this,Dt) returns process covariance matrix, upon 
        %   application of Dt.
        % * Qk = covariance(this) returns process covariance matrix, upon 
        %   application of this.TimeVariant. 
        % (i.e. Default Dt = this.TimeVariant)
        %
        % See also FEVAL, RANDOM, PDF.
            switch(nargin)
                case 1 
                    Dt = this.TimeVariant;
            end
            
            % Return process covariance
            cov = this.Q(Dt); % Time variant
        end
        
        function wk = random(this, Ns, Dt)
        % RANDOM Generates random samples from the dynamic model's noise
        % distribution.
        %
        % Parameters
        % ----------
        % Ns: scalar, optional
        %   The number of samples to be generated.
        %   (default = 1)
        % Dt: scalar, optional
        %   A time variant. 
        %   (default=this.TimeVariant)
        %
        % Returns
        % -------
        % wk: (NumStateDims x Ns)
        %   A matrix, whose columns correspond to indivual/independent
        %   noise samples.
        %
        % Usage
        % -----
        % * wk = random(this,Ns,Dt) generates and returns a set of Ns samples
        %   generated from the noise distribution of the dynamic model, i.e.
        %   wk ~ N(0,Q(Dt)), where Q is the noise covariance upon application 
        %   of the time variant Dt.
        % * wk = random(this,Ns) generates and returns a set of Ns samples
        %   generated from the noise distribution of the dynamic model, i.e.
        %   wk ~ N(0,Q), where Q is the noise covariance upon application 
        %   of the time variant this.TimeIndex.
        % (i.e. Default Dt = this.TimeVariant)
        %
        % See also FEVAL, PDF, COVARIANCE
       
            switch(nargin)
                case 1
                    Ns = 1;
                    Dt = this.TimeVariant;
                case 2
                    Dt = this.TimeVariant;
            end
              
            wk = mvnrnd(zeros(this.NumStateDims,Ns)',this.Q(Dt))';
        end
        
        function prob = pdf(this, xk, xkm1, Dt)
        % EVAL Evaluates the probability/likelihood p(x_k|x_{k-1}) of a 
        % (set of) new state vector(s), given a (set of) old state vector(s)  
        %
        % Parameters
        % ----------
        % xk: (NumStateDims x Ns) matrix
        %   A matrix, whose columns correspond to individual new state vectors.
        % xkm1: (NumStateDims x Np) matrix
        %   A matrix, whose columns correspond to individual old state vectors
        % Dt: scalar, optional
        %   A time variant. 
        %   (default=this.TimeVariant)
        %
        % Returns
        % -------
        % prob: (Np x Ns) matrix
        %   A matrix, where each element (i,j) corresponds to the evaluated
        %   probability p(xk(:,j)|xkm1(:,i))
        %
        % Usage
        % -----
        % * prob = pdf(xk, xkm1, Dt) evaluates and returns a (Np x Ns)
        %   probability/likelihood matrix given the (NumStatesDim x Ns) xk
        %   and (NumStatesDim x Np) xkm1 state matrices. Dt is an optional 
        %   argument (Default = this.TimeVariant) time variable, which is 
        %   used for computing the model's covariance.
        %
        % See also FEVAL, RANDOM, COVARIANCE
            
            if(nargin<4)
                Dt = this.TimeVariant;
            end
            
            xk_km1 = this.feval(xkm1,0,Dt);
            prob = zeros(size(xk,2), size(xkm1,2));
            if(size(xkm1,2)>size(xk,2))
                for i=1:size(xk,2)
                    prob(i,:) = gauss_pdf(xk(:,i), xk_km1, this.Q(Dt));
                end
            else
                for i=1:size(xkm1,2)
                    prob(:,i) = gauss_pdf(xk, xk_km1(:,i), this.Q(Dt))';  
                end
            end
                        
        end
        
    end
end