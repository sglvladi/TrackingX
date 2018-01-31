classdef LinGaussObsModelX_2D < ObservationModelX
% LinGaussObsModelX_2D class
%
% Summary of LinGaussObsModelX_2D
% This is a class implementation of a time-invariant 2D Linear-Gaussian 
% Observation Model.
%
% The model is described by the following equations:
%   y(t) = H*x(t-1) + v(t),  v(t)~ N(0,R)
%
% LinGaussObsModelX_2D Properties:
%   - ObsErrVariance  Measurement error variance on each co-ordinate
%
% LinGaussObsModelX_2D Methods:
%   - feval(~)         Equivalent to applying the model transition equations 
%   - random(~)        Measurement noise sample generator function
%   - pdf(~)           Function to evaluate the probability/likelihood 
%                      p(y_k|x_k) of between a (set of) measurement(s)and a
%                      (set of) state vector(s) 
%   - covariance(~)    Returns the measurement noise covariance matrix 
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool   

    properties (Access = private)
        H
        R
    end
        
    properties
        Mapping = [1 3]
        ObsErrVariance
    end
    
    methods
        function this = LinGaussObsModelX_2D(varargin)
        % LINGAUSSOBSMODEL_2D Constructor method
        %   
        % DESCRIPTION: 
        % * LinGaussObsModelX_2D(xDim,r) instantiates an object handle 
        %   configured with the provided number of state dimentions xDim
        %   and the measurement noise variance r.
        % * LinGaussObsModelX_2D(___,Mapping) instantiates an object handle 
        %   where Mapping is a (1 x 2) row vector used to specify the index
        %   of the state dimension to which each measurement dimension relates.
        %   e.g. Assume we a 4 dimentional state vector x = [a; b; c; d],
        %   a measurement y = [z;v] and we want to define the mapping such
        %   that z -> d and v -> a. Then Mapping = [4,1], i.e. the first
        %   index of the measurement vector (z) relates to the 4th index of 
        %   the state vector (d) and so on... 
        % * LinGaussObsModelX_2D(config) instantiates an object handle 
        %   configured with options specified as fields under the config
        %   structure. 
        %   e.g. config.NumStateDims, config.ObsErrVariance, config.Mapping
        % * ConstantVelocityModel_2D(___,Name,Value,___) instantiates an object 
        %   handle, configured with options specified by one or more
        %   Name,Value pair arguments.
        %
        %  See also feval, random, pdf, covariance.   
        
            % Call SuperClass method
            this@ObservationModelX();
            
            this.NumObsDims = 2;
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.NumStateDims = config.NumStateDims;
                    this.ObsErrVariance = config.ObsErrVariance;
                    if(isfield(config,'Mapping'))
                        this.Mapping = config.Mapping;
                    end
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addOptional('NumStateDims',[]);
            parser.addOptional('ObsErrVariance',[]);
            parser.addOptional('Mapping',[]);
            parser.parse(varargin{:});
            
            this.NumStateDims =  parser.Results.NumStateDims;
            this.ObsErrVariance =  parser.Results.ObsErrVariance;
            if(~isempty(parser.Results.Mapping))
                this.Mapping = parser.Results.Mapping;
            end
                       
            % Define .h
            H = zeros(2, this.NumStateDims);
            if(this.Mapping(1)~=0)
                H(1,this.Mapping(1)) = 1;
            end
            if(this.Mapping(2)~=0)
                H(2,this.Mapping(2)) = 1;
            end
            this.H = H;
          
            this.R = [this.ObsErrVariance, 0; 0, this.ObsErrVariance];
      
        end
        
        function yk = heval(this, xk, vk)
        % HEVAL Project a given state vector into the measurement space  
        %
        % DESCRIPTION:
        % * yk = HEVAL(this,xk,vk) returns the vector yk, produced by
        %   propagating the given state xk through the observation model
        %   and adding random noise vk. xk must be a (NumStateDims x Ns) 
        %   matrix and vk must be a (NumObsDims x Ns) matrix,  where Ns  
        %   is the number of states/columns in xk.
        % * yk = HEVAL(this,xk,vk) returns the vector yk, produced by
        %   propagating the given state xk through the observation model
        %   without the addition of noise.
        %   (i.e. Default wk = 0)
        % * yk = HEVAL(this) returns the model's transition matrix 
        %
        % See also PDF, RANDOM, COVARIANCE.
            
            switch(nargin)
                case 1
                    xk = 1;
                    vk = 0;
                case 2
                    vk = 0;
            end
            
            if(islogical(vk) && vk)
                vk = this.random(size(xk,2));
            end
            
            % Compute predicted measurement
            yk = this.H*xk + vk ;           
        end
        
        function Rk = covariance(this)
        % COVARIANCE Returns measurement noise covariance matrix.
        %
        % DESCRIPTION:
        % * R = covariance(this) returns the measurement noise covariance 
        %   matrix
        %
        % See also HEVAL, RANDOM, PDF.
                    
            % Return measurement covariance
            Rk = this.R;  
        end
        
        function vk = random(this, Ns)
        % RANDOM Generates random samples from the observation model's noise
        % distribution.
        %
        % DESCRIPTION:
        % * vk = random(this,Ns) generates and returns a set of Ns samples
        %   generated from the noise distribution of the observation model, 
        %   i.e. vk ~ N(0,R), where R is the noise covariance of the model.
        %   (Default Ns = 1)
        % See also HEVAL, PDF, COVARIANCE
            
            if nargin<2
                Ns = 1;
            end
        
            % Generate noise samples
            vk = mvnrnd(zeros(this.NumObsDims,Ns)',this.R)';
        end
        
        function prob = pdf(this, yk, xk)
        % PDF Evaluates the probability/likelihood p(y_k|x_k) between a 
        % (set of) measurement(s) and a (set of) state vector(s)  
        % 
        % DESCRIPTION:
        % * prob = pdf(y_k, x_k) evaluates and returns a (a x b)
        %   probability/likelihood matrix given the (NumObsDim x a) y_k
        %   and (NumStatesDim x b) x_k state matrices.
        %
        % See also HEVAL, RANDOM, COVARIANCE
        
            prob = zeros(size(yk,2), size(xk,2));
            % Increase speed by iterating over the smallest matrix 
            if(size(xk,2)>size(yk,2))
                for i=1:size(yk,2)
                    prob(i,:) = mvnpdf(yk(:,i)', this.heval(xk)', this.R)';
                end
            else
                for i=1:size(xk,2)
                    prob(:,i) = mvnpdf(yk', this.heval(xk(:,i))', this.R)';  
                end
            end
        end
    end
end