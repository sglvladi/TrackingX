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
% where:
%
%   y(t) = [x_pos; % Position in X coordinate
%           y_pos] % Position in Y coordinate
%
% LinGaussObsModelX_2D Properties:
%   + NumStateDims - The number of state dimensions.
%   + NumObsDims - The number of observation dimensions.
%   + ObsErrVariance - Measurement error variance on each co-ordinate
%
% LinGaussObsModelX_2D Methods:
%   + heval - Equivalent to applying the measurement model function 
%   + random - Measurement noise sample generator function, i.e. v_k ~ random(~)
%   + pdf - Function to evaluate the probability p(y_t|x_t) of 
%                a (set of) measurements, given a (set of) state vector(s)
%                   e.g. pdf = @(yt,xt) mvnpdf(yt,xt,R);
%   + covariance - Returns the measurement noise covariance matrix 
%
% (+) denotes puplic properties/methods
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
        % Parameters
        % ----------
        % NumStateDims: scalar
        %   The number of state dimensions
        % ObsErrVariance: scalar
        %   Measurement error variance on each co-ordinate
        % Mapping: row vector 
        %   A (1 x 2) row vector used to specify the index of the state 
        %   dimension to which each measurement dimension relates
        %  
        % Usage
        % -----
        % * LinGaussObsModelX_2D(NumStateDims,ObsErrVariance,___) instantiates 
        %   an object handle configured with the provided number of state 
        %   dimensions NumStateDims and the measurement noise variance ObsErrVariance.
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
        %   structure. Each field of the provided structure should have the
        %   same name with the corresponding method parameter.
        %   e.g. config.NumStateDims, config.ObsErrVariance, config.Mapping
        % * ConstantVelocityModel_2D(___,Name,Value,___) instantiates an object 
        %   handle, configured with options specified by one or more
        %   Name,Value pair arguments.
        %
        %  See also heval, random, pdf, covariance.   
        
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
        % Parameters
        % ----------
        % xk: matrix, optional
        %   - A (NumStateDims x Ns) matrix, where each column corresponds to
        %     a (NumStateDims x 1) state vector. 
        %   - If not provided, then the measurement matrix (H) will be returned.
        % vk: (NumObsDims x Ns) matrix or boolean, optional
        %   - If vk is a boolean and set True, then the default model noise 
        %     generation function (this.random()) will be used to generate the 
        %     noise. 
        %   - Otherwise, vk should be a matrix of equal dimensions to xk, 
        %     where each column corresponds to the random noise which will
        %     be added to each state vector. 
        %   - If not provided, then no noise will be added to the state vectors.
        %
        % Returns
        % -------
        % yk: (NumObsDims x Ns) or (NumObsDims x NumObsDims) matrix
        %   - If no parameters are passed to the function, then yk will be the
        %     (NumObsDims x NumObsDims) measurement matrix.
        %   - Else, yk will be a (NumObsDims x Ns) matrix, where each column 
        %     corresponds to a measurement projection of the respective column in xk.
        %
        % Usage
        % -----
        % * yk = HEVAL(this,xk,vk) returns the vector yk, produced by
        %   propagating the given state xk through the observation model
        %   and adding random noise vk. xk must be a (NumStateDims x Ns) 
        %   matrix and vk must be a (NumObsDims x Ns) matrix,  where Ns  
        %   is the number of states/columns in xk. If vk is a logical (i.e.
        %   boolean) and set True, then the default model noise generation
        %   function (this.random()) will be used to generate the noise.
        % * yk = HEVAL(this,xk) returns the vector yk, produced by
        %   propagating the given state xk through the observation model
        %   without the addition of noise.
        %   (i.e. Default vk = False)
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
        % Returns
        % -------
        % R: matrix
        %   The (NumObsDims x NumObsDims) noise covariance matrix.
        %
        % Usage
        % -----
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
        % Parameters
        % ----------
        % Ns: scalar, optional
        %   Number of samples to generate 
        %   (default = 1)
        %
        % Returns
        % -------
        % vk: matrix
        %   A (NumObsDims x Ns) matrix, where each column corresponds to a
        %   sample drawn from the observation model's noise distribution.
        %
        % Usage
        % -----
        % * vk = random(this,Ns) generates and returns a set of Ns samples
        %   generated from the noise distribution of the observation model, 
        %   i.e. vk ~ N(0,R), where R is the noise covariance of the model.
        %   (default = 1)
        %
        % See Also
        % --------
        % HEVAL, PDF, COVARIANCE
            
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
        % Parameters
        % ----------
        % yk: matrix
        %   A (NumObsDims x Nm) matrix, whose columns correspond to
        %   Nm individual measurements
        % xk: matrix
        %   A (NumStateDims x Ns) matrix, whose columns correspond to Ns
        %   individual state vectors
        %
        % Returns
        % -------
        % prob: matrix
        %   A (Nm x Ns) matrix, where :math:`prob(i,j)=p(y_k^i|x_k^j) 
        %
        % Usage
        % -----
        % * prob = pdf(y_k, x_k) evaluates and returns a (a x b)
        %   probability/likelihood matrix given the (NumObsDim x a) y_k
        %   and (NumStatesDim x b) x_k state matrices.
        %
        % See also HEVAL, RANDOM, COVARIANCE
        
            prob = zeros(size(yk,2), size(xk,2));
            % Increase speed by iterating over the smallest matrix 
            if(size(xk,2)>size(yk,2))
                for i=1:size(yk,2)
                    prob(i,:) = gauss_pdf(yk(:,i), this.heval(xk), this.R);
                end
            else
                for i=1:size(xk,2)
                    prob(:,i) = gauss_pdf(yk, this.heval(xk(:,i)), this.R);  
                end
            end
        end
    end
end