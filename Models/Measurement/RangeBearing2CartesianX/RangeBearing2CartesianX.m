classdef RangeBearing2CartesianX < MeasurementModelX
    % RangeBearing2Cartesian class
    %
    % Summary of RangeBearing2CartesianX
    % This is a class implementation of a time & state invariant Polar-to-Cartesian
    % Position Measurement Model, with additive Gaussian noise.
    %
    % The model is described by the following equations:
    %   
    %   y = h(x) + v,  v ~ N(0,R)
    % 
    % where:
    %   
    %   y = [ bearing ; 
    %         range ]
    %
    % LinearGaussianX Properties:
    %   + NumStateDims           - The number of state dimensions.
    %   + NumMeasDims            - The number of measurement dimensions.
    %   + MeasurementErrVariance - Measurement error variance on each co-ordinate
    %
    % LinearGaussianX Methods:
    %   + feval     - Equivalent to applying the measurement model function 
    %   + random    - Measurement noise sample generator function, i.e. v_k ~ random(~)
    %   + pdf       - Function to evaluate the probability p(y_t|x_t) of 
    %                 a (set of) measurements, given a (set of) state vector(s)
    %                   e.g. pdf = @(yt,xt) mvnpdf(yt,xt,R);
    %   + matrix    - Returns the model transformation matrix
    %   + covar     - Returns the measurement noise covariance matrix 
    %
    % (+) denotes puplic properties/methods
    %
    %  February 2018 Lyudmil Vladimirov, University of Liverpool   

    properties (Access = private)
        R
    end
        
    properties
        NumStateDims = 2
        MeasurementErrVariance
        Mapping = [1 2]
    end
    
    properties (Dependent)
      NumMeasDims
    end
    
    methods (Access = private)
        function yk = h(this,xk)
           yk = zeros(2,size(xk,2));
           [yk(1,:),yk(2,:)] = cart2pol(xk(this.Mapping(1),:),xk(this.Mapping(2),:));
        end
        
        function xk = h_inv(this,yk)
           xk = zeros(this.NumStateDims,size(yk,2));
           [xk(this.Mapping(1),:),xk(this.Mapping(2),:)] = pol2cart(yk(1,:),yk(2,:));
        end
        
        function initialise_(this, config)
            if(isfield(config,'NumStateDims'))
                this.NumStateDims = config.NumStateDims;
            end
            if(isfield(config,'MeasurementErrVariance'))
                this.MeasurementErrVariance = config.MeasurementErrVariance;
                [numRows,numCols] = size(this.MeasurementErrVariance);
                if(numRows==1 || numCols==1)
                    len = length(this.MeasurementErrVariance);
                    if(len ~= this.NumMeasDims)
                        error("Invalid MeasurementErrVariance size! Type 'help RangeBearing2CartesianX/RangeBearing2CartesianX' to get more information");
                    end        
                    this.R = diag(this.MeasurementErrVariance);
                else
                    if(numRows == this.NumMeasDims && numCols==this.NumMeasDims)
                        this.R = this.MeasurementErrVariance;
                    end
                end
            end
            if(isfield(config,'Mapping'))
                this.Mapping = config.Mapping;
            end
        end
    end
        
    methods
        function this = RangeBearing2CartesianX(varargin)
        % RangeBearing2CartesianX Constructor method
        % 
        % Parameters
        % ----------
        % NumStateDims: scalar
        %   The number of state dimensions (default: 1)
        % Mapping: row vector 
        %   A (1 x 2) row vector used to specify the index of the 
        %   state dimension which contributes to each measurement dimension
        % MeasurementErrVariance: scalar, vector or matrix
        %   Measurement error (co)variance, that can take the following values:
        %   * Specifying MeasurementErrVariance as a scalar will produce a noise 
        %     covariance matrix whose diagonal elements are equal to the value 
        %     of MeasurementErrVariance, with the off-diagonals equal to 0;
        %   * Specifying MeasurementErrVariance as a vector of length NumMeasDims, 
        %     will place all elements of the vector on the respective diagonals 
        %     of the noise covariance, with the off-diagonals equal to 0;
        %   * Finally, if MeasurementErrVariance is specified as a matrix 
        %     of shape (NumMeasDims x NumMeasDims) then its value will be used 
        %     as the full noise covariance.
        % Usage
        % ----- 
        % * RangeBearing2CartesianX(___,Name,Value) instantiates an object 
        %   handle, configured the parameters specified by one or
        %   more Name,Value pair arguments. 
        % * RangeBearing2CartesianX(config) instantiates an object handle   
        %   configured with the parameters specified inside the 'config' 
        %   structure, whose  fieldnames correspond to the parameter names 
        %   as given in the Parameters section above.
        %
        % See also feval, random, pdf, covariance.   
        
            % Call SuperClass method
            this@MeasurementModelX();
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
                  
        end
        
        function yk = feval(this, xk, vk)
        % feval Project a given state vector on a 2D polar plane  
        %
        % Parameters
        % ----------
        % xk: matrix, optional
        %   - A (NumStateDims x Ns) matrix, where each column corresponds to
        %     a (NumStateDims x 1) state vector. 
        % vk: (NumMeasDims x Ns) matrix or boolean, optional
        %   - If vk is a boolean and set True, then the default model noise 
        %     generation function (this.random()) will be used to generate the 
        %     noise. 
        %   - Otherwise, vk should be a matrix of equal number of columns to xk, 
        %     where each column corresponds to the random noise which will
        %     be added to each transformed state vector. 
        %   - If not provided, then no noise will be added to the state vectors.
        %
        % Returns
        % -------
        % yk: (NumMeasDims x Ns) or (NumMeasDims x NumMeasDims) matrix
        %   - If no parameters are passed to the function, then yk will be the
        %     (NumMeasDims x NumMeasDims) measurement matrix.
        %   - Else, yk will be a (NumMeasDims x Ns) matrix, where each column 
        %     corresponds to a measurement projection of the respective column in xk.
        %
        % See also pdf, random, covar.
            
            switch(nargin)
                case 2
                    vk = 0;
            end
            
            if(islogical(vk) && vk)
                vk = this.random(size(xk,2));
            end
            
            % Compute predicted measurement
            
            yk = this.h(xk) + vk ;             
        end
        
        function xk = finv(this, yk, vk)
        % finv Project a given measurement vector back into Cartesian
        %
        % Parameters
        % ----------
        % yk: (NumMeasDims x Ns) matrix
        %   - If no parameters are passed to the function, then yk will be the
        %     (NumMeasDims x NumMeasDims) measurement matrix.
        %   - Else, yk will be a (NumMeasDims x Ns) matrix, where each column 
        %     corresponds to a measurement projection of the respective column in xk.
        % vk: (NumMeasDims x Ns) matrix , optional
        %   - vk should be a matrix of equal dimensions to yk, where each 
        %     column corresponds to the random noise which will be substracted
        %     from each measurement vector (i.e. yk-vk) prior to projecting back. 
        %   - If not provided, then no noise will be considered.
        %
        % Returns
        % -------
        % xk: matrix, optional
        %   - A (NumStateDims x Ns) matrix, where each column corresponds to
        %     a (NumStateDims x 1) state vector. 
        %   - If not provided, then the measurement matrix (H) will be returned.
        %
        %   See also OBS_COV, OBS_NOISE, SAMPLE, EVAL.
            
            switch(nargin)
                case 2
                    vk = 0;
            end
            
            % Compute predicted measurement
            xk = this.h_inv(yk - vk) ;           
        end
        
        function Rk = covar(this)
        % covar Returns measurement noise covariance matrix.
        %
        % Returns
        % -------
        % R: matrix
        %   The (NumMeasDims x NumMeasDims) noise covariance matrix.
        %
        % See also feval, random, pdf.
                    
            % Return measurement covariance
            Rk = this.R;  
        end
        
        function vk = random(this, Ns)
        % random Generates random samples from the measurement model's noise
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
        %   A (NumMeasDims x Ns) matrix, where each column corresponds to a
        %   sample drawn from the measurement model's noise distribution.
        %
        % See Also
        % --------
        % feval, pdf, covar
            
            if nargin<2
                Ns = 1;
            end
        
            % Generate noise samples
            vk = mvnrnd(zeros(this.NumMeasDims,1)',this.R, Ns)';
        end
        
        function prob = pdf(this, yk, xk, Pk)
        % pdf Evaluates the probability/likelihood N(yk | xk, Pk) between a 
        % (set of) measurement(s) and a (set of) state vector(s)  
        %
        % Parameters
        % ----------
        % yk: matrix
        %   A (NumMeasDims x Nm) matrix, whose columns correspond to
        %   Nm individual measurements
        % xk: matrix
        %   A (NumStateDims x Ns) matrix, whose columns correspond to Ns
        %   individual state vectors
        % Pk: matrix, optional
        %   A (NumStateDims x NumStateDims) or (NumStateDims x NumStateDims x Ns)
        %   covariance matrix/matrices.
        %   - If Pk is a (NumStateDims x NumStateDims x Ns) matrix, then
        %     each (NumStateDims x NumStateDims x i) matrix corresponds to 
        %     the covariance used when evaluating of sample xk(:,i).
        %   - If Pk is a (NumStateDims x NumStateDims) matrix, then it is
        %     implied that all samples in xk have the same covariance.
        %   Note that the measurement noise covariance will be added to the
        %   transformed covariance matrix prior to evaluation.
        %   (default Pk = 0, which implies that only the measurement noise
        %    will be taken into account)
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
        % See also feval, random, covar
            
            Ns = size(xk,2);
            Nm = size(yk,2);
            if(nargin<4)
                Pk = 0;
            end
            if(size(Pk,3)>1 && size(Pk,3)~=Ns)
                error("Pk size error! Check function documentation for more details");
            end
            prob = zeros(Nm, Ns);
            % Increase speed by iterating over the smallest matrix 
            if(Ns>Nm && size(Pk,3)==1)
                %Sk = this.H*Pk*this.H' + this.R;
                for i=1:size(yk,2)
                    prob(i,:) = gauss_pdf(yk(:,i), this.feval(xk), this.R);
                end
            else
                for i=1:size(xk,2)
                    %Sk = this.H*Pk(:,:,i)*this.H' + this.R;
                    prob(:,i) = gauss_pdf(yk, this.feval(xk(:,i)), this.R);  
                end
            end
        end
        
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        
        function numMeasDims = get.NumMeasDims(this)
            numMeasDims = 2;
        end
        
    end
end