classdef GaussianStateX < GaussianDistributionX & StateX
% GaussianStateX class
%
% Summary of GaussianStateX:
% Class implementation of the primitive Gaussian State type. GaussianStateX
% is essentially a StateX that comes pre-initialised with a Gaussian
% Distribution.
%
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    methods
        function this = GaussianStateX(varargin)
        % GaussianStateX Construct a multivariate Gaussian state
        %
        % Parameters
        % ----------
        % numVariables: scalar, optional
        %   The number of random variables. Onlt necessary if neither mean
        %   nor covar specified.
        % timestamp: datetime, optional
        %   The state timestamp
        %
        % -- OR ---
        %
        % distribution: ProbabilityDistributionX
        %   Another distribution, based on which the new Gaussian
        %   distribution should be created.
        % timestamp: datetime, optional
        %   The state timestamp
        %
        % -- OR ---
        %
        % mean: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % timestamp: datetime, optional
        %   The state timestamp
            
            [timestamp, other] = StateX.extract_timestamp(varargin);
            
            this@StateX(timestamp);
            this@GaussianDistributionX (other{:});
        end
    end
    
    methods (Access=protected)
        function stateVector = getVector(this)
            stateVector = this.Mean;
        end
        function Vector = setVector(this, stateVector)
            Vector = stateVector;
            this.Mean = stateVector;
        end
        function numDims = getNumDims(this)
            numDims = this.NumVariables;
        end
    end

end

