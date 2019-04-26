classdef GaussianMixtureStateX < GaussianMixtureX & StateX
% GaussianMixtureStateX class
%
% Summary of GaussianMixtureStateX:
% Class implementation of the primitive Gaussian Mixture State type. 
% GaussianMixtureStateX is essentially a StateX that comes pre-initialised with a Gaussian
% Distribution.
% 
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    methods
        function this = GaussianMixtureStateX(varargin)
        % GaussianMixtureStateX Construct a multivariate Gaussian 
        %   distribution
        %
        % Parameters
        % ----------
        % components: (1 x NumComponents) structure array
        %   An array of component structures, where each structure contains
        %   the following 2-3 fields:
        %       - Mean: (NumVariables x 1) column vector
        %       - Covar: (NumVariables x NumVariables) matrix
        %       - Weight: scalar
        % timestamp: datetime, optional
        %   The state timestamp
        %                       
        % -- OR ---
        %
        % means: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covars: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % weights: (1 x NumVariables) row vector, optional
        %   The weights for all Mixture components
        %   Default is 1/NumComponents for all components
        % timestamp: datetime, optional
        %   The state timestamp
        
            [timestamp, other] = StateX.extract_timestamp(varargin);
            
            this@StateX(timestamp);
            this@GaussianMixtureX(other{:});
        end
    end
    methods (Access=protected)
        function stateVector = getVector(this)
            stateVector = this.Means;
        end
        function Vector = setVector(this, stateVector)
            error("Setting the Vector property is not currently supported for GaussianMixtureStateX objects!");
        end
        function numDims = getNumDims(this)
            numDims = this.NumVariables;
        end
    end
end

