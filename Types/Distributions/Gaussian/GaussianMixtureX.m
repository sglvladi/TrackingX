classdef GaussianMixtureX < ProbabilityDistributionX
% GaussianMixtureX class
%
% Summary of GaussianMixtureX:
% This is a class implementation of a Multivariate Gaussian distribution
%
% GaussianMixtureX Properties:
%   + Means         - Mean vectors of all Gaussian Mixture Components
%   + Covars    	- The covariance matrices of all Gaussian Mixture Copmponents
%   + Weights       - The Mixture weights of all Gaussian Mixture Components
%   + Components    - Cell array of Mixture Component objects
%   + NumVariables  - The number of random variables in each of the
%   + NumComponents - The number of mixture components in the Gaussian
%
% GaussianMixtureX Methods:
%   + GaussianMixtureX  - Constructor method
%   + reset - Reset the distribution with a new number of random variables
%   + random - Draw random samples from the multivariate normal distribution
%   + pdf - Evaluate the density of the multivariate normal distribution
%   + fitToSamples - Fit a normal distribution to a set of weighted samples
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties (Dependent)
        % Means: (NumVariables x NumComponents) matrix
        %   Mean vectors of the Gaussian Mixture Components
        Means
        
        % Covars: (NumVariables x NumVariables x NumComponents) matrix
        %   Covariance matrices of the Gaussian Mixture Components
        %   All covariance matrices must be symmetric and positive semi-definite.
        Covars
        
        % Weights: (1 x NumComponents) row vector
        %   Mixture Weights of the Gaussian Mixture Components
        Weights
        
        % NumComponents: scalar
        %   The number of mixture components in the Gaussian Mixture 
        NumComponents
        
        % Components: (1 x NumComponents) structure array
        %   A structure array whose elements represent individual mixture
        %   components
        Components = {};
    end
    
    properties   
        Means_ = 0
        Covars_ = 1
        Weights_ = 1
        NumComponents_ = 1
    end
    
    methods
        function this = GaussianMixtureX(varargin)
        % MultivariateNormalDistribution Construct a multivariate Gaussian 
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
        %                       
        %                     - OR -
        %
        % means: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covars: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % weights: (1 x NumVariables) row vector, optional
        %   The weights for all Mixture components
        %   Default is 1/NumComponents for all components
            
            this.reset(varargin{:});
        end
        
        function samples = random(this, numSamples)
        % Samples Draw random samples from the Gaussian Mixture
        %
        % Parameters
        % ----------
        % numSamples: scalar
        %   The number of samples to be drawn from the mixture
        %
        % Returns
        % -------
        % samples: (NumVariables x numSamples) matrix
        %   The set of samples drawn from the mixture.
            
            assert(numSamples >= 1);
            gm = gmdistribution(this.Means',this.Covars,this.Weights);
            samples = gm.random(numSamples)';
        end  
        
        function prob = pdf(this, samples)
        % pdf Evaluate the density of the multivariate normal distribution
        %
        % Parameters
        % ----------
        % samples: (NumVariables x numSamples) matrix
        %   
        % Returns
        % -------
        % prob: (1 x numSamples) matrix
        %   The evaluated density for the corresponding set of samples
            
            gm = gmdistribution(this.Means',this.Covars,this.Weights);
            prob = gm.pdf(samples')';
        end
        
        function varargout = cluster(this, data)
        % cluster Cluster data for Gaussian mixture distribution.
        %
        %   Essentially a re-implementation of 'gmdistribution/cluster',
        %   where data = X'. Type 'help gmdistribution/cluster' for more
        %   details.
        %
        % Parameters
        % ----------
        % data: (NumVariables x NumData) matrix
        %
        % Returns
        % -------
        % idx: (NumData x 1) column vector
        %   Cluster index, where idx(i) is the cluster index of observation 
        %   i and indicates the Gaussian mixture component with the largest 
        %   posterior probability given the observation i.
        % nlogL: scalar
        %   Negative loglikelihood value of the Gaussian mixture model given 
        %   the data, returned as a numeric value.
        % P: (NumData x NumComponents) matrix
        %   Posterior probability of each Gaussian mixture component given 
        %   each observation in the data.
        % logpdf: (NumData x 1) vector
        %   Logarithm of the estimated pdf, evaluated at each observation in
        %   the data.
        % d2: (NumData x NumComponents)
        %   Squared Mahalanobis distance of each observation in data to each 
        %   Gaussian mixture component. 
        
            gm = gmdistribution(this.Means',this.Covars,this.Weights);
            [varargout{1:nargout}] = gm.cluster(data);
        end
        
        function d2 = mahal(this, data)
        % mahal Mahalanobis distance to Gaussian mixture component
        %
        %   Essentially a re-implementation of 'gmdistribution/mahal',
        %   where data = X'. Type 'help gmdistribution/mahal' for more
        %   details.
        %
        % Parameters
        % ----------
        % data: (NumVariables x NumData) matrix
        %   Input data
        %
        % Returns
        % -------
        % d2: (NumData x NumData) matrix
        %   Squared Mahalanobis distance of each observation in data to 
        %   each Gaussian mixture component
        %
            gm = gmdistribution(this.Means',this.Covars,this.Weights);
            d2 = gm.mahal(data');
        end
        
        function reset(this, varargin)
        % reset Reset the distribution with a new number of random variables
        % 
        % Parameters
        % ----------
        % components: (1 x NumComponents) structure array
        %   An array of component structures, where each structure contains
        %   the following 2-3 fields:
        %       - Mean: (NumVariables x 1) column vector
        %       - Covar: (NumVariables x NumVariables) matrix
        %       - Weight: scalar
        %
        %                  - OR -
        %
        % means: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covars: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % weights: (1 x NumVariables) row vector, optional
        %   The weights for all Mixture components
        %   Default is 1/NumComponents for all components
            
            switch(nargin)
                case(2)
                    if (iscell(varargin{1}))
                        components = varargin{1};
                        this.Means_ = cell2mat(cellfun(@(c)c.Mean,components,'UniformOutput',false));
                        this.Covars_ = cell2mat(permute(cellfun(@(c)c.Covar,components,'UniformOutput',false),[1,3,2]));
                        this.Weights_ = cell2mat(cellfun(@(c)c.Weight,components,'UniformOutput',false));
                    end
                otherwise
                    this.Means = varargin{1};
                    this.Covars = varargin{2};
                    [numVariables,numComponents] = size(this.Means);
                    if(nargin>3)
                        this.Weights = varargin{3};
                    else
                        this.Weights = repmat(1/numComponents,1,numComponents);
                    end
                    this.NumVariables_ = numVariables;
                    this.NumComponents_ = numComponents;
            end
        end
        
        function this = plus(this,other)
            this.Means_ = [this.Means_, other.Means];
            this.Covars_ = cat(3, this.Covars_, other.Covars);
            this.Weights_ = [this.Weights_, other.Weights];
        end
    end
    
    methods
        
        function components = get.Components(this)
            numComponents = this.NumComponents;
            components = cell(1,numComponents);
            for i = 1:numComponents
                components{i}.Mean = this.Means(:,i);
                components{i}.Covar = this.Covars(:,:,i);
                components{i}.Weight = this.Weights(:,i);
            end
        end
        
        function means = get.Means(this)
        %get.Means Getter for Means property
            means = this.Means_;
        end
        
        function set.Means(this, means)
        %set.Means - Setter for Mean property
            this.Means_ = means;
        end
        
        function covars = get.Covars(this)
        %get.Covars Getter for Covars property
            covars = this.Covars_;
        end
        
        function set.Covars(this, covars)
        %set.Covars Setter for Covars property
            this.Covars_ = covars;
        end
        
        
        function weights = get.Weights(this)
        %get.Means Getter for Means property
            weights = this.Weights_;
        end
        
        function set.Weights(this, weights)
        %set.Means - Setter for Mean property
            this.Weights_ = weights;
        end
        
        function numComponents = get.NumComponents(this)
            numComponents = size(this.Means_,2);
        end
    end
    
end

