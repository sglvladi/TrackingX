classdef UniformDistributionX < ProbabilityDistributionX
% UniformDistributionX class
% 
% Summary of UniformDistributionX:
%   This is a class implementation of a Multivariate Uniform distribution
%
% UniformDistributionX Properties:
%   + Limits - Lower and upper limits for each of the random variables
%
% UniformDistributionX Properties:
%   + UniformDistributionX - Constructor method
%   + reset - Reset the distribution with a new number of random variables
%   + random - Draw random samples from the multivariate uniform distribution
%   + pdf - Evaluate the density of the multivariate uniform distribution
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties (Dependent)
        % Limits - Lower and upper limits for each of the random variables
        %    Each row of this NumRandomVariables-by-2 array corresponds
        %    to the lower and upper limits of a single random variable. The
        %    limits define an open interval.
        %    If the underlying distribution is k-variate, the array has k
        %    rows.
        %
        %    Default: [zeros(NumRandomVariables, 1) ones(NumRandomVariables, 1)]
        Limits
    end
    
    properties (Access = private)
        % Limits_ - Internal storage for random variable limits
        Limits_ = [0 1];
    end
    
    methods
        function this = UniformDistributionX(varargin)
        % MultivariateUniformDistribution Construct a multivariate Uniform 
        % distribution
        %
        % Parameters
        % ----------
        % Limits: (NumVariables x 2) matrix, optional
        %   Lower and upper limits for each of the random variables, where
        %   each row corresponds to a new variable
        %   (default: [0,1])
            this.reset(varargin{:});
        end
        
        function samples = random(this, numSamples)
        % Samples Draw random samples from the Uniform Distribution
        %
        % Parameters
        % ----------
        % numSamples: scalar
        %   The number of samples to be drawn from the distribution
        %
        % Returns
        % -------
        % samples: (NumVariables x numSamples) matrix
        %   The set of samples drawn from the mixture.
        
            assert(numSamples >= 1);
            samples = unifrnd(repmat(this.Limits(:,1),1,numSamples), ... 
                              repmat(this.Limits(:,2),1,numSamples), ...
                              this.NumVariables, numSamples);
        end
        
        function prob = pdf(this, samples)
        % pdf Evaluate the density of the multivariate Uniform distribution
        %
        % Parameters
        % ----------
        % samples: (NumVariables x NumSamples) matrix
        %   
        % Returns
        % -------
        % prob: (1 x numSamples) matrix
        %   The evaluated density for the corresponding set of samples
            
            numSamples = size(samples,2);
            assert(numSamples >= 1);
            prob = prod(unifpdf(samples,this.Limits(:,ones(1,numSamples)),this.Limits(:,2*ones(1,numSamples))),1);
        end
        
        function reset(this, varargin)
        %reset Reset the distribution with a new number of random variables
        %
        % Parameters
        % ----------
        % Limits: (NumVariables x 2) matrix, optional
        %   Lower and upper limits for each of the random variables, where
        %   each row corresponds to a new variable
        %   (default: [0,1])
            
            if(nargin>1)
                this.Limits = varargin{1};
            end
            this.NumVariables = size(this.Limits,1);
        end
        
        function limits = get.Limits(this)
            %get.Limits Getter for Limits property
            limits = this.Limits_;
        end
        
        function set.Limits(this, limits)
            %set.Limits Setter for Limits property
            
            % Assert that lower limits are all <= to upper limits
            assert(all(limits(:,1) <= limits(:,2)));
            
            this.Limits_ = double(limits);
            this.NumVariables = size(this.Limits_,1);
        end
    end    
end

