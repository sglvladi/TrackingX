classdef PoissonDistributionX < ProbabilityDistributionX
% PoissonDistributionX class
%
% Summary of PoissonDistributionX:
% This is a class implementation of a Poisson distribution
%
% PoissonDistributionX Properties:
%   + Mean - Mean of the Poisson distribution
%   + Covar - Variance of the Poisson distribution
%
% PoissonDistributionX Methods:
%   + PoissonDistributionX  - Constructor method
%   + reset - Reset the distribution with a new number of random variables
%   + random - Draw random samples from the multivariate normal distribution
%   + pdf - Evaluate the density of the multivariate normal distribution
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties (Dependent)
        % Mean: scalar 
        %   Mean of the Poisson distribution
        Mean
        
        % Covar: scalar 
        %   Variance of the Poisson distribution
        Covar
    end
    
    properties (Access = private)
        % Underlying poisson distribution
        pd_ = makedist('Poisson');
    end
    
    methods
        function this = PoissonDistributionX(varargin)
        % PoissonDistributionX Construct a Poisson distribution
        %
        % Parameters
        % ----------
        % mean: non-negative scalar
        %   Mean of the Poisson distribution
            this.reset(varargin{:});
        end
        
        function samples = random(this, numSamples)
        % Samples Draw random samples from the Gaussian 
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

            % Incorporate actual mean
            samples = this.pd_.random(numSamples,1)';
        end  
        
        function prob = pdf(this, samples)
        % pdf Evaluate the density of the poisson distribution
        %   PROB = PDF(OBJ, SAMPLES) returns the density of the poisson distribution,
        %   evaluated at each column of SAMPLES.
        %
            prob = this.pd_.pdf(samples);
        end
        
        function reset(this, varargin)
            if(nargin >= 2)
                mean = varargin{1};

                % Reset number of random variables. This will also take care of
                % input validation.
                this.NumVariables = size(mean,1);

                % Re-initialize the mean and covariance            
                this.Mean = mean;
            end
        end
    end
    
    methods
        function meanValue = get.Mean(this)
        %get.Mean Getter for Mean property
            meanValue = this.pd_.mean();
        end
        
        function set.Mean(this, meanValue)
        %set.Mean - Setter for Mean property
            this.pd_ = makedist('Poisson','lambda',  meanValue);
        end
        
        function covariance = get.Covar(this)
        %get.Covariance Getter for Covariance property
            covariance = this.pd_.var();
        end
    end
    
end

