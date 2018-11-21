classdef GaussianDistributionX < ProbabilityDistributionX
% GaussianDistributionX class
%
% Summary of GaussianDistributionX:
% This is a class implementation of a Multivariate Gaussian distribution
%
% GaussianDistributionX Properties:
%   + Mean - Mean vector of the Gaussian distribution
%   + Covar - The covariance matrix of the Gaussian distribution
%
% GaussianDistributionX Methods:
%   + GaussianDistributionX  - Constructor method
%   + reset - Reset the distribution with a new number of random variables
%   + random - Draw random samples from the multivariate normal distribution
%   + pdf - Evaluate the density of the multivariate normal distribution
%   + fitToSamples - Fit a normal distribution to a set of weighted samples
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX template classes
    
    properties (Dependent)
        % Mean: (NumVariables x 1) column vector 
        %   Mean vector of the Gaussian distribution
        Mean
        
        % Covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        %   The covariance matrix has to be symmetric and positive semi-definite.
        Covar
    end
    
    properties (Access = private)
        % Mean_ - Internal storage for mean vector
        Mean_
        
        % Covar_ - Internal storage for covariance matrix
        Covar_
        
        % Std_ - Internal storage for standard deviation matrix
        Std_
    end
    
    methods
        function this = GaussianDistributionX(varargin)
        % MultivariateNormalDistribution Construct a multivariate Gaussian distribution
        %
        % Parameters
        % ----------
        % numVariables: scalar, optional
        %   The number of random variables. Onlt necessary if neither mean
        %   nor covar specified.
        % mean: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
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
            samples = mvnrnd(this.Mean',this.Covar,numSamples)';
        end  
        
        function prob = pdf(this, samples)
        % pdf Evaluate the density of the multivariate normal distribution
        %   PROB = PDF(OBJ, SAMPLES) returns the density of the multivariate 
        %   normal distribution, evaluated at each column of SAMPLES.
        %
            prob = mvnpdf(samples', this.Mean', this.Covar)';
        end
        
        function reset(this, varargin)
            
            if(nargin == 2)
                this.NumVariables = varargin{1};
            elseif(nargin>2)
                mean = varargin{1};
                covar =  varargin{2};
                %reset Reset the distribution with a new number of random variables

                % Reset number of random variables. This will also take care of
                % input validation.
                this.NumVariables = size(mean,1);

                % Re-initialize the mean and covariance            
                this.Mean = mean;
                this.Covar = covar;
            end
        end
    end
    
    methods (Static)
        function [sampleMean, sampleCovar] = fitToSamples(samples, weights)
        % fitToSamples Fit a normal distribution to a set of weighted samples
        %
        % Parameters
        % ----------
        % samples: (NumVariables x numSamples)
        %   A set of numSamples samples
        % weights: (1 x numSamples)
        %   Weights corresponding to each sample.
        %
        % Returns
        % -------
        % sampleMean: (NumVariables x 1)
        %   The estimated mean
        % sampleCovar: (NumVariables x NumVariables)
        %   The estimated covariance
            
            assert(isvector(weights));
            assert(size(samples,2) == size(weights,2));
            % Calculate weighted mean (assuming weights are normalized)
            %
            % sampleMean = sum(bsxfun(@times, samples, weights),1);
            sampleMean = samples*weights';
            
            % Calculating the covariance matrix is about 3 times
            % as expensive as calculating the mean, so skip it if output
            % is not needed.
            if nargout <= 1
                return;
            end
            
            % Calculate unbiased, weighted covariance matrix
            % (assuming weights are normalized)
            weightSqr = sum(weights.*weights);
            if abs(weightSqr - 1.0) < sqrt(eps)
                % To avoid NaN values, use factor of 1.0 if the sum of
                % squared weights is close to 1.0.
                % For example, this case can occur if there is only a
                % single particle in the input, or if all particles except
                % one have a weight of 0.
                % In this case, the covariance matrix will be all zeros,
                % which is the same behavior as the cov built-in.
                factor = 1.0;
            else
                factor = 1/(1 - weightSqr);
            end
            meanDiff = bsxfun(@minus, samples, sampleMean);
            
            sampleCovar = factor * (bsxfun(@times, meanDiff, weights) * meanDiff');
        end
    end
    
    methods
        function meanValue = get.Mean(this)
        %get.Mean Getter for Mean property
            meanValue = this.Mean_;
        end
        
        function set.Mean(this, meanValue)
        %set.Mean - Setter for Mean property
            this.Mean_ = double(meanValue);
            validateattributes(meanValue, {'numeric'}, {'size', [this.NumVariables,1]}, ...
                                           'GaussianDistributionX', 'Mean');
        end
        
        function covariance = get.Covar(this)
        %get.Covariance Getter for Covariance property
            covariance = this.Covar_;
        end
        
        function set.Covar(this, covariance)
        %set.Covariance Setter for Covariance property
            
            validateattributes(covariance, {'numeric'}, ...
                               {'size', [this.NumVariables, this.NumVariables]}, ...
                                'GaussianDistributionX', 'Covar');
            
            % Verify properties of covariance matrix. It's supposed to be
            % symmetric and positive semi-definite.
            % If that's the case, the number of negative eigenvalues,
            % numNegEigenValues, will be 0.
            % numNegEigenValues will be NaN if the input is not symmetric.
%             [T,numNegEigenValues] = cholcov(double(covariance));
%             if(numNegEigenValues~=0)
%                 warning('Applying correction to ensure symmetric and positive semi-definite covariance!');
%                 covariance = (covariance+covariance')/2;
%                 [T,numNegEigenValues] = cholcov(double(covariance));
%             end
%             assert(numNegEigenValues == 0);
            covariance = (covariance+covariance')/2;
            
            % Store Cholesky factor to speed up subsequent sampling
            % operations.
            this.Covar_ = double(covariance);
            %this.Std_ = real(T);
        end
    end
    
end

