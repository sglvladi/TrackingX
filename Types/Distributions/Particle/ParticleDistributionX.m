classdef ParticleDistributionX < ProbabilityDistributionX
% ParticleDistributionX class
%
% Summary of ParticleDistributionX:
% This is a class implementation of a Multivariate Gaussian distribution
%
% ParticleDistributionX Properties:
%   + NumParticles - The number of particles used to represent the distribution
%   + Particles - Particle States
%   + Weights - Particle Weights
%   + Mean - Mean vector of the Gaussian distribution
%   + Covar - The covariance matrix of the Gaussian distribution
%
% ParticleDistributionX Methods:
%   + ParticleDistributionX  - Constructor method
%   + random - Draw random samples from the multivariate normal distribution
%   + pdf - Evaluate the density of the multivariate normal distribution
%   + reset - Reset the distribution with a new number of random variables
%
% (+) denotes puplic properties/methods
% 
% See also TransitionMatrixX, MeasurementModelX and ControlModelX template classes

    properties (Dependent)
        % Particles: (NumVariables x NumParticles) column vector 
        %   Set of particles used to approximate the distribution
        Particles
        
        % Weights: (1 x NumParticles) matrix
        %   Particle importance weights
        Weights
        
        % NumParticles: scalar
        %   The number of particles used to represent the distribution
        NumParticles
        
        % Mean: (NumVariables x 1) column vector 
        %   Mean vector of the Particle distribution
        Mean
        
        % Covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Particle distribution
        Covar
    end
    
    properties (Access = protected)
        Particles_
        Weights_
        NumParticles_ = 1000;
        resampler_ = SystematicResamplerX();
    end
    
    methods
        function this = ParticleDistributionX(varargin)
        % ParticleDistributionX Construct a particle distribution
        %
        % Parameters
        % ----------
        % distribution: ProbabilityDistributionX
        %   Another distribution, based on which the new particle
        %   distribution should be created.
        % numParticles: scalar
        %   The number of particles to use to represent the new
        %   distribution
        %
        % -- OR ---
        %
        % mean: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % numParticles: scalar
        %   The number of particles to use to represent the new
        %   distribution
        %
        % - OR -
        %
        % particles: (NumVariables x NumParticles) matrix
        %   An initial set of particles
        % weights: (1 x NumParticles) matrix
        %   An initial set of weights
            this.reset(varargin{:});
        end
        
        function [samples, weights] = random(this, numSamples)
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
            
            weights_cpy = this.Weights;
            intensity = sum(weights_cpy);
            weights_cpy = weights_cpy./intensity;
            
            % Draw samples from the Particle distribution
            [samples, weights] = this.resampler_.resample(this.Particles,weights_cpy,numSamples);
            weights = weights.* intensity;
        end  
        
        function prob = pdf(this, samples)
        % pdf Evaluate the density of the multivariate normal distribution
        %   PROB = PDF(OBJ, SAMPLES) returns the density of the multivariate 
        %   normal distribution, evaluated at each column of SAMPLES.
        %
            prob = mvnpdf(samples', this.Mean', this.Covar)';
        end
        
        function [newParticles, newWeights, idx] = resample(this, numSamples)
        % resample Resample the underlying particles/weights. Sum of the
        %   weights is equal before and after resampling. 
        %
        % Parameters
        % ----------
        % numSamples: scalar, optional
        %   If specified, the distribution will resample itself to the provided
        %   number of samples numSamples.
        %
        % Returns
        % -------
        % newParticles: (NumStateDims x numSamples) double
        %   The newly resampled set of particles.
        % newWeights: (1 x numSamples) double
        %   The newly resampled set of weights.
        % idx: (1 x numSamples) int
        %   Mapping between resampled and old particles.
            
            % Ensure weights in resampling sum to 1
            weights_cpy = this.Weights;
            intensity = sum(weights_cpy);
            weights_cpy = weights_cpy./intensity;
            
            [newParticles, newWeights, idx] = this.resampler_.resample(this.Particles, weights_cpy, numSamples);
            newWeights = newWeights * intensity;
            this.reset(newParticles, newWeights);
        end
        
        function reset(this, varargin)
        % reset Reset the distribution with a new number of random variables
        %
        % Parameters
        % ----------
        % distribution: ProbabilityDistributionX
        %   Another distribution, based on which the new particle
        %   distribution should be created.
        % numParticles: scalar
        %   The number of particles to use to represent the new
        %   distribution
        %
        % -- OR ---
        %
        % mean: (NumVariables x 1) column vector
        %   Mean vector of the Gaussian distribution
        % covar: (NumVariables x NumVariables) matrix
        %   The covariance matrix of the Gaussian distribution
        % numParticles: scalar
        %   The number of particles to use to represent the new
        %   distribution
        %
        % - OR -
        %
        % particles: (NumVariables x NumParticles) matrix
        %   An initial set of particles
        % weights: (1 x NumParticles) matrix
        %   An initial set of weights
            
            switch(nargin)
                case(2)
                    if(isnumeric(varargin{1}))
                        this.NumParticles_ = varargin{1};
                        this.NumVariables = 1;
                    elseif(isa(varargin{1}, 'ParticleDistributionX'))
                        otherDist = varargin{1};
                        this.Particles = otherDist.Particles;
                        this.Weights = otherDist.Weights;
                    elseif(isa(varargin{1}, 'ProbabilityDistributionX'))
                        % Other distribution
                        otherDist = varargin{1};
                        numParticles = this.NumParticles_;
                        this.Particles = otherDist.random(numParticles);
                        this.Weights = repmat(1/numParticles,1,numParticles);
                    end
                otherwise
                    if(isa(varargin{1}, 'ParticleDistributionX'))
                        otherDist = varargin{1};
                        this.NumParticles_ = varargin{2};
                        [this.Particles, this.Weights] = otherDist.resample(this.NumParticles_);
                    elseif(isa(varargin{1}, 'ProbabilityDistributionX'))
                        % Other distribution
                        otherDist = varargin{1};
                        this.NumParticles_ = varargin{2};
                        this.Particles = otherDist.random(this.NumParticles_);
                        this.Weights = repmat(1/this.NumParticles_,1,this.NumParticles_);
                    else
                        % Particles
                        this.Particles = varargin{1};
                        this.Weights = varargin{2};
                    end
            end            
        end
    end
    
    methods
        function numParticles = get.NumParticles(this)
            if(isempty(this.NumParticles_))
                this.NumParticles_ = size(this.Particles,2);
            end
            numParticles = this.NumParticles_;
        end
        
        function particles = get.Particles(this)
            particles = this.Particles_;
        end
        
        function set.Particles(this, particles)
            this.Particles_ = particles;
            [this.NumVariables_, this.NumParticles_] = size(particles);
        end
        
        function weights = get.Weights(this)
            weights = this.Weights_;
        end
        
        function set.Weights(this, weights)
            this.Weights_ = weights;
        end
        
        function meanValue = get.Mean(this)
        %get.Mean Getter for Mean property
            meanValue = this.Particles*this.Weights';
        end
        
        function covar = get.Covar(this)
        %get.Covariance Getter for Covariance property
            covar = weightedcov(this.Particles,this.Weights);
        end
    end
    
end

