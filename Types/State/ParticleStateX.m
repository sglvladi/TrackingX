classdef ParticleStateX < StateX
% ParticleStateX class
%
% Summary of ParticleStateX:
% Class implementation of the primitive Particle State type. ParticleStateX
% is essentially a StateX that comes pre-initialised with a Particle
% Distribution.
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties (Dependent)
        NumParticles
        Particles
        Weights
    end
    methods
        function this = ParticleStateX(varargin)
        % ParticleStateX Constructor method
        %   
        % DESCRIPTION: 
        % * s = ParticleStateX() returns an unconfigured object handle.
        %
        %  See also addprop  
            
            switch(nargin)
                case(0)
                    % No arguments
                    this.Distribution = ParticleDistributionX();
                case(1)
                    % Single-argument/Distribution
                    this.Distribution = varargin{1};
                case(2)
                    % Single-argument/Distribution
                    particles = varargin{1};
                    weights = varargin{2};
                    this.Distribution = ParticleDistributionX(particles,weights);
                case(3)
                    % Two arguments: 1) mean, 2) covar
                    mean = varargin{1};
                    covar = varargin{2};
                    numParticles = varargin{3};
                    this.Distribution = ParticleDistributionX(mean,covar,numParticles);
            end
        end
        
        function numParticles = get.NumParticles(this)
        %get.NumParticles Getter for NumParticles property
            numParticles = size(this.Distribution_.Particles,2);
        end
        function particles = get.Particles(this)
        %get.Particles Getter for Particles property
            particles = this.Distribution_.Particles;
        end
        function set.Particles(this, particles)
        %set.Particles Setter for Particles property
           this.Distribution_.Particles = particles;
        end
        function weights = get.Weights(this)
        %get.Particles Getter for Weights property
            weights = this.Distribution_.Weights;
        end
        function set.Weights(this, weights)
        %set.Particles Setter for Weights property
            this.Distribution_.Weights= weights;
        end
    end
    
    methods (Access = protected)
        function Dist = setDistribution(this,dist)
            assert(isa(dist,'ParticleDistributionX'));
            Dist = dist;
        end
    end
end

