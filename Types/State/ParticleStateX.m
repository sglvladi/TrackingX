classdef ParticleStateX < ParticleDistributionX & StateX
% ParticleStateX class
%
% Summary of ParticleStateX:
% Class implementation of the primitive Particle State type. ParticleStateX
% is essentially a StateX that comes pre-initialised with a Particle
% Distribution.
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    methods
        function this = ParticleStateX(varargin)
        % ParticleStateX Construct a particle state
        %
        % Parameters
        % ----------
        % distribution: ProbabilityDistributionX
        %   Another distribution, based on which the new particle
        %   distribution should be created.
        % numParticles: scalar
        %   The number of particles to use to represent the new
        %   distribution
        % timestamp: datetime, optional
        %   The state timestamp
        %
        % -- OR ---
        %
        % particles: (NumVariables x NumParticles) matrix
        %   An initial set of particles
        % weights: (1 x NumParticles) matrix
        %   An initial set of weights
        % timestamp: datetime, optional
        %   The state timestamp
            
            [timestamp, other] = StateX.extract_timestamp(varargin);
            
            this@StateX(timestamp);
            this@ParticleDistributionX(other{:});
        end
    end
    
    methods (Access=protected)
        function stateVector = getVector(this)
            stateVector = this.Mean;
        end
        function Vector = setVector(this, stateVector)
            error("Setting the Vector property is not currently supported for ParticleStateX objects!");
        end
        function numDims = getNumDims(this)
            numDims = this.NumVariables;
        end
    end
    
end

