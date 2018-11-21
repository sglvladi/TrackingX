classdef (Abstract) ProbabilityDistributionX < BaseX
% ProbabilityDistributionX class
%
% Summary of ProbabilityDistributionX:
% This is the base class for all distribution objects.
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
 
    properties (Dependent, SetAccess = protected)
        % NumVariables - Number of random variables for this distribution
        NumVariables
    end       

    properties (Access = protected)
        % NumVariables_ - Internal storage for number of random variables
        NumVariables_
    end

    methods (Abstract)
        % Random Draw random samples from the distribution
        %
        % Parameters
        % ----------
        % numSamples: scalar
        %   Number of samples to draw from distribution
        %
        samples = random(this, numSamples)

        % Reset Reset the distribution 
        % 
        % Parameters
        % ----------
        % numVariables: scalar
        %   Number of random variables that are described by the
        %   distribution
        reset(this, numVariables);            
    end

    methods
        function this = ProbabilityDistributionX()
        %ProbabilityDistribution Construct a numVariables-variate probability distribution
            
        end

        function numVariables = get.NumVariables(this)
        %get.NumVariables Custom getter for NumVariables property
            numVariables = this.NumVariables_;
        end

        function set.NumVariables(this, numVariables)
        %set.NumVariables Custom setter for NumVariables property
            validateattributes(numVariables, {'numeric'}, {'scalar', 'integer', '>=', 1}, ...
                'ProbabilityDistribution', 'numVariables');
            this.NumVariables_ = double(numVariables);
        end
    end
end
