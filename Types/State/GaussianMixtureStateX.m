classdef GaussianMixtureStateX < StateX
% GaussianMixtureStateX class
%
% Summary of GaussianMixtureStateX:
% Class implementation of the primitive Gaussian Mixture State type. 
% GaussianMixtureStateX is essentially a StateX that comes pre-initialised with a Gaussian
% Distribution.
%
% 
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
      
    properties (Dependent)
        NumComponents
        Components
        Means
        Covars
        Weights
    end
    
    methods
        function this = GaussianMixtureStateX(varargin)
        % GaussianMixtureStateX Constructor method
        %   
        % DESCRIPTION: 
        % * s = GaussianMixtureStateX() returns an unconfigured object handle.
        %
        %  See also addprop  
            
            switch(nargin)
                case(0)
                    % No arguments
                    this.Distribution = GaussianMixtureX();
                case(1)
                    % Single-argument/Distribution
                    if(isa(varargin{1}, 'ProbabilityDistributionX'))
                        this.Distribution = copy(varargin{1});
                    elseif (iscell(varargin{1}))
                        components = varargin{1};
                        this.Distribution = GaussianMixtureX(components);
                    end
                case(2)
                    % Two arguments: 1) mean, 2) covar
                    means = varargin{1};
                    covars = varargin{2};
                    this.Distribution = GaussianMixtureX(means,covars);
                case(3)
                    % Two arguments: 1) mean, 2) covar
                    means = varargin{1};
                    covars = varargin{2};
                    weights = varargin{3};
                    this.Distribution = GaussianMixtureX(means,covars,weights);
            end
        end
        
        function numComponents = get.NumComponents(this)
        %get.NumComponents Getter for NumComponents property
            numComponents = size(this.Distribution_.Components,2);
        end
        function components = get.Components(this)
        %get.Components Getter for Components property
            components = this.Distribution_.Components;
        end
        function set.Components(this, components)
        %set.Components Setter for Components property
           this.Distribution_.Components = components;
        end
        function means = get.Means(this)
        %get.Components Getter for Weights property
            means = this.Distribution_.Means;
        end
        function set.Means(this, means)
        %set.Components Setter for Weights property
            this.Distribution_.Means= means;
        end
        function covars = get.Covars(this)
        %get.Components Getter for Weights property
            covars = this.Distribution_.Covars;
        end
        function set.Covars(this, covars)
        %set.Components Setter for Weights property
            this.Distribution_.Covars= covars;
        end
        function weights = get.Weights(this)
        %get.Components Getter for Weights property
            weights = this.Distribution_.Weights;
        end
        function set.Weights(this, weights)
        %set.Components Setter for Weights property
            this.Distribution_.Weights= weights;
        end
    end
    
    methods (Access = protected)
        function Dist = setDistribution(this,dist)
            assert(isa(dist,'GaussianMixtureX'));
            Dist = dist;
        end
    end
end

