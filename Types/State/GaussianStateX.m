classdef GaussianStateX < StateX
% GaussianStateX class
%
% Summary of GaussianStateX:
% Class implementation of the primitive Gaussian State type. GaussianStateX
% is essentially a StateX that comes pre-initialised with a Gaussian
% Distribution.
%
% 
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
       
    methods
        function this = GaussianStateX(varargin)
        % GaussianStateX Constructor method
        %   
        % DESCRIPTION: 
        % * s = GaussianStateX() returns an unconfigured object handle.
        %
        %  See also addprop  
            
            switch(nargin)
                case(0)
                    % No arguments
                    this.Distribution = GaussianDistributionX();
                case(1)
                    % Single-argument/Distribution
                    this.Distribution = varargin{1};
                case(2)
                    % Two arguments: 1) mean, 2) covar
                    mean = varargin{1};
                    covar = varargin{2};
                    this.Distribution = GaussianDistributionX(mean,covar);
            end
        end
    end
    
    methods (Access = protected)
        function Dist = setDistribution(this,dist)
            assert(isa(dist,'GaussianDistributionX'));
            Dist = dist;
        end
    end
end

