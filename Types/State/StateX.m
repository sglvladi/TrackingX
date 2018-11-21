classdef StateX < BaseX & dynamicprops
% StateX class
%
% Summary of StateX:
% Class implementation of the primitive State type
%
% StateX Properties:
%   + Distribution - A ProbabilityDistributionX object
% 
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties (Access = protected)
        Distribution_
    end
    
    properties (Dependent)
        Distribution
        Mean
        Covar
        NumStateDims
    end
    
    methods
        function mean = get.Mean(this)
        %get.Mean - Getter for Mean property
            mean = this.Distribution_.Mean;
        end
        function set.Mean(this, meanValue)
        %set.Mean - Setter for Mean property
            this.Distribution.Mean = double(meanValue);
        end
        function covar = get.Covar(this)
        %get.Covar Getter for Covar property
            covar = this.Distribution_.Covar;
        end
        function set.Covar(this, covariance)
        %set.Covariance Setter for Covariance property
            this.Distribution.Covar =  double(covariance);
        end       
        function dist = get.Distribution(this)
        %get.Distribution Getter for Distribution property
            dist = this.Distribution_;
        end
        function set.Distribution(this, dist)
        %set.Distribution Setter for Distribution property
            this.Distribution_ = this.setDistribution(dist);     
        end
    end
    
    methods
        function this = StateX(varargin)
        % StateX Constructor method
        %   
        % DESCRIPTION: 
        % * s = StateX() returns an unconfigured object handle.
        %
        %  See also addprop  
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                this.Distribution = varargin{1};
                return;
            end
        end
    end
    
    
    methods (Access = protected)
        function Dist = setDistribution(this,dist)
            Dist = dist;
        end
    end
end

