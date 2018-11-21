classdef DistributionBasedBirthModelX < BirthModelX
% DistributionBasedBirthModelX class
%
% Summary of DistributionBasedBirthModelX
% This is a class implementation of a generic distribution-based birth model.
%
% DistributionBasedBirthModelX Properties:
%   + Distribution - The underlying model intensity distribution.
%
% DistributionBasedBirthModelX Methods:
%   + random    - Birth intensity sample generator function, i.e. v_k ~ random(~)
%   + pdf       - Function to evaluate the instensity p(b_t) at a given point in
%                 the state space.
%
% (+) denotes puplic properties/methods
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool   
       
    properties
        BirthIntensity = 0;
        Distribution
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if(isfield(config,'Distribution'))
                this.Distribution = config.Distribution;
            end
            if(isfield(config,'BirthIntensity'))
                this.BirthIntensity = config.BirthIntensity;
            end
        end
    end
    
    methods
        function this = DistributionBasedBirthModelX(varargin)
        % DistributionBasedBirthModelX Constructor method
        % 
        % Parameters
        % ----------
        % Distribution: ProbabilityDistributionX object handle
        %   The underlying intensity distribution of the birth model.
        % BirthIntensity: scalar
        %   The birth probability of the model.
        %
        % See also feval, random, pdf, covariance.   
        
            % Call SuperClass method
            this@BirthModelX();
                     
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
                             
        end
        
        function [vk,wk] = random(this, Ns)
        % RANDOM Generates i.i.d. random samples from the birth model intensity
        %
        % Parameters
        % ----------
        % Ns: scalar, optional
        %   Number of samples to generate 
        %   (default = 1)
        %
        % Returns
        % -------
        % vk: matrix
        %   A (NumStateDims x Ns) matrix, where each column corresponds to a
        %   sample drawn from the birth model intensity.
        %
        % See Also
        % --------
        % pdf
            
            if nargin<2
                Ns = 1;
            end
        
            % Generate noise samples
            vk = this.Distribution.random(Ns);
            
            if(nargout==2)
                wk = this.BirthIntensity(ones(1,Ns))/Ns;
            end
        end
        
        function int = pdf(this, xk)
        % pdf Evaluates the birth intensity p(b_t) at given points xk of the
        %   state space
        %
        % Parameters
        % ----------
        % xk: matrix
        %   A (NumStateDims x Ns) matrix, whose columns correspond to Ns
        %   individual state vectors
        %
        % Returns
        % -------
        % prob: matrix
        %   A (Nm x Ns) matrix, where :math:`prob(i,j)=p(y_k^i|x_k^j) 
        %
        % See also random
            
            int = this.BirthIntensity*this.Distribution.pdf(xk);
        end
    end
end