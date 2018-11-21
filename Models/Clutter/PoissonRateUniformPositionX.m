classdef PoissonRateUniformPositionX < ClutterModelX
% PoissonRateUniformPositionX class
%
% Summary of PoissonRateUniformPositionX
% This is a class implementation of a clutter model, described by a Poisson
% distributed false alarm rate and a Uniform spatial distribution.
%
% PoissonRateUniformPositionX Properties:
%   + ClutterRate - The clutter rate.
%   + Limits - Limits of the uniform spatial distribution of clutter.
%
% PoissonRateUniformPositionX Methods:
%   + random    - Birth intensity sample generator function, i.e. v_k ~ random(~)
%   + pdf       - Function to evaluate the instensity p(b_t) at a given point in
%                 the state space.
%
% (+) denotes puplic properties/methods
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool   
       
    properties (Dependent)
        ClutterRate
        Limits
    end
    
    properties %(Access=protected)
        poisson_ = PoissonDistributionX();
        uniform_ = UniformDistributionX();
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if(isfield(config,'ClutterRate'))
                this.ClutterRate = config.ClutterRate;
            end
            if(isfield(config,'Limits'))
                this.Limits = config.Limits;
            end
        end
    end
    
    methods
        function this = PoissonRateUniformPositionX(varargin)
        % PoissonRateUniformPositionX Constructor method
        % 
        % Parameters
        % ----------
        % ClutterRate: scalar
        %   The mean clutter rate of the underlying Poisson distribution.
        % Limits: (NumMeasDims x 2) matrix
        %   The limits of the underlying uniform spatial distribution.
        %
        % See also random, pdf.   
        
            % Call SuperClass method
            this@ClutterModelX();
                     
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
        
        function vk = random(this, varargin)
        % RANDOM Generates i.i.d. random samples from the clutter model intensity
        %
        % Parameters
        % ----------
        % Ns: scalar, optional
        %   Number of samples to generate.
        %   (default is 1 for 'cardinality' type, while for 'spatial' type it
        %    is determined by sampling from the underlying cardinality distribution)
        % type: string, optional
        %   Set to 'cardinality', to sample from the cardinality distribution, 
        %   or 'spatial' to sample from the spatial distribution.
        %   (default = 'spatial')
        %
        % Returns
        % -------
        % vk: (1 x Ns) vector, If type of 'cardinality' is used.
        %   A vector containing iid samples generated from the cardinality 
        %   distribution 
        % vk: (NumStateDims x Ns)  matrix, If type not set or set to 'spatial'. 
        %   A matrix, where each column corresponds to a sample drawn from 
        %   the clutter model intensity.
        %
        % See Also
        % --------
        % pdf
            
            switch(nargin)
                case(1)
                    type = 'spatial';
                    Ns = this.poisson_.random(1);
                case(2)
                    if(isnumeric(varargin{1}))
                        type = 'spatial';
                        Ns = varargin{1};
                    else
                        type = varargin{1};
                        if(strcmp(type,'spatial'))
                            Ns = this.poisson_.random();
                        else
                            Ns = 1;
                        end
                    end
                case(3)
                    type = varargin{2};
                    Ns = varargin{1};
            end
            
            if(strcmp(type,'spatial'))
                vk = this.uniform_.random(Ns);
            else
                vk = this.poisson_.random(Ns);
            end
        end
        
        function int = pdf(this, varargin)
        % pdf Evaluates the clutter intensity at given points xk of the state space
        %
        % Parameters
        % ----------
        % xk: matrix
        %   If type is 'spatial', xk should be a (NumStateDims x Ns) matrix, whose columns correspond to Ns
        %   individual spatial state vectors.
        %   If type is 'cardinality', xk should be a (1 x Ns) matrix, whose columns correspond to Ns
        %   individual cardinality state vectors.
        % type: string, optional
        %   Set to 'cardinality', to sample from the cardinality distribution, 
        %   or 'spatial' to sample from the spatial distribution.
        %   (default = 'spatial')
        %
        % Returns
        % -------
        % prob: matrix
        %   A (Nm x Ns) matrix, where :math:`prob(i,j)=p(y_k^i|x_k^j) 
        %
        % See also random
            switch(nargin)
                case(2)
                    type = 'spatial';
                    xk = varargin{1};
                otherwise
                    type = varargin{2};
                    xk = varargin{1};
            end
    
            if(strcmp(type,'spatial'))
                int = this.ClutterRate*this.uniform_.pdf(xk);
            else
                int = this.poisson_.pdf(xk);
            end
        end
        
        function clutterRate = get.ClutterRate(this)
            clutterRate = this.poisson_.Mean;
        end
        function set.ClutterRate(this, clutterRate)
            this.poisson_.Mean = clutterRate;
        end
        function limits = get.Limits(this)
            limits = this.uniform_.Limits;
        end
        function set.Limits(this, limits)
            this.uniform_.Limits = limits;
        end
    end
end