classdef ConstantDetectionProbabilityX < DetectionModelX
% ConstantDetectionProbabilityX class
%
% Summary of ConstantDetectionProbabilityX
% This is a class implementation of a detection model, described by a Poisson
% distributed false alarm rate and a Uniform spatial distribution.
%
% ConstantDetectionProbabilityX Properties:
%   + DetectionProbability - The detection probability.
%
% ConstantDetectionProbabilityX Methods:
%   + pdf       - Function to evaluate the instensity p(b_t) at a given point in
%                 the state space.
%
% (+) denotes puplic properties/methods
%
%  February 2018 Lyudmil Vladimirov, University of Liverpool   
       
    properties
        DetectionProbability
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if(isfield(config,'DetectionProbability'))
                this.DetectionProbability = config.DetectionProbability;
            end
        end
    end
    
    methods
        function this = ConstantDetectionProbabilityX(varargin)
        % ConstantDetectionProbabilityX Constructor method
        % 
        % Parameters
        % ----------
        % DetectionProbability: scalar
        %   The detection probability over the search space.
        %
        % See also pdf.   
        
            % Call SuperClass method
            this@DetectionModelX();
                     
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
        
        function int = pdf(this, varargin)
        % pdf Evaluates the detection intensity at given points xk of the state space
        %
        % Parameters
        % ----------
        % xk: matrix
        %   A (NumStateDims x Ns) matrix, whose columns correspond to Ns individual state vectors.
        %
        % Returns
        % -------
        % prob: matrix
        %   A (1 x Ns) matrix, where :math:`prob(i,j)=p(D_k^i|x_k^j) 
            
            Ns = 1;
            if(nargin>1)
                xk = varargin{1};
                Ns = size(xk,2);
            end
            int = this.DetectionProbability*ones(1,Ns);
            x = [-3422.85261310741,-3191.58441229017,-2993.55608122595,-3228.08558459284, -3422.85261310741];
            y = [1.472966334316646e+03,1.123217741935625e+03,1.243952901078573e+03,1.579667558371468e+03, 1.472966334316646e+03];
            a = find(inpolygon(xk(1,:),xk(3,:),x,y));
            if numel(a)>0
                int(a) = 0.1;
            end
        end
    end
end