classdef RMSEX < MetricGeneratorX
% RMSEX Class 
%
% Summary of RMSEX:
% This is a class implementation of the Optimal Sub Pattern Assignment (OSPA)
% metric generator.
%
% RMSEX Properties:
%   + CutOffThreshold
%   + Order
%
% MultinomialResamplerX Methods:
%    RMSEX  - Constructor method
%    evaluate - Perform multinomial resampling 
% 
% See also SystematicResamplerX
    properties
        Mapping = [];
    end
    
    methods (Access=protected)
        function initialise_(this, config)
        end
    end
    
    methods
        function this = RMSEX(varargin)
        % RMSEX Constructor method
        %   
        % DESCRIPTION: 
        % * RMSEX(C,P) returns a "RMSEX" object handle configured with the
        %   cut-off threshold C and order paramater P.
        %
        % See also RMSEX/evaluate
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function metric = evaluate(this,GroundTruth,Track)
        % EVALUATE Evaluate the OSPA metric of a set of tracks compared to
        % the groundtrouth.
        %   
        % DESCRIPTION: 
        % * [METRIC,LOCDIST,CARDDIST] = EVALUATE(THIS, TRACKLIST, GROUNDTRUTH) 
        %   evaluates the OSPA metric METRIC, including its localisation LOCDIST
        %   and cardinality CARDDIST components, of the TRACKLIST compared to the 
        %   GROUNDTRUTH.
        %
        % See also SystematicResamplerX/resample
            
            gndVectors = [GroundTruth.Trajectory.Vector];
            trackVectors = [Track.Trajectory.Vector];
            n = size(gndVectors,2);
            mse = zeros(1,n);
            for i=1:n
                mse(i) = ((gndVectors(1,i) - trackVectors(1,i))^2 + (gndVectors(3,i) - trackVectors(3,i))^2);
            end
            metric = mean(sqrt(sum(mse)));
        end
    end
end

