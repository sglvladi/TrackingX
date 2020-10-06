classdef GOSPAX < MetricGeneratorX
% GOSPAX Class 
%
% Summary of GOSPAX:
% This is a class implementation of the Generalized Optimal Sub Pattern 
% Assignment (GOSPA) metric generator.
%
% GOSPAX Properties:
%   + CutOffThreshold
%   + Order
%
% GOSPAX Methods:
%    GOSPAX  - Constructor method
%    evaluate - Perform multinomial resampling 
% 
% See also SystematicResamplerX
    properties
        CutOffThreshold = 1;
        Order = 1;
        Alpha = 2;
        Mapping = [];
    end
    
    methods (Access=protected)
        function initialise_(this, config)
            if (isfield(config,'CutOffThreshold'))
                this.CutOffThreshold  = config.CutOffThreshold;
            end
            if (isfield(config,'Order'))
                this.Order  = config.Order;
            end
            if (isfield(config,'Alpha'))
                this.Alpha  = config.Alpha;
            end
        end
    end
    
    methods
        function this = GOSPAX(varargin)
        % GOSPAX Constructor method
        %   
        % DESCRIPTION: 
        % * GOSPAX(C,P,A) returns a "GOSPAX" object handle configured with 
        %   the cut-off threshold C, order paramater P and cardinality 
        %   penalty factor A.
        %
        % See also GOSPAX/evaluate
            
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
        
        function [metric,locDist,missDist,falseDist] = evaluate(this,GroundTruth,TrackList)
        % EVALUATE Evaluate the OSPA metric of a set of tracks compared to
        % the groundtrouth.
        %   
        % DESCRIPTION: 
        % * [METRIC,LOCDIST,MISSDIST,FALSEDIST] = EVALUATE(THIS, TRACKLIST, GROUNDTRUTH) 
        %   evaluates the GOSPA metric METRIC. including its localisation LOCDIST
        %   and cardinality CARDDIST components, of the TRACKLIST compared to the 
        %   GROUNDTRUTH.
        %
        % See also SystematicResamplerX/resample
            
            gndVectors = [GroundTruth.Vector];
            trackVectors = zeros(size(gndVectors,1),0);
            numTracks = numel(TrackList);
            for i = 1:numTracks
                trackVectors(:,i) = TrackList{i}.State.Vector;
            end
            
            [metric,~,decomp] = GOSPA(gndVectors,trackVectors,...
                                      this.Order,...
                                      this.CutOffThreshold,...
                                      this.Alpha);
            locDist = decomp.localisation;
            missDist = decomp.missed;
            falseDist = decomp.false;
        end
    end
end

