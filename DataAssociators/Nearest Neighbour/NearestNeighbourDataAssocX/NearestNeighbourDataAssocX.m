classdef NearestNeighbourDataAssocX < ProbabilisticDataAssocX
% NearestNeighbourDataAssocX class
%
% Summary of NearestNeighbourDataAssocX:
% This is a class implementation of a Joint Probabilistic Data Association Filter.
%
% NearestNeighbourDataAssocX Properties:
%   + TrackList - A (1 x Tracks) vector of TrackX objects
%   + MeasurementList - A MeasurementListX object, containing a (1 xTracks) 
%       array of MeasurementX objects
%   + Gater - A GaterX object used to perform gating
%   + Clusterer - A ClustererX object used to perform clustering of tracks
%   + ClutterModel - A ClutterModelX subclass used to model clutter
%
% NearestNeighbourDataAssocX Methods:
%   + NearestNeighbourDataAssocX  - Constructor method
%   + associate - Performs JPDAF association step
%   + updateTracks - Performs JPDAF update step
%
% (+) denotes puplic properties/methods
% (¬) denotes dependent properties
%
% See also ParticleFilterX, KalmanFilerX.
    
    methods
        function this = NearestNeighbourDataAssocX(varargin)
        % NearestNeighbourDataAssocX - Constructor method
        % 
        % Parameters
        % ----------
        % Gater: GaterX subclass, optional
        %   A gater object which should be used to perform gating of measurements. 
        %   (default = None, meaning that no gating is performed)
        % Clusterer: ClustererX subclass, optional
        %   A clusterer object which should be used to perform clustering of tracks. 
        %   (default = None, meaning that no clustering is performed)
        % Hypothesiser: HypothesiserX subclass, optional
        %   A hypothesiser object used generate and evaluate association
        %   hypotheses. (default = EfficientHypothesisManagementX());
        % DetectionProbability: scalar, optional
        %   The target detection probability (default = 1)
        % ClutterModel: ClutterModelX subclass
        %   A clutter model
        %
        % Usage
        % -----
        % * nn = NearestNeighbourDataAssocX(___,Name,Value) instantiates an  
        %   object handle, configured with the options specified by one or 
        %   more Name,Value pair arguments.
        %
        %  See also NearestNeighbourDataAssocX/associate, NearestNeighbourDataAssocX/updateTracks.   
                    
            this.Hypothesiser = AuctionHypothesiserX();
            
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
    end
end