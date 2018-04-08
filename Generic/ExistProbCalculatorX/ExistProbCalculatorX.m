classdef ExistProbCalculatorX < BaseX 
% ExistProbCalculatorX class
%
% Summary of ExistProbCalculatorX:
% This class can be used to compute (recursively) the existence probability
% of tracks.
%
% ExistProbCalculatorX Properties:
%   TrackList
%   AssocWeightsMatrix
%   ProbOfDetection
%   ProbOfDeath
%   ProbOfGating
%
% ExistProbCalculatorX Methods:
%   + ExistProbCalculatorX - Constructor method
%   + calculate - Calculate the existence probabilities
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       TrackList
       AssocWeightsMatrix
       ProbOfDetection
       ProbOfDeath
       ProbOfGating
    end
        
    methods
        function this = ExistProbCalculatorX(varargin)
        % EXISTPROBCALCULATORX Constructor method
        %
        % Parameters
        % ----------
        % ProbOfDetection: scalar
        %   A scalar value in the range [0,1] which specifies the probability
        %   of detection for all targets
        % ProbOfDeath: scalar
        %   A scalar value in the range [0,1] which specifies the probability
        %   of death for all targets
        % ProbOfGating: scalar
        %   A scalar value in the range [0,1] which specifies the probability
        %   of gating for all targets
        % 
        % Usage
        % -----
        % * ExistProbCalculatorX() returns a ExistProbCalculatorX object handle.
        % * ExistProbCalculatorX(ProbOfDetection,ProbOfDeath,ProbOfGating) returns 
        %   a ExistProbCalculatorX( object handle configured with the provided
        %   probability of detection ProbOfDetection, probability of death ProbOfDeath 
        %   and probability of gating ProbOfGating. This is particularly
        %   useful when all targets are characterised by the same parameters.
        % * ExistProbCalculatorX(config) can also be used, where config is a
        %   structure with fields config.ProbOfDetection, config.ProbOfDeath and
        %   config.ProbOfGating.
        % * PhdTrackInitiatorX(___,Name,Value) instantiates an object 
        %   handle, configured with additional options specified by one or
        %   more Name,Value pair arguments.
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.ProbOfDetection = config.ProbOfDetection;
                    this.ProbOfDeath = config.ProbOfDeath;
                    this.ProbOfGating = config.ProbOfGating;
                end
                return;
            elseif(nargin==3)
                this.ProbOfDetection = varargin{1};
                this.ProbOfDeath = varargin{2};
                this.ProbOfGating = varargin{3};
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.ProbOfDetection = config.ProbOfDetection;
            this.ProbOfDeath = config.ProbOfDeath;
            this.ProbOfGating = config.ProbOfGating;
        end
        
        function argout = calculate(this,varargin)
        % CALCULATE Calculates the probability(ies) of existence
        % 
        % Parameters
        % ----------
        % varargin{1}: (1 x NumTracks) cell array or TrackX object
        %   Either, a cell array of TrackX objects or a single TrackX
        %   object.
        % varargin{2}: (NumTracks x Nm+1) or (1 x Nm+1) matrix
        %   In the case that varargin{1} is a cell array, varargin{2}
        %   should be a (NumTracks x Nm+1) matrix, else it should be a 
        %   (1 x Nm+1) row vector, where in both cases, each column contains 
        %   the association weight between a measurement and a track. The
        %   first column corresponds to the dummy measurement.
        %
        % Returns
        % -------
        % argout: (1 x NumTracks) cell array or TrackX object
        %   Either, a cell array of TrackX objects or a single TrackX
        %   object, with their/its ProbOfExist field updated
        %
        % Usage
        % -----
        % * TrackList = calculate(this,TrackList,AssocWeightsMatrix) 
        %   computes the new probability of existence for each TrackX
        %   object in TrackList and returns the updated TrackList
        % * Track = calculate(this,Track,AssocWeightsMatrix) 
        %   computes the new probability of existence for the TrackX
        %   object Track and returns the updated Track
            if(isempty(varargin{1}))
                argout = varargin{1};
                return
            end
            if(iscell(varargin{1}))
                % Compute Existence Probabilities
                this.TrackList = varargin{1};
                this.AssocWeightsMatrix = varargin{2};
                numTracks = numel(this.TrackList);
                for trackInd = 1:numTracks
                    this.TrackList{trackInd}.ProbOfExist = (1 - this.ProbOfDeath)*this.TrackList{trackInd}.ProbOfExist;
                    denom = sum(this.AssocWeightsMatrix(trackInd,:))*this.TrackList{trackInd}.ProbOfExist + this.AssocWeightsMatrix(trackInd,1)*(1-this.TrackList{trackInd}.ProbOfExist)/((1-this.ProbOfDetection*this.ProbOfGating));
                    this.TrackList{trackInd}.ProbOfExist = (sum(this.AssocWeightsMatrix(trackInd,:))*this.TrackList{trackInd}.ProbOfExist)/denom;
                end
                argout = this.TrackList;
            elseif(isa(varargin{1},'TrackX'))
                Track = varargin{1};
                AssocWeights = varargin{2};
                Track.ProbOfExist = (1 - this.ProbOfDeath)*Track.ProbOfExist;
                denom = sum(AssocWeights)*Track.ProbOfExist + AssocWeights(1)*(1-Track.ProbOfExist)/((1-this.ProbOfDetection*this.ProbOfGating));
                Track.ProbOfExist = (sum(AssocWeights)*Track.ProbOfExist)/denom;
                argout = Track;
            end
        end 
    end
end