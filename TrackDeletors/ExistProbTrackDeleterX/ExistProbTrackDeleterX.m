classdef ExistProbTrackDeleterX < TrackInitiatorX 
% ExistProbTrackDeleterX Abstract class
%
% Summary of ExistProbTrackDeleterX:
% This is an implementation of a Track Initiator which uses a PHD Filter
% to detect and initiate tracks.
%
% ExistProbTrackDeleterX Properties:
%   DeleteThreshold
%
% ExistProbTrackDeleterX Methods:
%   + ExistProbTrackDeleterX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       DeleteThreshold
    end
        
    methods
        function this = ExistProbTrackDeleterX(varargin)
        % ExistProbTrackDeleterX Constructor method
        %
        % Parameters
        % ----------
        % DeleteThreshold: scalar
        %   A scalar value in the range [0,1] which specifies the threshold
        %   probability of existance to be used to delete tracks. 
        % 
        % Usage
        % -----
        % * ExistProbTrackDeleterX(DeleteThreshold) returns 
        %   a ExistProbTrackDeleterX object handle configured with the provided
        %   deletion threshold DeleteThreshold.
        % * PhdTrackInitiatorX(config) can also be used, where config is a
        %   structure with fields config.DeleteThreshold
        % * PhdTrackInitiatorX(___,Name,Value) instantiates an object 
        %   handle, configured with additional options specified by one or
        %   more Name,Value pair arguments.
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.DeleteThreshold = config.DeleteThreshold;
                else
                    this.DeleteThreshold = varargin{1};
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.DeleteThreshold = config.DeleteThreshold;
        end
        
        function [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks,DeleteThreshold)
        % DELETETRACKS Perform track deletion
        % 
        % Parameters
        % ----------
        % Tracks: (1 x NumTracks) cell array
        %   A cell array of TrackX objects, from which tracks will be
        %   potentially deleted
        % DeleteThreshold: scalar, optional
        %   A scalar value in the range [0,1] which specifies the threshold
        %   probability of existance to be used to delete tracks.
        %
        % Returns
        % -------
        % SurvivingTracks: (1 x NumSurvivingTracks) cell array
        %   A cell array of TrackX objects, each corresponding to a TrackX
        %   object from TrackList, which has survived the deletion process.
        % DeletedTracks: (1 x NumDeletedTracks) cell array
        %   A cell array of TrackX objects, each corresponding to a TrackX
        %   object from TrackList, which has not survived the deletion process.
        %
        % Usage
        % -----
        % * [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks,DeleteThreshold)
        %   utilises the provided parameters to perform track deletion.
        %   Note that the provided DeleteThreshold parameter shall be used 
        %   to update the instance's respective property this.DeleteThreshold.
        % * [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks)
        %   performs track deletion using the instance's this.DeleteThreshold. 
        
            if(nargin==3)
               this.DeleteThreshold = DeleteThreshold;
            end
            
            numTracks = numel(Tracks);
            DeletedTracks = [];
            for trackInd = 1:numTracks
                if(Tracks{trackInd}.ProbOfExist<0.1)
                    DeletedTracks{end+1} = Tracks{trackInd};
                    Tracks{trackInd} = [];
                end
            end
            if(numel(Tracks))
                Tracks = Tracks(~cellfun('isempty',Tracks));
            end
            SurvivingTracks = Tracks;
        end 
    end
end