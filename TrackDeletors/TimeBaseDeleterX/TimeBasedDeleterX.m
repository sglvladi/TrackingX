classdef TimeBasedDeleterX < TrackInitiatorX 
% TimeBasedDeleterX Abstract class
%
% Summary of TimeBasedDeleterX:
% This is an implementation of a Track Deleter that deletes tracks based on
% a threshold set on the time since a track was last updated.
%
% TimeBasedDeleterX Properties:
%   DeleteThreshold
%
% TimeBasedDeleterX Methods:
%   + TimeBasedDeleterX - Constructor method
%
% (+) denotes public properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       DeleteThreshold
    end
        
    methods
        function this = TimeBasedDeleterX(varargin)
        % TimeBasedDeleterX Constructor method
        %
        % Parameters
        % ----------
        % DeleteThreshold: timedelta
        %   Specifies the threshold time difference to be used to delete tracks. 
        % 
        % Usage
        % -----
        % * TimeBasedDeleterX(DeleteThreshold) returns 
        %   a TimeBasedDeleterX object handle configured with the provided
        %   deletion threshold DeleteThreshold.
            
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
        
        function [SurvivingTracks, DeletedTracks] = deleteTracks(this, Tracks, timestamp)
        % DELETETRACKS Perform track deletion
        % 
        % Parameters
        % ----------
        % Tracks: (1 x NumTracks) TrackX array
        %   An array of TrackX objects, from which tracks will be
        %   potentially deleted
        %
        % Returns
        % -------
        % SurvivingTracks: (1 x NumSurvivingTracks) TrackX array
        %   An array of TrackX objects, each corresponding to a TrackX
        %   object from TrackList, which has survived the deletion process.
        % DeletedTracks: (1 x NumDeletedTracks) TrackX array
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
            
%             numTracks = numel(Tracks);
%             SurvivingTracks = TrackX.empty();
%             DeletedTracks = TrackX.empty();
%             for trackInd = 1:numTracks
%                 if(timestamp - Tracks(trackInd).TimeOfLastUpdate<this.DeleteThreshold)
%                     DeletedTracks(end+1) = Tracks(trackInd);
%                 else
%                     SurvivingTracks(end+1) = Tracks(trackInd);
%                 end
%             end
            numTracks = numel(Tracks);
            DeletedTracks = [];
            for trackInd = 1:numTracks
                if(timestamp - Tracks{trackInd}.TimeOfLastUpdate>this.DeleteThreshold)
%                     disp('Deleting TIME based!');
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