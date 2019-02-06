classdef TrackListX < BaseX
% TrackListX class
% 
% Summary of TrackListX:
% Class implementation of the primitive Track List type.
    
    properties (Dependent)
        Trajectories
        Tags
        Timestamps
        NumTracks
    end

    properties
        Tracks
    end
    
    methods
        function this = TrackListX(varargin)
        % TrackListX Base constructor
        %
        % Parameters
        % ----------
        % tracks: cell or TrackX
        %   A cell or object array of TrackX objects
        %   
            if nargin && ~isempty(varargin{1})
                if isa(varargin{1},'TrackListX')
                    tracks = varargin{1}.Tracks;
                    this.Tracks = tracks.copy();
                else
                    tracks = varargin{1};
                    numTracks = numel(tracks);
                    if iscell(tracks)
                        for i=1:numTracks
                            dTracks(i) = tracks{i};
                        end
                    else
                        dTracks = tracks;
                    end
                    this.Tracks = dTracks;
                end
            else 
                this.Tracks = [];
            end
        end
        
        function Vectors = get.Vectors(this)
            if(isempty(this.Tracks))
                Vectors = [];
            else
                Vectors = [this.Tracks.Vector];
            end
        end
        
        function Timestamps = get.Timestamps(this)
            if(isempty(this.Tracks))
                Timestamps = [];
            else
                Timestamps = [this.Tracks.Timestamp];
            end
        end
        
        function numTracks = get.NumTracks(this)
            if(isempty(this.Tracks))
                numTracks = 0;
            else
                numTracks = numel(this.Tracks);
            end
        end
    end
    
    methods (Static)
        function [model, other] = extract_model(varargs)
            model = [];
            other = varargs;
            if numel(varargs) && isa(varargs{end}, 'StateSpaceModelX')
                model = varargs{end};
                other = varargs(1:end-1);
            end
        end
    end
end

