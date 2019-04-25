classdef FieldBasedDeleterX < TrackInitiatorX 
% FieldBasedDeleterX class
%
% Summary of FieldBasedDeleterX:
% This is an implementation of a Track Deleter that deletes tracks based on
% the value of a given track field.
%
% FieldBasedDeleterX Properties:
%   ~ Fieldname - The name of the field whose value will be compared
%                 against the specified reference point 
%   ~ ReferenceValue - The value used as a deletion reference point
%   ~ ReferenceOperand - The binary operand used to compare the field and
%                        reference values
%
% FieldBasedDeleterX Methods:
%   + FieldBasedDeleterX - Constructor method
%
% (+) denotes public properties/methods
% (~) denotes protected properties/methods
%
% March 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
       Fieldname
       ReferenceValue
       ReferenceOperand = 'lt'
    end
        
    methods
        function this = FieldBasedDeleterX(varargin)
        % FieldBasedDeleterX Constructor method
        %
        % Parameters
        % ----------
        % Fieldname: string
        %   The name of the field whose value will be compared
        %   against the specified reference point.
        % ReferenceValue: any binary comparable type
        %   The value used as a deletion reference point. The type can be 
        %   any type that can be compared using binary arithmetic
        %   comparisons (i.e. <,>,<=,>=,==,~=)
        % ReferenceOperand: string, optional
        %   A string specifying the operand to be used when comparing the
        %   field (LHS) and reference point (RHS) values. The possible 
        %   values are:
        %       i) 'lt'  : Less than (<)
        %      ii) 'leq' : Less than or equal (<=)
        %     iii) 'gt'  : Greater than (>)
        %      iv) 'geq' : Greater than or equal (>=)
        %       v) 'eq'  : Equal (==)
        %      vi) 'neq' : Not equal (~=)
        %   The default is 'lt'.
        %
        % Usage
        % -----
        % * FieldBasedDeleterX(config) returns 
        %   a FieldBasedDeleterX object handle configured with the provided
        %   config.
            
            if(nargin==0)
                error('Not enough input arguments');
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.Fieldname = config.Fieldname;
                    this.ReferenceValue = config.ReferenceValue;
                    if(isfield(config,'ReferenceOperand'))
                        this.ReferenceOperand = config.ReferenceOperand;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.Fieldname = config.Fieldname;
            this.ReferenceValue = config.ReferenceValue;
            if(isfield(config,'ReferenceOperand'))
                this.ReferenceOperand = config.ReferenceOperand;
            end
        end
        
        function [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks)
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
        % * [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks)
        %   utilises the provided parameters to perform track deletion.
        %   Note that the provided DeleteThreshold parameter shall be used 
        %   to update the instance's respective property this.DeleteThreshold.
        % * [SurvivingTracks, DeletedTracks] = deleteTracks(this,Tracks)
        %   performs track deletion using the instance's this.DeleteThreshold. 
                    
            numTracks = numel(Tracks);
            DeletedTracks = [];
            for trackInd = 1:numTracks
                if isprop(Tracks{trackInd}, this.Fieldname)
                    value = get(Tracks{trackInd}, this.Fieldname);
                    if(this.evalOperand(value, this.ReferenceValue,...
                                        this.ReferenceOperand))
                        DeletedTracks{end+1} = Tracks{trackInd};
                        Tracks{trackInd} = [];
                    end 
                end
            end
            if(numel(Tracks))
                Tracks = Tracks(~cellfun('isempty',Tracks));
            end
            SurvivingTracks = Tracks;
        end
        
        function result = evalOperand(~, val, ref, operand)
            result = false;
            switch(operand)
                case 'lt'
                    result = val < ref;
                case 'gt'
                    result = val > ref;
                case 'leq'
                    result = val <= ref;
                case 'geq'
                    result = val >= ref;
                case 'eq'
                    result = val == ref;
                case 'neq'
                    result = val ~= ref;
            end
        end
    end
end