classdef TrackX < BaseX & dynamicprops
% TrackX class
%
% Summary of TrackX:
% This is a class implementation of a Track object. TrackX objects can be 
% used to store the trajectories, as well as other relevant information  
% about targets. 
%
% TrackX Properties:
%  + Trajectory - A (1 x N) StateX array representing the trajectory of a 
%    target.
%  + Tag - A TagX identifier object
%
% Notes:
% * TrackX is derived from the dynamicprops MATLAB class, meaning that
%   new class properties can be added to any TrackX instance as and when
%   needed. Type "help dynamicprops.addprop" for details on how to add new
%   properties.
% * Deleting a dynamic property can be achieved by calling the function
%   delete(findprop(t,'PropName')), where t is a TrackX object handle.
% 
% See also dynamicprops
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties (Dependent)
        State
        TimeOfLastUpdate
        TimeOfInitiation
    end

    properties
        Tag
        Trajectory
    end
    
    methods
        function this = TrackX(varargin)
        % TrackX Constructor method
        %
        % Parameters
        % ----------
        % Trajectory: (1 x N) array of StateX objects
        %   A sequence of (timestamped) State objects
        % Tag: TagX
        %   A (unique) track identifier
        %
        % Description
        % -----------
        %   TrackX() returns an unconfigured/placeholder TrackX object.
        %   
        %   TrackX(T) returns a TrackX object, where T can either be a StateX
        %   array, i.e. a Trajectory, or a TagX object.
        %
        %   TrackX(T,t) returns a TrackX object where T is a trajectory 
        %   and t is a tag/identifier.
        %
        %   TrackX(Name,Value,...) where parameters are supplied in the form
        %   of (Name,Value) argument pairs. If the constructor comes across 
        %   any arguments whose Name do not match an existing class property, 
        %   then new (public) properties will be created for the new object.  
        %
        %   TrackX(S) where S is a structure with appropriately set fields.
        %   Any fieldnames that cannot be matched to existing class properties
        %   will result in new properties being declared for the object.
        %
        %  See also addprop  
            
            switch(nargin)
                case 0
                    return;
                case 1
                    switch(class(varargin{1}))
                    % 
                        case 'cell'
                            fields = varargin{1};
                            for field = fields
                                if ~isprop(this,field{:})
                                    p = this.addprop(field{:});
                                    p.NonCopyable = false;
                                end
                            end
                        case 'struct'
                            config = varargin{1};
                            fields = fieldnames(config);
                            for i = 1:numel(fields)
                                fieldname = fields{i};
                                switch(fieldname)
                                    case("Trajectory")
                                        this.Trajectory = config.(fieldname);
                                    case('Tag')
                                        this.Tag = config.(fieldname);
                                    otherwise
                                        if ~isprop(this,fieldname)
                                            p = this.addprop(fieldname);
                                            p.NonCopyable = false;
                                            this.(fieldname) = config.(fieldname);
                                        end
                                end
                            end
                        otherwise
                            if isa(varargin{1},'StateX')
                                this.Trajectory = varargin{1};
                            elseif isa(varargin{1},'TagX')
                                this.Tag = varargin{1};
                            end
                            return;
                    end
                    return;
                case 2
                    if isa(varargin{1},'StateX')
                        this.Trajectory = varargin{1};
                        this.Tag = varargin{2};
                        return;
                    end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            fields = fieldnames(config);
            for i = 1:numel(fields)
                fieldname = fields{i};
                switch(fieldname)
                    case("Trajectory")
                        this.Trajectory = config.(fieldname);
                    case('Tag')
                        this.Tag = config.(fieldname);
                    otherwise
                        if ~isprop(this,fieldname)
                            p = this.addprop(fieldname);
                            p.NonCopyable = false;
                            this.(fieldname) = config.(fieldname);
                        end
                end
            end
        end
        
        
        function state = get.State(this)
        % get.State - Getter for dependent State property
            state = this.Trajectory(end);
        end
        
        function toi = get.TimeOfInitiation(this)
        % get.TimeOfInitiation - Getter for dependent TimeOfInitiation property
            toi = this.Trajectory(1).Timestamp;
        end
        
        function lupd = get.TimeOfLastUpdate(this)
        % get.TimeOfLastUpdate - Getter for dependent TimeOfLastUpdate property
            lupd = this.Trajectory(end).Timestamp;
        end
    end
end

