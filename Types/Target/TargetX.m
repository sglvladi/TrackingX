classdef TargetX < BaseX & dynamicprops
% TargetX class
%
% Summary of TargetX:
% This is a class implementation of a target object. TargetX objects can be 
% used to store the trajectories, as well as other relevant information, of 
% targets when developing/implementing MTT algorithms within TrackingX. 
%
% TargetX Properties:
%   - Trajectory = A (xDim x N) matrix representing the trajectory of a 
%     target 
%
% Notes:
% * TargetX is derived from the dynamicprops MATLAB class, meaning that
%   new class properties can be added to any TargetX instance as and when
%   needed. Type "help dynamicprops.addprop" for details on how to add new
%   properties.
% * Deleting a dynamic property can be achieved by calling the function
%   delete(findprop(t,'PropName')), where t is a TargetX object handle.
% 
% See also dynamicprops
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        Trajectory
    end
    
    methods
        function this = TargetX(varargin)
        % TARGETX Constructor method
        %   
        % DESCRIPTION: 
        % * t = TargetX() returns an unconfigured object handle.
        % * t = TargetX(traj) returns an object handle, who's "Trajectory"
        %   property contains the trajectory stored in the matrix traj.
        %   traj should be a (xDim x N) matrix, where xDim denotes the
        %   number of state dimensions and N denotes the length of the
        %   trajectory.
        % * t = TargetX(___,Name,Value,___) returns an object handle, with
        %   properties created for each Name,Value pair arguments. Since the 
        %   "Trajectory" property already exists by default, any Value
        %   passed to the constructor will be stored within the property.
        % * t = TargetX(config) returns an object handle, with
        %   properties created for each field of the config structure,
        %   using the field name as the name of the property and the field
        %   value as the value to be stored within that property. Since the 
        %   "Trajectory" property already exists by default, any Value
        %   passed to the constructor under "config.Trajectory" will be 
        %   stored within the property.
        %
        %  See also addprop  
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    fields = fieldnames(config);
                    for i = 1:numel(fields)
                        fieldname = fields{i};
                        if(strcmp(fieldname,"Trajectory"))
                            this.Trajectory = config.(fieldname);
                        continue;
                        end
                        p = this.addprop(fieldname);
                        p.NonCopyable = false;
                        this.(fieldname) = config.(fieldname);
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            fields = fieldnames(config);
            for i = 1:numel(fields)
                fieldname = fields{i};
                if(strcmp(fieldname,"Trajectory"))
                    this.Trajectory = config.(fieldname);
                continue;
                end
                p = this.addprop(fieldname);
                p.NonCopyable = false;
                this.(fieldname) = config.(fieldname);
            end
        end
    end
end

