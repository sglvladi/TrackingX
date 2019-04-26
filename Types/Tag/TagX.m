classdef TagX < BaseX
% TagX class
% 
% Summary of TagX:
%   Class implementation of a generic TagX type.

    properties (Dependent)
        % ID: The id stored by the intentifier
        ID
    end       

    properties (Access = protected)
        % ID_: Internal storage for the id
        ID_ = []
    end
    
    methods
        function this = TagX(id)
        % TagX Base constructor
        %
        % Parameters
        % ----------
        % id: any
        %   The state vector
            
            this.ID_ = id;
        end
        
        function id = get.ID(this)
        %get.ID Custom getter for Vector property
            id = this.ID_;
        end
        
        function equals = eq(obj1, obj2)
        % Overrides equality operator
            equals = obj1.ID == obj2.ID;
        end
    end
    
end