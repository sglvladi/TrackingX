classdef MeasurementX < StateX
% MeasurementX class
% 
% Summary of MeasurementX:
% Class implementation of the primitive Measurement type.
    
    properties
        % Model: StateSpaceModelX
        %   The ssm from which the measurement was generated
        Model = []
        
        % Tag: TagX
        %   An optional tag
        Tag = []
    end
    
    methods
        function this = MeasurementX(varargin)
        % Measurement Base constructor
        %
        % Parameters
        % ----------
        % Vector: numeric
        %   The state vector
        % Timestamp: datetime, optional
        %   An optional timestamp that can be attached to a StateX object.
        % Model: StateSpaceModelX, optional
        %   The ssm from which the measurement was generated
            
            
            [model, other] = MeasurementX.extract_model(varargin);
            [tag, other] = MeasurementX.extract_tag(varargin);
            this@StateX(other{:});
            
            this.Model = model;
            this.Tag = tag;
        end
        
    end
    
    methods (Static)
        function [model, other] = extract_model(varargs)
            model = [];
            other = varargs;
            for i = 1:nargin
                if isa(varargs{i}, 'StateSpaceModelX')
                     model = varargs{i};
                     varargs{i} = [];
                     other = varargs;
                end
            end
        end
        
        function [tag, other] = extract_tag(varargs)
            tag = [];
            other = varargs;
            for i = 1:numel(varargs)
                if isa(varargs{i}, 'TagX')
                     tag = varargs{i};
                     varargs{i} = [];
                     other = varargs(~cellfun('isempty',varargs));
                end
            end
        end
    end
end

