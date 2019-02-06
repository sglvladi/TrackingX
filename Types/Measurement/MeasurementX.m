classdef MeasurementX < StateX
% MeasurementX class
% 
% Summary of MeasurementX:
% Class implementation of the primitive Measurement type.
    
    properties
        % Model: StateSpaceModelX
        %   The ssm from which the measurement was generated
        Model = []
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
            
            this@StateX(other{:});
            
            this.Model = model;
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

