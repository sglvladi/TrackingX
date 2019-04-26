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
            
            
            [model, other] = MeasurementX.extract_model(varargin{:});
            [tag, other] = MeasurementX.extract_tag(other{:});
            this@StateX(other{:});
            
            this.Model = model;
            this.Tag = tag;
        end
        
    end
    
    methods (Static)
        function [model, other] = extract_model(varargin)
            model = [];
            other = varargin;
            for i = 1:nargin
                if isa(varargin{i}, 'StateSpaceModelX')
                     model = varargin{i};
                     varargin{i} = [];
                     other = varargin;
                end
            end
        end
        
        function [tag, other] = extract_tag(varargin)
            tag = [];
            other = varargin;
            for i = 1:numel(varargin)
                if isa(varargin{i}, 'TagX')
                     tag = varargin{i};
                     varargin{i} = [];
                     other = varargin(~cellfun('isempty',varargin));
                end
            end
        end
    end
end

