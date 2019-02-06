classdef MeasurementListX < BaseX
% MeasurementListX class
% 
% Summary of MeasurementListX:
% Class implementation of the primitive Measurement List type.
    
    properties (Dependent)
        Measurements
        NumMeasurements
    end

    properties
        Vectors
        Timestamp
    end
    
    methods
        function this = MeasurementListX(varargin)
        % Measurement Base constructor
        %
        % Parameters
        % ----------
        % measurements: vector of MeasurementX objects
            
            if nargin>1
                this.Vectors = varargin{1};
                this.Timestamp = varargin{2};
            elseif nargin && ~isempty(varargin{1})
                if isa(varargin{1},'MeasurementListX')
                    this.Vectors = varargin{1}.Vectors;
                    this.Timestamp = varargin{1}.Timestamp;
                elseif isa(varargin{1},'MeasurementX')
                    measurements = varargin{1};
                    this.Vectors = [measurements.Vector];
                    this.Timestamp = measurements.Timestamp;
                else
                    this.Vectors = varargin{1};
                    this.Timestamp = [];
                end
            else
                this.Vectors = [];
                this.Timestamp = [];
            end
        end
        
        function numMeasurements = get.NumMeasurements(this)
            numMeasurements = size(this.Vectors,2);
        end
        
        function measurements = get.Measurements(this)
            if(isempty(this.Vectors))
                measurements = MeasurementX.empty();
            else
                for k = 1 : this.NumMeasurements
                    measurements(k) = MeasurementX(this.Vectors(:,k),...
                                                   this.Timestamp);
                end
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

