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
        Models
        Tags
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
                    this.Tags = varargin{1}.Tags;
                    this.Models = varargin{1}.Models;
                    this.Timestamp = varargin{1}.Timestamp;
%                     if(isfield(varargin{1},'Tags'))
%                         this.Tags
                elseif isa(varargin{1},'MeasurementX')
                    measurements = varargin{1};
                    this.Vectors = [measurements.Vector];
                    this.Tags = [measurements.Tag];
                    this.Models = [measurements.Model];
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
                    if ~isempty(this.Models)
                        model = this.Models(k);
                    else
                        model = [];
                    end
                    if ~isempty(this.Tags)
                        tag = this.Tags(k);
                    else
                        tag = [];
                    end
                    measurements(k) = MeasurementX(this.Vectors(:,k),...
                                                   this.Timestamp,...
                                                   model,...
                                                   tag);
                end
            end
        end
    end
end

