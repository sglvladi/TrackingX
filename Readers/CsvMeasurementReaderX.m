classdef CsvMeasurementReaderX < MeasurementReaderX
% CsvMeasurementReaderX Abstract class
%
% Summary of CsvMeasurementReaderX
%  This is the base class for all TrackingX measurement readers. Must be 
%  used as the superclass of any custom defined TrackingX measurement readers. 
%
% CsvMeasurementReaderX Properties:
%   + Measurements - The measurements at a given timestep.
%
% CsvMeasurementReaderX Methods: 
%   ~ read - Read and return measurements at current timestep
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility with the TrackingX framework.
    
    properties %(Access = protected)
        Table_
        RowIndex_ = 1
        TimeValues_
    end
    properties
        Filename
        StateFields
        TimeField
        TimeFormat
        TimeZone = 'UTC'
        MetadataFields
    end

    methods (Access = protected)
        function initialise_(this, config)
            this.Filename = config.Filename;
            this.Table_ = readtable(this.Filename,'Delimiter', ',');
            this.StateFields = config.StateFields;
            this.TimeField = config.TimeField;
            this.TimeValues_ = this.Table_{:,{this.TimeField}};
            if(isfield(config,'TimeFormat'))
                this.TimeFormat = config.TimeFormat;
            end
            if(isfield(config,'TimeZone'))
                this.TimeZone = config.TimeZone;
            end
            if(isfield(config,'MetadataFields'))
                this.MetadataFields = config.MetadataFields;
            end
        end
    end
    
    methods
        
        function this = CsvMeasurementReaderX(varargin)
        % CsvMeasurementReaderX Constructor method
        %
        % Parameters
        % ----------
        % Filename: string
        %   Path to csv file
        % StateFields: (1 x N) cell array
        %   List of columns names to be used in state vector
        % TimeField: string
        %   Name of column to be used as time field
        % TimeFormat: string, optional
        %   Datetime format. Default is empty in which case an attempt will
        %   be made to auto-identify the format
        % TimeZone: string, optional
        %   Optional timezone. Default is 'UTC'
        % MetadataFields: (1 x M) string array, optional
        %   List of column names to be added to the measurement's metadata
            
            % Call SuperClass method
            this@MeasurementReaderX();
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
        
        function [measurements, timestamp] = read(this)
            
            rowInds = this.RowIndex_;
            time = this.TimeValues_{this.RowIndex_};
            if(~isempty(this.TimeFormat))
              timestamp = datetime(this.TimeValues_{this.RowIndex_},...
                                   'InputFormat',this.TimeFormat,...
                                   'TimeZone',this.TimeZone);
            else
              timestamp = datetime(this.TimeValues_{this.RowIndex_});
            end
            while(time == this.TimeValues_{this.RowIndex_})
                this.RowIndex_ = this.RowIndex_ + 1;
                rowInds(end+1) = this.RowIndex_;
            end
            stateVectors = this.Table_{rowInds, this.StateFields};              
            
            measurements = MeasurementListX(stateVectors', timestamp);
            this.RowIndex_ = this.RowIndex_ + 1;
        end
    end
end