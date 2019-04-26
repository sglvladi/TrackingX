classdef (Abstract) MeasurementSimulatorX < BaseX
% MeasurementSimulatorX class
% 
% Summary of MeasurementSimulatorX:
%   This is the base class implementation of a measurement simulator
%
% MeasurementSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by 
%             the simulator
%
% MeasurementSimulatorX Methods:
%   + MeasurementSimulatorX - Constructor method
%   + simulate() - Generate measurements from ground-truth data
%
% (+) denotes puplic properties/methods
% 
% See also MeasurementModelX and ClutterModelX
    
    properties 
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Measurement - Object handle to MeasurementModelX SubClass 
        %       + Clutter - Object handle to ClutterModelX SubClass 
        %       + Detection - Object handle to DetectionModelX SubClass
        Model
    end
    
    methods (Access=protected)
        function initialise_(this, config)
            if(isfield(config,'Model'))
                this.Model = config.Model;
            end
        end
    end
    
    methods
        function this = MeasurementSimulatorX(varargin)
        % MeasurementSimulatorX Construct a measurement simulator object
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Measurement - Object handle to MeasurementModelX SubClass 
        %       + Clutter - Object handle to ClutterModelX SubClass 
        %       + Detection - Object handle to DetectionModelX SubClass
        %       + ...
            
            % First check to see if a structure was received
            if isstruct(varargin{1})
                config = varargin{1};
                this.initialise_(config);
                return;
            elseif isa(varargin{1},'StateSpaceModelX')
                this.Model = varargin{1};
                return
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
        end
    end
    
    methods (Abstract)
        simulate(this); 
    end
end
