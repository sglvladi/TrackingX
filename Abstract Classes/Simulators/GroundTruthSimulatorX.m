classdef (Abstract) GroundTruthSimulatorX < BaseX
% GroundTruthSimulatorX class
% 
% Summary of GroundTruthSimulatorX:
%   This is the base class implementation of a target ground-truth simulator
%
% GroundTruthSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by 
%             the simulator
%
% GroundTruthSimulatorX Methods:
%   + GroundTruthSimulatorX - Constructor method
%   + simulate() - Generate ground-truth data
%
% (+) denotes puplic properties/methods
% 
% See also GroundTruthModelX and ClutterModelX
    
    properties 
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Transition - Object handle to TransitionModelX SubClass 
        %       + Birth - Object handle to BirthModelX SubClass, optional 
        %       + Survival - Object handle to SurvivalModelX SubClass, optional
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
        function this = GroundTruthSimulatorX(varargin)
        % GroundTruthSimulatorX Construct a measurement simulator object
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Transition - Object handle to TransitionModelX SubClass 
        %       + Birth - Object handle to BirthModelX SubClass, optional 
        %       + Survival - Object handle to SurvivalModelX SubClass, optional
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
