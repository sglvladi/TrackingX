classdef SingleTargetGroundTruthSimulatorX < GroundTruthSimulatorX
% SingleTargetGroundTruthSimulatorX class
% 
% Summary of SingleTargetGroundTruthSimulatorX:
%   Class implementation of a single-target measurement simulator
%
% SingleTargetGroundTruthSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by the simulator
%
% SingleTargetGroundTruthSimulatorX Methods:
%   + SingleTargetGroundTruthSimulatorX - Constructor method
%   + simulate() - Generate measurements from ground-truth data
%
% (+) denotes puplic properties/methods
% 
% See also MeasurementModelX and ClutterModelX template classes
    
    properties
        InitialState
        TimestepDuration = duration(0,0,1);
        NumTimesteps = 100;
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            initialise_@GroundTruthSimulatorX(this,config);
            if (isfield(config,'InitialState'))
                this.InitialState = config.InitialState;
            end
            if (isfield(config,'TimestepDuration'))
                this.TimestepDuration = config.TimestepDuration;
            end
            if (isfield(config,'NumTimesteps'))
                this.NumTimesteps = config.NumTimesteps;
            end
        end
    end
    methods
        function this = SingleTargetGroundTruthSimulatorX(varargin)
        % SingleTargetGroundTruthSimulatorX Construct a measurement simulator object
        %
        % Parameters
        % ---------- 
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Measurement - Object handle to MeasurementModelX SubClass 
        %       + Clutter - Object handle to ClutterModelX SubClass 
        %       + Detection - Object handle to DetectionModelX SubClass
        % InitialState: StateX object
        %   Initial state to use to initiate ground truth track
        % TimestepDuration: duration, optional
        %   Time spent between timesteps (default = duration(0,0,1))
        % NumTimesteps: scalar, optional
        %   The number of timestep to run the simulation over. (default = 100)
            
            % Call super-class
            this@GroundTruthSimulatorX(varargin{:});
             
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    config = varargin{1};
                    this.initialise_(config);
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            config = parser.Unmatched;
            this.initialise_(config);
            
        end
        
        function Track = simulate(this, varargin)
        % simulate Simulate measurements given a set of ground truth data
        %
        % Parameters
        % ----------
        % Track: GroundTruthTrackX, optional
        %   A prior ground-truth track. If provided, the simulation
        %   simulate the track forward.
        % 
        % Returns
        % -------
        % Track: GroundTruthTrackX
        %   The simulated track.
            
            switch (nargin)
                case 1
                    Track = GroundTruthTrackX(this.InitialState);
                case 2
                    Track = varargin{1};
            end
            
            if ~isdatetime(Track.State.Timestamp)
                Track.Trajectory(end).Timestamp = datetime();
            end
            
            for i = 1:this.NumTimesteps
                new_timestamp = Track.State.Timestamp + this.TimestepDuration;
                state_vector = ...
                    this.Model.Transition.feval(Track.State.Vector, true, ...
                                                this.TimestepDuration);
                Track.Trajectory(end+1) = GroundTruthStateX(state_vector,... 
                                                            new_timestamp);
            end
        end
    end
end
