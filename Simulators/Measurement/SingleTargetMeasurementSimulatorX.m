classdef SingleTargetMeasurementSimulatorX < MeasurementSimulatorX
% SingleTargetMeasurementSimulatorX class
% 
% Summary of SingleTargetMeasurementSimulatorX:
%   Class implementation of a single-target measurement simulator
%
% SingleTargetMeasurementSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by the simulator
%
% SingleTargetMeasurementSimulatorX Methods:
%   + SingleTargetMeasurementSimulatorX - Constructor method
%   + simulate() - Generate measurements from ground-truth data
%
% (+) denotes puplic properties/methods
% 
% See also MeasurementModelX and ClutterModelX template classes
    
    methods
        function this = SingleTargetMeasurementSimulatorX(varargin)
        % SingleTargetMeasurementSimulatorX Construct a measurement simulator object
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Measurement - Object handle to MeasurementModelX SubClass 
        %       + Clutter - Object handle to ClutterModelX SubClass 
        %       + Detection - Object handle to DetectionModelX SubClass
            
            % Call super-class
            this@MeasurementSimulatorX(varargin{:});
            
        end
        
        function MeasurementScans = simulate(this, Track, varargin)
        % simulate Simulate measurements given a set of ground truth data
        %
        % Parameters
        % ----------
        % Track: GroundTruthTrackX
        %   An object containing ground-truth information for a single-target. 
        %   Detections will be generated on the basis of the Trajectory 
        %   property for each track.
        % 
        % Returns
        % -------
        % MeasurementScans: (1 x NumScans) MeasurementListX array
        %   An object array, whose elements contain the simulated measurement
        %   scan at each timestep.
            
            % Dummy Detection probability
            MeasurementScans = MeasurementListX.empty();
            
            % Initialise storage
            numTimesteps = numel(Track.Trajectory);
            
            for k = 1:numTimesteps
                
                targetState = Track.Trajectory(k);
                timestamp = targetState.Timestamp;
                
                % Compute number of tracks and clutter measurements
                detectionProbability = 1;
                if ~isempty(this.Model.Clutter)
                    detectionProbability = this.Model.Detection.pdf(targetState.Vector);
                end
                    
                targetDetected = binornd(1,detectionProbability);
                numClutter = 0;
                if ~isempty(this.Model.Clutter)
                    numClutter = this.Model.Clutter.random('cardinality');
                end
                numMeasurements = targetDetected + numClutter;
                
                % Empty MeasurementX array
                measurements = MeasurementX.empty(0, numMeasurements);

                % Generate true measurement
                if targetDetected
                    measurementVector = this.Model.Measurement.feval(targetState.Vector,true);
                    measurements(end+1) = MeasurementX(measurementVector, timestamp, this.Model);
                end

                % Generate clutter measurements
                if(numClutter>0)
                    clutterVectors = this.Model.Clutter.random(numClutter);
                    for clutterVector = clutterVectors
                        measurements(end+1) = MeasurementX(clutterVector, timestamp, this.Model);
                    end
                end
                
                MeasurementScans(k) = MeasurementListX(measurements);
            end

        end
    end
end
