classdef MultiTargetMeasurementSimulatorX < BaseX
% MultiTargetMeasurementSimulatorX class
% 
% Summary of MultiTargetMeasurementSimulatorX:
%   This is a class implementation of a simple measurement simulator
%
% MultiTargetMeasurementSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by the simulator
%
% MultiTargetMeasurementSimulatorX Methods:
%   + MultiTargetMeasurementSimulatorX - Constructor method
%   + simulate() - Generate measurements from ground-truth data
%
% (+) denotes puplic properties/methods
% 
% See also MeasurementModelX and ClutterModelX template classes
    
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
        function this = MultiTargetMeasurementSimulatorX(varargin)
        % MultiTargetMeasurementSimulatorX Construct a measurement simulator object
        %
        % Parameters
        % ----------
        % Model: StateSpaceModelX object
        %   A state-space model that should define the following models:    
        %       + Measurement - Object handle to MeasurementModelX SubClass 
        %       + Clutter - Object handle to ClutterModelX SubClass 
        %       + Detection - Object handle to DetectionModelX SubClass
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
        
        function MeasurementScans = simulate(this, varargin)
        % simulate Simulate measurements given a set of ground truth data
        %
        % Parameters
        % ----------
        % TrackStateSequence: (1 x NumTimesteps) cell array
        %   A cell array, whose elements contain the ground-truth state 
        %   information for each target, at each timestep. 
        % 
        % Returns
        % -------
        % MeasurementScans: (1 x NumScans) MeasurementListX array
        %   An object array, whose elements contain the simulated measurement
        %   scan at each timestep.
            
            % Dummy Detection probability
            detectionProbability = 1; %this.DetectionProbability;
            
            % Initialise storage
            GroundTruth = varargin{1};
            numTimesteps = numel(GroundTruth);
            MeasurementScans = MeasurementListX.empty();

            for k = 1:numTimesteps
                
                % Empty MeasurementX array
                measurements = MeasurementX.empty();
                
                % Compute number of tracks and clutter measurements
                trackStates = GroundTruth{k};
                numTracks = numel(trackStates);
                if numTracks
                    trackStateVectors = [trackStates.Vector];
                    timestamp = trackStates(1).Timestamp;

                    detectedTargetInds = find(binornd(1,detectionProbability,1,numTracks));
                    numDetectedTargets = size(detectedTargetInds,2);

                    numClutter = 0;
                    if ~isempty(this.Model.Clutter)
                        numClutter = this.Model.Clutter.random('cardinality');
                    end
                    numMeasurements = numDetectedTargets + numClutter;

                    % Generate true measurements
                    detectedTargetStateVectors = trackStateVectors(:, detectedTargetInds);

                    % Generate true measurement
                    measurementVectors = zeros(this.Model.Measurement.NumMeasDims,numMeasurements);
                    if numDetectedTargets
                        measurementVectors(:, 1:numDetectedTargets) = this.Model.Measurement.feval(detectedTargetStateVectors,true);
                    end

                    % Generate clutter measurements
                    if(numClutter>0)
                        measurementVectors(:, numDetectedTargets+1:end) = this.Model.Clutter.random(numClutter);
                    end

                    MeasurementScans(k) = MeasurementListX(measurementVectors, timestamp);
                end
            end

        end
    end
end
