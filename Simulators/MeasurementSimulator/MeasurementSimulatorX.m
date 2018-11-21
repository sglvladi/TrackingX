classdef MeasurementSimulatorX < BaseX
% MeasurementSimulatorX class
% 
% Summary of MeasurementSimulatorX:
%   This is a class implementation of a simple measurement simulator
%
% MeasurementSimulatorX Properties:
%   + Model - A StateSpaceModelX object defining the base models used by the simulator
%
% MeasurementSimulatorX Methods:
%   + MeasurementSimulatorX - Constructor method
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
        DetectionProbability
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
        
        function [DataList, newGroundTruth] = simulate(this, GroundTruth)
        % simulate Simulate measurements given a set of ground truth data
        %
        % Parameters
        % ----------
        % GroundTruth: (1 x NumTimesteps) cell array of (NumStateDims x NumTargets) matrices
        %   A cell array, whose elements contain the ground truth positions of targets
        %   at each timestep.
        % 
        % Returns
        % -------
        % DataList: (1 x NumTimesteps) cell array of (NumStateDims x NumMeasurements) matrices
        %   A cell array, whose elements contain the simulated measurement scan
        %   at each timestep.
        % newGroundTruth: (1 x NumTimesteps) cell array of (NumStateDims x NumTargets) matrices
        %   A copy of GroundTruth, where the state vectors of ground truth
        %   targets are transformed to conform to the specified state space
        %   model

            % Initialise storage
            numTimesteps      = numel(GroundTruth);
            DataList          = cell(1,numTimesteps);
            newGroundTruth = cell(1,numTimesteps);
            
            % Dummy Detection probability
            detectionProbability = this.DetectionProbability;

            for k = 1:numTimesteps

                % Compute number of tracks and clutter measurements
                numTrueTracks = size(GroundTruth{k},2);
                detectedTargetInds = find(binornd(1,detectionProbability,1,numTrueTracks));
                numDetectedTargets = size(detectedTargetInds,2);
                numClutter = this.Model.Clutter.random('cardinality');
                numMeasurements = numDetectedTargets + numClutter;

                % Initialise measurement list for k-th timestep
                DataList{k} = zeros(this.Model.Measurement.NumMeasDims, numMeasurements);

                % Generate true measurements
                trueTargetStates = zeros(this.Model.Measurement.NumStateDims, numDetectedTargets);
                if(size(GroundTruth{k},1) ~= this.Model.Measurement.NumStateDims)
                    trueTargetStates(this.Model.Measurement.Mapping,:) = GroundTruth{k}(:,detectedTargetInds);
                else
                    trueTargetStates = GroundTruth{k}(:,detectedTargetInds);
                end
                newGroundTruth{k}(this.Model.Measurement.Mapping,:) = GroundTruth{k}(:,detectedTargetInds);
                DataList{k}(:,1:numDetectedTargets) = this.Model.Measurement.feval(trueTargetStates,true);

                % Generate clutter measurements
                if(numClutter>0)
                    DataList{k}(:,numDetectedTargets+1:end) = this.Model.Clutter.random(numClutter);
                end
            end

        end
    end
end
