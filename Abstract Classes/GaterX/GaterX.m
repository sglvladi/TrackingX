classdef (Abstract) GaterX < BaseX  % Extends trackingX.BaseX
% GaterX Abstract class
%
% Summary of GaterX:
% This is the base class for all TrackingX gaters.
% Any custom defined Gater should be derived from this GaterX base class. 
%
% GaterX Properties:
%   + TrackList - A (1-by-NumTracks) cell array of TrackX objects, where
%                 NumTracks denotes the number of tracks
%   + MeasurementList - A (NumObsDims-by-NumMeasurements) matrix of measurements,
%                       where NumObsDims denotes the dimentionality of
%                       measurements and NumMeasurements is the number of
%                       measurements
%   + GateVolumes - A (1-by-NumTracks) vector of the gate volumes computed
%                   by the Gater when GaterX.gate() is called
%   + ValidationMatrix - A (NumTracks-by-NumMeasurements) matrix containing
%                        ones for entries representing valid association
%                        hypotheses between tracks and measurements
%
% GaterX Methods:
%   + GaterX - Constructor method
%
% (+) denotes puplic properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        TrackList
        MeasurementList
        ValidationMatrix
        GateVolumes
    end
    
    methods (Abstract)
        gate(this);
    end
    methods
        function this = GaterX(varargin)
        % GATERX Constructor method
        %   
        % DESCRIPTION: 
        % * GaterX() returns a "GaterX" object handle
            this@BaseX();
        end
    end
end