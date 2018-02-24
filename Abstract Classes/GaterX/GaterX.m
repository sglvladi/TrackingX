classdef (Abstract) GaterX < BaseX  % Extends trackingX.BaseX
% GaterX Abstract class
%
% Summary of GaterX:
% This is the base class for all TrackingX resamplers.
% Any custom defined Gater should be derived from this GaterX base class. 
%
% GaterX Properties:
%   None
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
            
        end
    end
end