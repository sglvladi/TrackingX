classdef (Abstract) MeasurementReaderX < ReaderX
% MeasurementReaderX Abstract class
%
% Summary of MeasurementReaderX
%  This is the base class for all TrackingX measurement readers. Must be 
%  used as the superclass of any custom defined TrackingX measurement readers. 
%
% MeasurementReaderX Properties:
%   + Measurements - The measurements at a given timestep.
%
% MeasurementReaderX Methods: 
%   ~ read - Read and return measurements at current timestep
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility with the TrackingX framework.
    
    properties
        Measurements
    end
    
    methods (Abstract)
        read(this);
    end
end