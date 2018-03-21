classdef (Abstract) TrackInitiatorX < BaseX 
% TrackInitiatorX Abstract class
%
% Summary of TrackInitiatorX:
% This is the base class for all TrackingX resamplers.
% Any custom defined Track Initiator should be derived from this TrackInitiatorX base class. 
%
% TrackInitiatorX Properties:
%   None
%
% TrackInitiatorX Methods:
%   + TrackInitiatorX - Constructor method
%
% (+) denotes public properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
    end
    
    methods
        function this = TrackInitiatorX(varargin)
        % TRACKINITIATORX Constructor method
        %   
        % DESCRIPTION: 
        % * TrackInitiatorX() returns a "TrackInitiatorX" object handle
            
            this@BaseX();
        end
    end
end