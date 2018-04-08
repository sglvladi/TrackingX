classdef (Abstract) TrackDeleterX < BaseX 
% TrackInitiatorX Abstract class
%
% Summary of TrackDeleterX:
% This is the base class for all TrackingX track deleters.
% Any custom defined track deleter should be derived from this TrackDeleterX base class. 
%
% TrackDeleterX Properties:
%   None
%
% TrackDeleterX Methods:
%   + TrackDeleterX - Constructor method
%
% (+) denotes public properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
    end
    
    methods
        function this = TrackDeleterX(varargin)
        % TRACKINITIATORX Constructor method
        %   
        % DESCRIPTION: 
        % * TrackDeleterX() returns a "TrackDeleterX" object handle
            
            this@BaseX();
        end
    end
end