classdef (Abstract) ResamplerX < BaseX & matlabshared.tracking.internal.Resampler  % Extends trackingX.BaseX
% ResamplerX Abstract class
%
% Summary of ResamplerX:
% This is the base class for all TrackingX resamplers.
% Any custom defined Resampler should be derived from this ResamplerX base class. 
%
% ResamplerX Properties:
%   None
%
% ResamplerX Methods:
%   + ResamplerX - Constructor method
%
% (+) denotes puplic properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
    end
    
    methods
        function this = ResamplerX(varargin)
        % RESAMPLERX Constructor method
        %   
        % DESCRIPTION: 
        % * ResamplerX() returns a "ResamplerX" object handle
            
        end
    end
end