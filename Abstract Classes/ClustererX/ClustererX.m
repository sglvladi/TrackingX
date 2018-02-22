classdef (Abstract) ClustererX < BaseX  % Extends trackingX.BaseX
% ClustererX Abstract class
%
% Summary of ClustererX:
% This is the base class for all TrackingX clusterers.
% Any custom defined Clusterer should be derived from this ClustererX base class. 
%
% ClustererX Properties:
%   None
%
% ClustererX Methods:
%   + ClustererX - Constructor method
%
% (+) denotes puplic properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        ClusterList
        UnassocTrackInds
    end
    
    methods (Abstract)
        cluster(this);
    end
    methods
        function this = ClustererX(varargin)
        % CLUSTERERX Constructor method
        %   
        % DESCRIPTION: 
        % * ClustererX() returns a "ClustererX" object handle
            
        end
    end
end