classdef (Abstract) HypothesiserX < BaseX 
% HypothesiserX Abstract class
%
% Summary of HypothesiserX:
% This is the base class for all TrackingX hypothesisers.
% Any custom defined hypothesiser should be derived from this HypothesiserX base class. 
%
% HypothesiserX Properties:
%   None
%
% HypothesiserX Methods:
%   + HypothesiserX - Constructor method
%   ~ hypothesise - Performs the generation an evaluation of measurement
%                   hypotheses
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
    end
    
    methods (Abstract)
        hypothesise(this);
    end
    methods
        function this = HypothesiserX(varargin)
        % HYPOTHESISERX Constructor method
        %   
        % DESCRIPTION: 
        % * HypothesiserX() returns a "HypothesiserX" object handle
            
            this@BaseX();
        end
    end
end