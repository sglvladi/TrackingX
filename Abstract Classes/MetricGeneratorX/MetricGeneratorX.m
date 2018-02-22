classdef (Abstract) MetricGeneratorX < BaseX 
% MetricGeneratorX Abstract class
%
% Summary of MetricGeneratorX:
% This is the base class for all TrackingX resamplers.
% Any custom defined Resampler should be derived from this ResamplerX base class. 
%
% MetricGeneratorX Properties:
%   None
%
% MetricGeneratorX Methods:
%   + MetricGeneratorX - Constructor method
%
% (+) denotes puplic properties/methods
%
% February 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
    end

    methods (Abstract)
        % evaluate Evaluate the correctness of a set of tracks compared to
        % the groundtrouth.
        %   METRIC = EVALUATE(THIS, TRACKLIST, GROUNDTRUTH) evaluates the
        %   correctness of the TRACKLIST compared to the GROUNDTRUTH.
        %   TRACKLIST and GROUNDTRUTH can be structured based on the needs
        %   of the respective metric generator. 

        metric = evaluate(~,TrackList,GroundTruth);
    end
end