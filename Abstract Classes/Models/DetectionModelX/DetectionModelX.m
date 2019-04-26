classdef (Abstract) DetectionModelX < BaseX
% DetectionModelX Abstract class
%
% Summary of DetectionModelX
%  This is the base class for all TrackingX detection models. Must be used as
%  the superclass of any custom defined TrackingX detection model. 
%
% DetectionModelX Properties:
%   + NumMeasDims - The number of observation dimensions.
%
% DetectionModelX Methods: 
%   ~ pdf    - Function to evaluate the probability of detection;
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility with the TrackingX framework.
%
% For examples see also the source of CONSTANTDETECTIONPROBABILITYX
    
    properties
    end
    
    methods (Abstract)
        pdf(this);      % Process noise sample generator function
    end
    
    methods
        function this = DetectionModelX()
        % DetectionModelX Constructor method
        %
        % See also CONSTANTDETECTIONPROBABILITYX
        
            this@BaseX();
        end        
    end
end