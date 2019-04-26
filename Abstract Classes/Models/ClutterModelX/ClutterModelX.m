classdef (Abstract) ClutterModelX < BaseX
% ClutterModelX Abstract class
%
% Summary of ClutterModelX
%  This is the base class for all TrackingX clutter models. Must be used as
%  the superclass of any custom defined TrackingX clutter model. 
%
% ClutterModelX Properties:
%   + NumMeasDims - The number of observation dimensions.
%
% ClutterModelX Methods: 
%   ~ random - Process noise sample generator function
%   ~ pdf    - Function to evaluate the probability p(x_k|x_{k-1}) of 
%              a set of new states, given a set of (particle) state vectors
%              e.g. eval = @(xk,xkm1) mvnpdf(xk,xkm1,Q);
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility with the TrackingX framework.
%
% For examples see also the source of POISSONNUMBERUNIFORMPOSITIONCLUTTERMODELX
    
    properties
        NumMeasDims
    end
    
    methods (Abstract)
        pdf(this);      % Process noise sample generator function
        random(this);   % Function to evaluate the probability p(x_k|x_{k-1}) of a set of new states, given a set of (particle) state vectors
    end
    
    methods
        function this = ClutterModelX()
        % ClutterModelX Constructor method
        %
        % See also CONSTANTVELOCITYMODEL_1D, CONSTANTHEADINGMODEL.
        
            this@BaseX();
        end        
    end
end