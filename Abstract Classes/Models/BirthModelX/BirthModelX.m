classdef (Abstract) BirthModelX < BaseX
% BirthModelX Abstract class
%
% Summary of BirthModelX
%  This is the base class for all TrackingX birth models. Must be used as
%  the superclass of any custom defined TrackingX birth model. 
%
% BirthModelX Properties:
%   + NumMeasDims - The number of observation dimensions.
%
% BirthModelX Methods: 
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
% For examples see also the source of GenericBirthModelX
    
    properties
    end
    
    methods (Abstract)
        pdf(this);      % Function to evaluate the intensity of the birth model 
        random(this);   % Generate iid samples from the birth model
    end
    
    methods
        function this = BirthModelX()
        % BirthModelX Constructor method
        %
        % See also GenericBirthModelX.
        
            this@BaseX();
        end        
    end
end