classdef MeasurementModelX < BaseX
% MeasurementModelX Abstract class
%
% Summary of MeasurementModelX
%  This is the base class for all TrackingX measurement models. Must be used as
%  the superclass of any custom defined TrackingX measurement model. 
%
% MeasurementModelX Properties:
%   + NumStateDims   The number of state dimensions.
%   + NumMeasDims     The number of measurement dimensions.
%
% MeasurementModelX Methods:
%   ~ feval     - Equivalent to applying the measurement model function 
%   ~ random    - Measurement noise sample generator function, i.e. v_k ~ random(~)
%   ~ pdf       - Function to evaluate the probability p(y_t|x_t) of 
%                 a (set of) measurements, given a (set of) state vector(s)
%                 e.g. pdf = @(yt,xt) mvnpdf(yt,xt,R);
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility of any custom MeasurementModelX with TrackingX.
%       
% For examples see also the source of LinearGaussianX
    
    properties (Abstract)
        NumStateDims
        NumMeasDims
        Mapping
    end
    
    methods (Abstract)
        feval(this);
        pdf(this);
        random(this);
    end
    
    methods
        function this = MeasurementModelX(varargin)
        % MeasurementModelX Constructor method
        %   
        % DESCRIPTION:
        % * MeasurementModelX() Returns an measurement model object handle
        %
        % See also CONSTANTVELOCITYMODELX, CONSTANTHEADINGMODELX.
        
            this@BaseX();     
        end  
    end
end