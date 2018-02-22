classdef ObservationModelX < BaseX
% ObservationModelX Abstract class
%
% Summary of ObservationModelX
%  This is the base class for all TrackingX observation models. Must be used as
%  the superclass of any custom defined TrackingX observation model. 
%
% ObservationModelX Properties:
%   + NumStateDims   The number of state dimensions.
%   + NumObsDims     The number of observation dimensions.
%
% ObservationModelX Methods:
%   + feval(~)   Equivalent to applying the model transition function 
%   + random(~)  Process noise sample generator function, i.e. w_k ~ noise(~)
%   + pdf(~)     Function to evaluate the probability p(x_t|x_{t-1}) of 
%                 a set of new states, given a set of (particle) state vectors
%                   e.g. eval = @(xt,xtm1) mvnpdf(xt,xtm1,Q);
%
% (+) denotes puplic properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility of any custom ObservationModelX with TrackingX.
%       
% For examples see also the source of LINEARGAUSSOBSMODEL_1D
    
    properties
        NumStateDims
        NumObsDims
    end
    
    methods (Abstract)
        heval(this);
        pdf(this);
        random(this);
    end
    
    methods
        function this = ObservationModelX(varargin)
        % OBSERVATIONMODEL Constructor method
        %   
        % DESCRIPTION:
        % * ObservationModelX() Returns an observation model object handle
        %
        % See also CONSTANTVELOCITYMODELX, CONSTANTHEADINGMODELX.
        
            this@BaseX();     
        end  
    end
end