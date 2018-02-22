classdef ControlModelX < BaseX
% ControlModelX Abstract class
%
% Summary of ControlModelX
%  This is the base class for all TrackingX control models. Must be used as
%  the superclass of any custom defined TrackingX control model. 
%
% ControlModelX Properties:
%   + NumStateDims   The number of state dimensions.
%   + NumControlDims     The number of observation dimensions.
%
% ControlModelX Methods:
%   + feval(~)   Equivalent to applying the model transition function 
%   + random(~)  Process noise sample generator function, i.e. w_k ~ noise(~)
%   + pdf(~)     Function to evaluate the probability p(x_t|x_{t-1}) of 
%                 a set of new states, given a set of (particle) state vectors
%                   e.g. eval = @(xt,xtm1) mvnpdf(xt,xtm1,Q);
%
%  (+) denotes puplic properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility of any custom ControlModelX with TrackingX.
%       
% For examples see also the source of See also LINEARNOISELESSCTRMODEL_1D
    
    properties
    end
    
    methods
        function this = ControlModelX(varargin)
        % CONTROLMODEL Constructor method
        %   
        % DESCRIPTION:
        % * ControlModelX() Returns a "ConrolModel" object handle
        %
        % See also CONSTANTVELOCITYMODELX, CONSTANTHEADINGMODELX.
      
            this@BaseX(varargin);      
        end
    end
end