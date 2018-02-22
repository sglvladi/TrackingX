classdef (Abstract) DynamicModelX < BaseX
% DynamicModelX Abstract class
%
% Summary of DynamicModelX
%  This is the base class for all TrackingX dynamic models. Must be used as
%  the superclass of any custom defined TrackingX dynamic model. 
%
% DynamicModelX Properties:
%   + NumStateDims   The number of state dimensions.
%
% DynamicModelX Methods:
%   + apply(~)   Equivalent to applying the model transition function 
%   + rnd(~)     Process noise sample generator function, i.e. w_k ~ noise(~)
%   + pdf(~)     Function to evaluate the probability p(x_t|x_{t-1}) of 
%                 a set of new states, given a set of (particle) state vectors
%                   e.g. eval = @(xt,xtm1) mvnpdf(xt,xtm1,Q);
%
% (+) denotes puplic properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility with the TrackingX framework.
%
% For examples see also the source of CONSTANTVELOCITYMODELX and CONSTANDHEADINGMODELX
    
    properties
        NumStateDims
    end
    
    methods (Abstract)
        feval(this);
        pdf(this);
        random(this);
    end
    
    methods
        function this = DynamicModelX()
        % DYNAMICMODEL Constructor method
        %   
        % DESCRIPTION:
        % * DynamicModelX() Returns a "DynamicModelX" object handle
        %
        % See also CONSTANTVELOCITYMODEL_1D, CONSTANTHEADINGMODEL.
        
            this@BaseX();
        end        
    end
end