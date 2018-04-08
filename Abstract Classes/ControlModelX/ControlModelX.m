classdef ControlModelX < BaseX
% ControlModelX Abstract class
%
% Summary of ControlModelX
%  This is the base class for all TrackingX control models. Must be used as
%  the superclass of any custom defined TrackingX control model. 
%
% ControlModelX Properties:
%   + NumStateDims - The number of state dimensions.
%   + NumControlDims - The number of control input dimensions.
%
% ControlModelX Methods:
%   ~ beval - Equivalent to applying the model control function 
%   ~ random - Process noise sample generator function, i.e. w_k ~ random(~)
%   ~ pdf - Function to evaluate the probability p(x_t|u_t) of a (set of)
%           state(s), given a (set of) control input(s)
%
% (+) denotes puplic properties/methods
% (~) denotes abstract properties/methods
%
% The above parameters and methods are accessed by the majority of 
% existing TrackingX library components and are COMPULSORY to guarantee 
% compatibility of any custom ControlModelX with TrackingX.
%       
% For examples see also the source of See also LINEARNOISELESSCTRMODELX_1D
    
    properties
        NumStateDims    % The number of state dimensions.
        NumControlDims  % The number of control input dimensions.
    end
    
    methods (Abstract)
        beval(this);    % Equivalent to applying the model control function
        random(this);   % Process noise sample generator function, i.e. w_k ~ random(~)
        pdf(this)       % Function to evaluate the probability p(x_t|u_t) of a (set of) state(s), given a (set of) control input(s)
    end
    
    methods
        function this = ControlModelX(varargin)
        % CONTROLMODEL Constructor method
        %
        % See also CONSTANTVELOCITYMODELX, CONSTANTHEADINGMODELX.
      
            this@BaseX(varargin);      
        end
    end
end