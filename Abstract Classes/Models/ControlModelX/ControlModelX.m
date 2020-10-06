%3CmxGraphModel%3E%3Croot%3E%3CmxCell%20id%3D%220%22%2F%3E%3CmxCell%20id%3D%221%22%20parent%3D%220%22%2F%3E%3CmxCell%20id%3D%222%22%20value%3D%22%2B%20NumStateDims%3A%20int%22%20style%3D%22text%3Bhtml%3D1%3BstrokeColor%3Dnone%3BfillColor%3Dnone%3Balign%3Dleft%3BverticalAlign%3Dtop%3BspacingLeft%3D4%3BspacingRight%3D4%3BwhiteSpace%3Dwrap%3Boverflow%3Dhidden%3Brotatable%3D0%3Bpoints%3D%5B%5B0%2C0.5%5D%2C%5B1%2C0.5%5D%5D%3BportConstraint%3Deastwest%3BfontStyle%3D2%22%20vertex%3D%221%22%20parent%3D%221%22%3E%3CmxGeometry%20x%3D%22-148.5%22%20y%3D%222198.136474609375%22%20width%3D%22160%22%20height%3D%2224%22%20as%3D%22geometry%22%2F%3E%3C%2FmxCell%3E%3C%2Froot%3E%3C%2FmxGraphModel%3Eclassdef ControlModelX < BaseX
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