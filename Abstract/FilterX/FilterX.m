classdef (Abstract) FilterX < BaseX  % Extends trackingX.BaseX
% FilterX Abstract class
%
% Summary of FilterX:
% This is the base class for all Tracking filters.
% Any custom defined Filters should be derived from this FilterX base class. 
%
% FilterX Properties:
%   - State   = Object handle to a TrackingX.State
%   - Model   = Object handle to StateSpaceModelX Sub/class
%
%   (*) Signifies properties necessary to instantiate a class object
%
% FilterX Methods:
%    FilterX    - Constructor method
%    predict    - Performs filter prediction step
%    update     - Performs filter update/correction step
%    iterate    - Performs a complete filter iteration (Predict & Update)
%    smooth     - Performs smoothing on a provided set of estimates
% 
% See also DynamicModel, ObservationModel and ControlModel
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        Model
    end
    
    methods (Abstract)
        initialise(this);
        predict(this);
        update(this);
    end
    methods
        function this = FilterX(varargin)
        % FILTER Constructor method
        %   
        % DESCRIPTION: 
        % * FilterX("Model", ssm) returns a "FilterX" object handle, configured
        %   to operate in the provided ssm "StateSpaceModel" object handle.
        % * FilterX(config), returns a "FilterX" object handle, configured
        %   to operate in the provided config.Model "StateSpaceModel" object handle.
        %   
        %  See also predict, update, iterate, smooth.
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    if (isfield(varargin{1},'Model'))
                        this.Model = varargin{1}.Model;
                    end
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addOptional('Model',[]);
            parser.parse(varargin{:});
            
            if(~isempty(parser.Results.Model))
                this.Model = parser.Results.Model;
            end         
        end
    end
end