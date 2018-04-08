classdef StateSpaceModelX < BaseX
% StateSpaceModelX Abstract class
%
% Summary of StateSpaceModelX:
%  This is the base class for all TrackingX state-space models.
%  State-space model brings together and provides a common interface to all
%  of 2/3 underlying models, i.e.:
%       * Dynamic/Motion model, making use of DynamicModel.
%       * Observation/Measurement model, making use of ObservationModel.
%       * ( Control model, making use of ControlModel ) - Optional!
%
% StateSpaceModelX Properties:
%   - Dyn     Object handle to a DynamicModel subclass
%   - Obs     Object handle to an ObservationModel subclass
%   - Ctr     Object handle to a ControlModel subclass
%
% StateSpaceModelX Methods:
%    StateSpaceModelX    - Constructor method
% 
% See also DynamicModel, ObservationModel and ControlModel template classes
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        Dyn
        Obs
        Ctr
    end
    
    methods
        function this = StateSpaceModelX(varargin)
        % STATESPACEMODEL Constructor method
        %   
        % DESCRIPTION:
        % * StateSpaceModelX(ssm,Dyn,Obs) configures and returns the object 
        %   ssm given the preconfigured dynamic Dyn and observation Obs models.       
        % * StateSpaceModelX(config); where "config" is a structure MUST contain 
        %   the required fields "config.Dyn" and "config.Obs" as objects 
        %   handles, while it can optionally contain "config.Ctr"
        % * StateSpaceModelX(___,Name,Value) instantiates an object 
        %   handle, configured with options specified by one or more Name,Value
        %   pair arguments.
        %
        % PARAMETERS:
        % * Dyn - (Required) DynamicModel subclass object handle
        % * Obs - (Required) ObservationModel subclass object handle
        % * Ctr - (Optional) ControlModel subclass object handle
        %   
        %  See also DynamicModel, ObservationModel and ControlModel.
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    this.Dyn = varargin{1}.Dyn;
                    this.Obs = varargin{1}.Obs;
                    if (isfield(varargin{1},'Ctr'))
                        this.Ctr = varargin{1}.Ctr;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addOptional('Dyn',@(x) isa(x,'DynamicModel'));
            parser.addOptional('Obs',@(x) isa(x,'ObservationModel'));
            parser.addParameter('Ctr',[]);
            parser.parse(varargin{:});
            
            this.Dyn = parser.Results.Dyn;
            this.Obs = parser.Results.Obs;
            
            if(this.Dyn.NumStateDims ~= this.Obs.NumStateDims)
                error('The state dimensions of the Dyn and Obs models do not agree!'); 
            end
            
            if(~isempty(parser.Results.Ctr))
                if(this.Dyn.StateDim ~= this.Ctr.StateDim)
                    error('The state dimensions of the Dyn and Ctr models do not agree!'); 
                end
                this.Ctr = parser.Results.Ctr;
            end       
        end
    end
end