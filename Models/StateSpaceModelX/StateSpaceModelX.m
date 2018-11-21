classdef StateSpaceModelX < BaseX
% StateSpaceModelX Abstract class
%
% Summary of StateSpaceModelX:
%  This is the base class for all TrackingX state-space models.
%  State-space model brings together and provides a common interface to all
%  of 2/3 underlying models, i.e.:
%       * Transition model, making use of TransitionModelX.
%       * Measurement model, making use of MeasurementModelX.
%       * ( Control model, making use of ControlModel ) - Optional!
%
% StateSpaceModelX Properties:
%   - Transition     Object handle to a TransitionModelX subclass
%   - Measurement    Object handle to an MeasurementModelX subclass
%   - Control        Object handle to a ControlModelX subclass
%   - Clutter        Object handle to a ClutterModelX subclass
%   - Detection      Object handle to an DetectionModelX subclass
%   - Birth          Object handle to a BirthModelX subclass
%
% StateSpaceModelX Methods:
%    StateSpaceModelX    - Constructor method
% 
% See also TransitionModelX, MeasurementModelX and ControlModel template classes
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        Transition
        Measurement
        Control
        Clutter
        Birth
        Detection
    end
    
    methods
        function this = StateSpaceModelX(varargin)
        % STATESPACEMODEL Constructor method
        %
        % Parameters
        % ==========
        % Transition: TransitionModelX subclass object handle
        %   A transition model
        % Measurement: MeasurementModelX subclass object handle
        %   A measurement model
        % Control: ControlModelX subclass object handle, optional
        %   
        % Usage
        % =====
        % * StateSpaceModelX(ssm,Transition,Measurement) configures and returns the object 
        %   ssm given the preconfigured transition and measurement models.       
        % * StateSpaceModelX(config); where "config" is a structure MUST contain 
        %   the required fields "config.Transition" and "config.Measurement" as objects 
        %   handles, while it can optionally contain "config.Control"
        % * StateSpaceModelX(___,Name,Value) instantiates an object 
        %   handle, configured with options specified by one or more Name,Value
        %   pair arguments.
        %   
        %  See also TransitionModelX, MeasurementModelX and ControlModel.
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    this.Transition = varargin{1}.Transition;
                    this.Measurement = varargin{1}.Measurement;
                    if (isfield(varargin{1},'Control'))
                        this.Control = varargin{1}.Control;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addOptional('Transition',[], @(x) isa(x,'TransitionModelX'));
            parser.addOptional('Measurement',[], @(x) isa(x,'MeasurementModelX'));
            parser.addParameter('Control',[]);
            parser.addParameter('Clutter',[]);
            parser.addParameter('Birth',[]);
            parser.addParameter('Detection',[]);
            parser.parse(varargin{:});
            
            this.Transition = parser.Results.Transition;
            this.Measurement = parser.Results.Measurement;
            
            if(this.Transition.NumStateDims ~= this.Measurement.NumStateDims)
                error('The state dimensions of the Transition and Measurement models do not agree!'); 
            end
            
            if(~isempty(parser.Results.Control))
                if(this.Transition.StateDim ~= this.Ctr.StateDim)
                    error('The state dimensions of the Transition and Ctr models do not agree!'); 
                end
                this.Control = parser.Results.Control;
            end
            this.Clutter = parser.Results.Clutter;
            this.Birth = parser.Results.Birth;
            this.Detection = parser.Results.Detection;
        end
    end
end