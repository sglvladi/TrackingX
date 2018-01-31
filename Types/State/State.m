classdef State < Base
% STATE class
%
% Summary of State:
% This is the  
%
% Filter Properties:
%   - Params  = structure to allow for the storage of user defined parameters
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
% See also DynamicModelX, ObservationModelX and ControlModelX template classes
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        Model
        Params
    end
    
    methods (Abstract)
        predict(this);
        update(this);
        smooth(this);
    end
    methods
        function this = FilterX(varargin)
        % FILTERX Constructor method
        %   
        % INPUTS:   
        %               
        %       Optional
        %       ========
        %       Model  StateSpaceModelX SubClass instance      
        %   
        % USAGE:
        %       f = FilterX();
        %       f = FilterX("Model", sys); Specify any or none of the 3 optional
        %           initialisation variables, in any given order;
        %       f = FilterX(config); where "config" is a structure may (not necessarily)
        %           contain a StateSpaceModelX instance under "config.Model";
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
                    if (isfield(varargin{1},'customParams'))
                        this.customParams = varargin{1}.customParams;
                    end
                end
                return;
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.addParameter('Model',NaN);
            parser.addParameter('customParams',NaN);
            parser.parse(varargin{:});
            
            if(~isnan(parser.Results.Model))
                this.Model = parser.Results.Model;
            end
            
            if(~isnan(parser.Results.Params))
                this.Model.customParams = parser.Results.Params;
            end
            
%             % Add DynModel
%             if(~isfield(Init,'DynModel'))
%                 error('No DynModel provided!');
%             else
%                 this.DynModel = Init.DynModel;
%             end
%             
%             % Add ObsModel
%             if(~isfield(Init,'ObsModel'))
%                 error('No ObsModel provided!');
%             else
%                 this.ObsModel = Init.ObsModel;
%             end
%             
%             % Validate CtrModel
%             if(isfield(Init,'CtrModel'))
%                 this.CtrModel = Init.CtrModel;
%             end
%            
        end
               
        function iterate(this,varargin)
        % ITERATE - Performs a complete filter iteration (Predict & Update)
        %   
        % INPUTS:
        %       N/A 
        %
        % OUTPUTS:
        %       f.iterate();
        %
        %   See also FilterX, predict, update, smooth.
        
            this.predict();  % Predict         
            this.update();   % Update
        end

    end
end

