classdef (Abstract) FilterX < BaseX % Extends trackingX.BaseX
% FilterX Abstract class
%
% Summary of FilterX:
% This is the base class for all TrackingX filters.
% Any custom defined Filter should be derived from this FilterX base class.
%
% By default, FilterX objects are designed with the intention to be utilised 
% as (Markov-Chain) State Estimators that opererate on a given problem as 
% deterministic Finite State Machines (FSMs), with self-contained memory.
% In other words, it is expected that for most problems, FilterX objects will 
% go through any heavy parameterisation once, after which they will be used
% to recursively execute a Prediction-Update State Estimation loop. 
% FilterX objects aim to preserve all the inter-state information they require 
% such that they can execute the next FSM step, without being provided with
% information they already know.    
%
%
% FilterX Properties:
%   + Model      - Object handle to StateSpaceModelX class, containing model
%                  parameterisations. Normally set once during filter initialisation and
%                  accessed later on.
%   + Prior      - Structure containing prior information (e.g. pdf)
%                  Normally set once during filter initialisation and
%                  accessed later on.
%   + Prediction - Structure containing prediction information (e.g. pdf)
%                  Data is over-written on every prediction step and is 
%                  subsequently read by the update (and later) step(s)
%   + Posterior  - Structure containing posterior information (e.g. pdf)
%                  Data is over-written on every update step and is subsequently
%                  utilised by 
%   + MeasurementList   - A (matrix of) column vector(s) representing the
%                         measurement(s) to be used during the filter update step
%
%   (*) Signifies properties necessary to instantiate a class object
%
% FilterX Methods:
%   + FilterX    - Constructor method
%   + predict    - Performs filter prediction step
%   + update     - Performs filter update/correction step
%
% (+) denotes puplic properties/methods
% 
% See also TransitionModelX, MeasurementModelX and ControlModelX
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.

    properties
        Model
        MeasurementList
    end
    
    properties (Access = protected)
        
    end
    
    methods (Abstract)
        predict(this);
        update(this);
    end
    
    methods (Access = protected)
        function initialise_(this, config)
            if (isfield(config,'Model'))
                this.Model = config.Model;
            end
        end
    end
          
    methods
        function this = FilterX(varargin)
        % FilterX Constructor method
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %  Object handle to a given state-space model
        %
        % Usage
        % -----
        % * FilterX(__, Name, Value) instantiates an object handle, configured 
        %   the parameters specified by one or more Name,Value pair arguments. 
        % * FilterX(config) instantiates an object handle configured with  
        %   the parameters specified inside the 'config' structure, whose  
        %   fieldnames correspond to the parameter names as given in the 
        %   Parameters section above.
        %   
        %  See also predict, update, iterate, smooth.
            
            if(nargin==0)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==1)
                if(isstruct(varargin{1}))
                    this.initialise_(varargin{1});
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            this.initialise_(parser.Unmatched); 
        end
        
        function initialise(this, varargin)
        % initialise Initialisation method. Reset and re-configure the filter.
        %   
        % Parameters
        % ----------
        % Model: StateSpaceModelX
        %  Object handle to a given state-space model.
        %
        % Usage
        % -----
        % * filter.initialise(__, Name, Value) instantiates an object handle,  
        %   configured the parameters specified by one or more Name,Value
        %   pair arguments. 
        % * filter.initialise(config) instantiates an object handle configured  
        %   with the parameters specified inside the 'config' structure, whose  
        %   fieldnames correspond to the parameter names as given in the 
        %   Parameters section above.
        %   
        %  See also predict, update, iterate, smooth.
            
            if(nargin==1)
                return;
            end
            
            % First check to see if a structure was received
            if(nargin==2)
                if(isstruct(varargin{1}))
                    this.initialise_(varargin{1});
                    return;
                end
            end
            
            % Otherwise, fall back to input parser
            parser = inputParser;
            parser.KeepUnmatched = true;
            parser.parse(varargin{:});
            this.initialise_(parser.Unmatched); 
        end
        
        % ===============================>
        % ACCESS METHODS
        % ===============================>       
        function set.Model(this,newModel)
            this.Model = setModel(this,newModel);
        end
        
        function set.MeasurementList(this,newMeasurementList)
            this.MeasurementList = setMeasurementList(this,newMeasurementList);
        end
    end

    methods (Access = protected)
        % ===============================>
        % ACCESS METHOD HANDLES
        % ===============================>
        function Prior = setPrior(this,newPrior)
            Prior = newPrior;
        end
        function Prediction = setPrediction(this,newPrediction)
            Prediction = newPrediction;
        end
        function Posterior = setPosterior(this,newPosterior)
            Posterior = newPosterior;
        end
        function Model = setModel(this,newModel)
            Model = newModel;
        end
        function measurementList = setMeasurementList(this,newMeasurementList)
            measurementList = newMeasurementList;
        end
    end
end