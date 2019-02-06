classdef StateListX < ObjectArrayWrapperX
% StateListX class
% 
% Summary of StateListX:
% Class implementation of the primitive State List type.
    
    properties (Dependent)
        Vectors
        Timestamps
        NumStates
    end

    properties
        States
    end
    
    methods
        function this = StateListX(varargin)
        % State Base constructor
        %
        % Parameters
        % ----------
        % states: vector of StateX objects
        %   
            this@ObjectArrayWrapperX('States');
            
            if nargin && ~isempty(varargin{1})
                if isa(varargin{1},'StateListX')
                    states = varargin{1}.States;
                    this.States = states.copy();
                else
                    states = varargin{1};
                    numStates = numel(states);
                    if iscell(states)
                        for i=1:numStates
                            dStates(i) = states{i};
                        end
                    else
                        dStates = states;
                    end
                    this.States = dStates;
                end
            else 
                this.States = [];
            end
        end
        
        
        function Vectors = get.Vectors(this)
            if(isempty(this.States))
                Vectors = [];
            else
                Vectors = this.States.Vector;
            end
        end
        
        function Timestamps = get.Timestamps(this)
            if(isempty(this.States))
                Timestamps = [];
            else
                Timestamps = this.States.Timestamp;
            end
        end
        
        function numStates = get.NumStates(this)
            if(isempty(this.States))
                numStates = 0;
            else
                numStates = numel(this.States);
            end
        end
    end
    
    methods (Static)
        function [model, other] = extract_model(varargs)
            model = [];
            other = varargs;
            if numel(varargs) && isa(varargs{end}, 'StateSpaceModelX')
                model = varargs{end};
                other = varargs(1:end-1);
            end
        end
    end
end

