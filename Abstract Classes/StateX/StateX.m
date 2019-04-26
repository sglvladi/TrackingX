classdef (Abstract) StateX < BaseX & dynamicprops
% StateX class
%
% Summary of StateX:
% Class implementation of the primitive State type.
% 
% See also dynamicprops
%
% November 2018 Lyudmil Vladimirov, University of Liverpool.
    
    properties
        % Timestamp: datetime 
        %   An optional timestamp that can be attached to a StateX object.
        Timestamp = []
        Metadata = struct()
    end

    properties (Dependent)
        % Vector: The state vector
        Vector
        
        % NumDims
        NumDims
    end       

    properties (Access = protected)
        % Vector_: Internal storage for the state vector
        Vector_ = []
    end
    
    methods
        function this = StateX(varargin)
        % StateX Base constructor
        %
        % Parameters
        % ----------
        % Vector: numeric
        %   The state vector
        % Timestamp: datetime
        %   An optional timestamp that can be attached to a StateX object.
            
            switch(nargin)
                case(1)
                    if isnumeric(varargin{1})
                        this.Vector_ = varargin{1};
                    else
                        this.Timestamp = varargin{1};
                    end
                otherwise
                    this.Vector_ = varargin{1};
                    this.Timestamp = varargin{2};
            end
        end
        
        function stateVector = get.Vector(this)
        %get.Vector Custom getter for Vector property
            stateVector = getVector(this);
        end

        function set.Vector(this, stateVector)
        %set.Vector Custom setter for Vector property
            this.Vector_ = setVector(this, stateVector);
        end
        
        function numDims = get.NumDims(this)
        %get.NumDims Custom getter for NumDims property
            numDims = getNumDims(this);
        end
    end
    
    methods (Access=protected)
        function stateVector = getVector(this)
            stateVector = this.Vector_;
        end
        function Vector = setVector(this, stateVector)
            Vector = stateVector;
        end
        function numDims = getNumDims(this)
            numDims = size(this.Vector,1);
        end
    end
    
    methods (Static)
        function [timestamp, other] = extract_timestamp(varargs)
            timestamp = [];
            other = varargs;
            if numel(varargs)
                if isdatetime(varargs{end})
                timestamp = varargs{end};
                other = varargs(1:end-1);
                elseif isa(varargs{1},'StateX')
                    timestamp = varargs{1}.Timestamp;
                    other = varargs;
                end
            end
        end
    end
end

