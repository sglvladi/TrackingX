classdef ObjectArrayWrapperX < BaseX
% ObjectArrayWrapperX class
%
% Summary of ObjectArrayWrapperX:
% Wrapper class for object arrays, that overloads the subsref and subsasgn 
% methods to provide specialised subscripted reference and assignment to a
% single class property.
% 
% More specifically, the modifications performed allow the ObjectArrayWrapperX
% class to provide access to a single "protected" Object Array, named "Array_"
% identical to using curly braces in cell arrays.  
%
% ObjectArrayWrapperX Methods:
%   + ObjectArrayWrapper  - Constructor method
%
% See also MesurementListX, KalmanFilerX.

    properties (Access = protected)
        % arrayName_ : char
        %   The name of the underlying object array
        arrayName_
    end

    methods

        function this = ObjectArrayWrapperX(varargin)
        % ObjectArrayWrapperX Constructor
        %
        % This base constructor MUST be invoked for the sub-class to perform as
        % intended.
        %
        % Parameters
        % ----------
        % arrayName: char
        %   The name of the (public) array to link 
            this.arrayName_ = varargin{1};
        end

        function varargout = subsref(this,S)
        % subsref Overoaded subscripted reference operator
        %
            switch S(1).type
                case '.'
                    % Use built-in for any other expression
                    [varargout{1:nargout}] = builtin('subsref',this,S);
                case '()'
                    % Use built-in for any other expression
                    try
                        [varargout{1:nargout}] = builtin('subsref',this,S);
                    catch
                    end
                case '{}'
                    P = S; P(1).type = '()';
                    if length(P) == 1 || length(P) == 2 && strcmp(S(2).type,'.')
                    % Overloaded referencing to this{indices} &
                    % this{indices}.PropertyName subscripts
                        [varargout{1:nargout}] = builtin('subsref',this.(this.arrayName_),P);
                    else
                    % Use built-in for any other expression
                        [varargout{1:nargout}] = builtin('subsref',this,S);
                    end
                otherwise
                    error('Not a valid indexing expression')
            end
        end
    end
end