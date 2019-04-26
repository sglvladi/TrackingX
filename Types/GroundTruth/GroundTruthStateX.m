classdef GroundTruthStateX < StateX
% GroundTruthStateX class
% 
% Summary of GroundTruthStateX:
% Class implementation of the primitive Measurement type.
    
    properties
        % Tag: TagX
        %   The ssm from which the measurement was generated
        Tag = []
    end
    
    methods
        function this = GroundTruthStateX(varargin)
        % GroundTruthStateX constructor
        %
        % Parameters
        % ----------
        % Vector: numeric
        %   The state vector
        % Timestamp: datetime, optional
        %   An optional timestamp that can be attached to a StateX object.
        % Tag: TagX, optional
        %   A tag linking the ground-truth state to a ground-truth track.
            
            [tag, other] = GroundTruthStateX.extract_tag(varargin);
            
            this@StateX(other{:});
            
            this.Tag = tag;
        end
    end
    
    methods (Static)
        function [tag, other] = extract_tag(varargs)
            tag = [];
            other = varargs;
            n = numel(varargs);
            
            if n>1
            % If nargin > 1, then expect an tag at index 2 or 3
                if isa(varargs{2}, 'TagX')
                    tag = varargs{2};
                    other = varargs(setdiff(1:n,2));
                elseif n>2 && isa(varargs{3}, 'TagX')
                    tag = varargs{3};
                    other = varargs(setdiff(1:n,3));
                end
            end
        end
    end
end

