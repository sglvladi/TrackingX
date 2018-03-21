classdef (Abstract) BaseX < matlab.mixin.Copyable & matlab.mixin.SetGetExactNames 
% BaseX class
%
% SUMMARY:
% Base class for all components of the TrackingX framework. 
% 
%  February 2018 Lyudmil Vladimirov, University of Liverpool
    
    properties
    end
    
    methods (Access = protected)
        function cp = copyElement(this)
        % COPYELEMENT Override matlab.mixin.Copyable copyElement function,
        %  such as to perform deep copy of all handle properties.
        
            % Shallow copy object
            cp = copyElement@matlab.mixin.Copyable(this);
            
            % Proceed to recursively deep copy all handles
            props = properties(this); 
            for i = 1: length(props) 
                prop = findprop(this,props{i});
                if(~prop.Dependent&&~strcmp(prop.SetAccess,'immutable'))
                    tmp = this.(props{i});
                    if(isa(tmp,'handle')) 
                        cp.(props{i}) = copy(tmp); 
                    else
                        cp.(props{i}) = tmp ; 
                    end 
                end
            end  
        end
    end

end