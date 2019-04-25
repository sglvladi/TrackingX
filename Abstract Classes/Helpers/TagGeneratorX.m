classdef (Abstract) TagGeneratorX < BaseX
% TagGeneratorX abstract class
%
% Summary of TagGeneratorX:
%   This is the base class for all Tag Generators

    methods (Abstract)
        generate(this); % Generate new tag(s)
%         forget(this);   % Remove tag(s) from list of allocated tags
    end

end

