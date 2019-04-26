classdef UuidTagGeneratorX < TagGeneratorX
% UuidTagGeneratorX class
%
% Summary of UuidTagGeneratorX:
%
%   This is a class implementation of a Tag Generator that yields (unique) 
%   TagX objects, by generating (unique) identifiers as random samples from 
%   a given population of identifiers, with or without replacement. 

    methods 
        
        function tags = generate(~, num_tags, varargin)
        % generate Generate a (number of) tag(s)
        %
        % Parameters
        % ----------
        % num_tags: scalar, optional
        %   The number of tags to be generated (default = 1)
        %
        % Returns
        % -------
        % TagX
        %   A (list of) generated tag object(s)
            
            num_tags(nargin < 2) = 1;
            tags = TagX.empty(0,num_tags);
            for i = 1:num_tags
                tags(i) = TagX(java.util.UUID.randomUUID);
            end
        end
    end

end

