classdef RandSampleTagGeneratorX < TagGeneratorX
% RandSampleTagGeneratorX class
%
% Summary of RandSampleTagGeneratorX:
%
%   This is a class implementation of a Tag Generator that yields (unique) 
%   TagX objects, by generating (unique) identifiers as random samples from 
%   a given population of identifiers, with or without replacement. 

    properties 
        allocatedIDs_ = []
        withReplacement_ = false;
        idSamples_
    end
    
    properties (Dependent)
        % NumAllocatedTags: scalar
        %   The number of tags already allocated.
        NumAllocatedTags
        
        % NumTotalTags: scalar
        %   The maximum number of tags that the generator can allocate.
        NumTotalTags
    end

    methods 
        
        function this = RandSampleTagGeneratorX(ids, with_replacement)
        % RandSampleTagGeneratorX constructor
        %
        % Parameters
        % ----------
        % ids: vector of objects
        %   A range of objects to be used as (unique) identifier values. Any
        %   type of object can be used, including custom classes, as long
        %   as they have appropriately defined (non-)equality operators.
        % with_replacement: logical, optional
        %   Set to 'true' to enable sample replacement. Note that this will
        %   lead to generation of non-unique identifiers. (defalut = false)
            
            this.idSamples_ = ids;
            if nargin > 1
                this.withReplacement_ = with_replacement;
            end
        end
        
        function tags = generate(this, num_tags, varargin)
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
            
            try
                ids = randsample(setdiff(this.idSamples_, this.allocatedIDs_),... 
                                 num_tags,...
                                 this.withReplacement_);
                this.allocatedIDs_ = union(this.allocatedIDs_,ids);
                tags = TagX.empty(0,num_tags);
                for i = 1:num_tags
                    tags(i) = TagX(ids(i));
                end
            catch ME
                if (strcmp(ME.identifier,'stats:randsample:SampleTooLarge'))
                    msg = ['The TagGeneratorX object has exhausted all possible unique tag keys.', ...
                           'You have requested for ', num2str(num_tags), ' new tags ',...
                           'to be generated, while only ', num2str(numel(setdiff(this.UsedValues_,this.Range))),...
                           '. \nConsider increasing the range and/or deteting unused tags'];
                    causeException = MException('trackingx:taggenerator:SampleTooLarge',msg);
                    ME = addCause(ME,causeException);
                end
                rethrow(ME)
            end
            
        end
        
        function success = forget(this, tags)
        % forget Remove (an array of) tag(s) from the list of known allocated tags.
        %
        % Parameters
        % ----------
        % tags: TagX object array
        %   A (list of) tag(s) to be forgotten.
        %
        % Returns
        % -------
        % logical
        %   A (list of) indicator(s), returning True if the tag was 
        %   found and removed, and False otherwise
        
            n = numel(tags);
            success = false(1,n);
            i = 1;
            for tag = tags
               if ismember(tag.ID, this.UsedValues_)
                   this.UsedValues_ = setdiff(this.UsedValues_,tag.ID);
                   success(i) = 1;
               end
               i=i+1;
            end
            
        end
    end

end

