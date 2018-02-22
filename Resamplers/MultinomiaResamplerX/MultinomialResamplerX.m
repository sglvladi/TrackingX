classdef MultinomialResamplerX < ResamplerX & matlabshared.tracking.internal.MultinomialResampler
% MULTINOMIALRESAMPLERX Class 
%
% Summary of MultinomialResamplerX:
% This is a class implementation of a multinomial resampler.
%
% MultinomialResamplerX Properties:
%   None
%
% MultinomialResamplerX Methods:
%    MultinomialResamplerX  - Constructor method
%    resample - Perform multinomial resampling 
% 
% See also SystematicResamplerX
    properties
    end
    
    methods
        function this = MultinomialResamplerX(varargin)
        % MULTINOMIALRESAMPLERX Constructor method
        %   
        % DESCRIPTION: 
        % * MultinomialResamplerX() returns a "MultinomialResamplerX" object 
        %   handle
        %
        % See also MultinomialResamplerX/resample
            
        end
        
        function [newSamples, newWeights,idx] = resample(this,samples,weights,Nnew)
        % RESAMPLE Perform multinomial resampling
        %   
        % DESCRIPTION: 
        % * idx = resample(this,weights) returns a row vector idx containing
        %   the indexes of samples generated through by performing multinomial
        %   resampling based on the vector weights.
        %
        % See also SystematicResamplerX/resample
        
            N = numel(weights);
            if(nargin<4)
                Nnew = N;
            end
            % idx = randsample(1:N, Nnew, true, weights);
            idx = resample@matlabshared.tracking.internal.MultinomialResampler(this,weights,Nnew);
            newSamples = samples(:,idx);           
            newWeights = repmat(1/Nnew, 1, Nnew); 
        end
    end
end

