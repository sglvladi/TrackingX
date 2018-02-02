classdef MultinomialResamplerX < ResamplerX
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
% [1] E. A. Wan and R. Van Der Merwe, "The unscented Kalman filter for nonlinear estimation," 
%     Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and 
%     Control Symposium (Cat. No.00EX373), Lake Louise, Alta., 2000, pp. 153-158.
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
        
        function [newSamples, newWeights] = resample(this,samples,weights)
        % RESAMPLE Perform multinomial resampling
        %   
        % DESCRIPTION: 
        % * idx = resample(this,weights) returns a row vector idx containing
        %   the indexes of samples generated through by performing multinomial
        %   resampling based on the vector weights.
        %
        % See also SystematicResamplerX/resample
        
            N = numel(weights);
            idx = randsample(1:N, N, true, weights);
            newSamples = samples(:,idx);           
            newWeights = repmat(1/N, 1, N); 
        end
    end
end

