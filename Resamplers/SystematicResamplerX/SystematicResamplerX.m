classdef SystematicResamplerX < ResamplerX & matlabshared.tracking.internal.SystematicResampler
% SYSTEMATICRESAMPLERX Class 
%
% Summary of SystematicResamplerX:
% This is a class implementation of a systematic resampler.
%
% SystematicResamplerX Properties:
%   None
%
% SystematicResamplerX Methods:
%    SystematicResamplerX  - Constructor method
%    resample - Perform systematic resampling 
% 
% See also MultinomialResamplerX

    properties
    end
    
    methods
        function this = SystematicResamplerX(varargin)
        % SYSTEMATICRESAMPLERX Constructor method
        %   
        % DESCRIPTION: 
        % * SystematicResamplerX() returns a "SystematicResamplerX" object 
        %   handle
        %
        % See also SystematicResamplerX/resample
            
        end
        
        function [newSamples, newWeights, idx] = resample(this, samples, weights, Nnew)
        % RESAMPLE Perform systematic resampling
        %   
        % Parameters
        % ----------
        % samples: (NumDims x NumSamples) matrix
        %   A matrix, whose columns correspond to individual samples
        % weights: (1 x NumSamples) row vector
        %   A row vector, whose columns correspond to the weights of the
        %   respective columns/samples of the samples matrix
        % Nnew: scalar, optional
        %   The number of samples to be resampled.
        %   (default = NumSamples, which implies that the same ammount of
        %   samples will be generated, as the number of samples provided)
        %
        % Returns
        % -------
        % newSamples: (NumDims x Nnew) matrix
        %   A matrix, whose columns correspond to resampled samples.
        % newWeights: (1 x Nnew) row vector
        %   A row vector, whose columns correspond to the weights of the
        %   respective columns/samples of the newSamples matrix
        % idx: (1 x Nnew) row vector
        %   A row vector of indices, mapping each sample in newSamples to the
        %   corresponding sample in samples.
        %
        % Usage
        % -----
        % * [newSamples, newWeights, idx] = resample(this,weights) returns a 
        %   row vector idx containing the indexes of samples generated through 
        %   by performing systematic resampling based on the vector weights.
        %
        % See also MultinomialResamplerX/resample
            
            N = numel(weights);
            if(nargin<4)
                Nnew = N;
            end
            
            idx = resample@matlabshared.tracking.internal.SystematicResampler(this,weights,Nnew);
            newSamples = samples(:,idx);           
            newWeights = repmat(1/Nnew, 1, Nnew); 
        end


    end
end

