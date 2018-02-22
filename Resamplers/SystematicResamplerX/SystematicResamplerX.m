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
%    resample - Perform multinomial resampling 
%
% [1] E. A. Wan and R. Van Der Merwe, "The unscented Kalman filter for nonlinear estimation," 
%     Proceedings of the IEEE 2000 Adaptive Systems for Signal Processing, Communications, and 
%     Control Symposium (Cat. No.00EX373), Lake Louise, Alta., 2000, pp. 153-158.
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
        % DESCRIPTION: 
        % * idx = resample(this,weights) returns a row vector idx containing
        %   the indexes of samples generated through by performing systematic
        %   resampling based on the vector weights.
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

