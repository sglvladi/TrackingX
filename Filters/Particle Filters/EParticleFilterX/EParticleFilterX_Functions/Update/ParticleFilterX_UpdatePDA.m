function [newParts,newWeights,x] = ParticleFilterX_UpdatePDA(lik,y,parts,weights,resampling_strategy,W,LikelihoodMatrix)
% PARTICLEFILTERX_UPDATE Perform the discrete-time PF weight update
% step, under the assumption of additive process noise.
%
% INPUTS:   lik(y,x) - Measurement likelihood function handle.
%           y        - A (yDim x Nm) measurement matrix
%           parts    - A (xDim x Np) particle matrix.
%           weights  - A (1 x Np) weights matrix.
%           resampling_strategy - A char vector indicating the desired 
%                                 resampling strategy. Options are:
%                                   1) 'multinomial_resampling'
%                                   2) 'systematic_resampling'
%           W        - The (1 x nY+1) measurement association/mixture weights 
%                      vector. (dummy measurement assumed at index 1)
%           LikelihoodMatrix - A (Nm x Np) measurement likelihood matrix
%                              per particle.
%                              (Optional, If not provided then it gets
%                               computed internally)
%
% OUTPUTS:  newParts   - The (xDim x Np) resampled particle matrix.
%           newWeights - The (1 x Np) updated weights matrix matrix.
%           x          - The (xDim x 1) state mean estimate vector
%                        (Optional, only computed if requested)
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    if(nargin<7)
        LikelihoodMatrix = lik(y, parts);
    end
    
    Np = length(weights);  % Np = number of particles
    
    % Compute expected likelihoods
    ExpectedLikelihoodMatrix = LikelihoodMatrix.*repmat(W(2:end)',1,Np);
    expectedLikelihoods = sum(ExpectedLikelihoodMatrix,1);
    el = sum(expectedLikelihoods,2)/sum(W(2:end)) - sum(expectedLikelihoods,2);
    ExpectedLikelihoodMatrix = [repmat(el/Np,1,Np); ExpectedLikelihoodMatrix];
    expectedLikelihoods = sum(ExpectedLikelihoodMatrix,1);

    % Update particle weights
    newWeights = weights .* expectedLikelihoods;

    % Normalize weight vector
    newWeights = newWeights./sum(newWeights,2);

    % Resampling
    [newParts, newWeights] = ParticleFilterX_Resample(parts, newWeights, resampling_strategy);
    
    if nargout==3
        % Compute estimated state
        x = sum(newWeights.*newParts,2);      
    end
end
