function [newWeights] = ParticleFilterX_UpdateWeights(lik,y,parts,weights)
% PARTICLEFILTERX_UPDATEWEIGHTS Perform the discrete-time PF weight update
% step, under the assumption of additive process noise.
%
% INPUTS:   lik(y,x) - A measurement likelihood function.
%           y        - A (yDim x Nm) measurement matrix
%           parts    - A (xDim x Np) particle matrix.
%           weights  - A (1 x Np) weights matrix.
%
% OUTPUTS:  newWeights - The (1 x Np) updated weights matrix matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Update and normalise weights
    newWeights = weights.*lik(y,parts);
    newWeights = newWeights./sum(newWeights,2);       
end
