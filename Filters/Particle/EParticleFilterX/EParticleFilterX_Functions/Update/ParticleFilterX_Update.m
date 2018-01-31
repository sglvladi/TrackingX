function [newParts,newWeights,x] = ParticleFilterX_Update(lik,y,parts,weights,resampling_strategy)
% PARTICLEFILTERX_UPDATE Perform the discrete-time PF weight update
% step, under the assumption of additive process noise.
%
% INPUTS:   lik(y,x) - Measurement likelihood function handle.
%           y        - A (yDim x Nm) measurement matrix
%           parts    - A (xDim x Np) particle matrix.
%           weights  - A (1 x Np) weights matrix.
%
% OUTPUTS:  newParts   - The (xDim x Np) resampled particle matrix.
%           newWeights - The (1 x Np) updated weights matrix matrix.
%           x          - The (xDim x 1) state mean estimate vector
%                        (Optional, only computed if requested)
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
   % Update particle weights
    newWeights = ParticleFilterX_UpdateWeights(lik,y,parts,weights);

    % Resample
    [newParts, newWeights] = ParticleFilterX_Resample(parts, newWeights, resampling_strategy);
    
    if nargout==3
        % Compute estimated state
        x = sum(newWeights.*newParts,2);      
    end
end
