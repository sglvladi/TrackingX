function [newWeights,x] = ParticleFilterX_Update(lik,y,parts,weights)
% PARTICLEFILTERX_UPDATE Perform the discrete-time PF weight update
% step, under the assumption of additive process noise.
%
% Parameters
% ----------
% lik: function handle 
%   Measurement likelihood function handle of the form lik(y,parts).
% y: matrix
%   A (yDim x Nm) measurement matrix
% parts: matrix
%   A (xDim x Np) particle matrix.
% weights: row vector
%   A (1 x Np) particle weights vector.
%
% Returns
% -------
% newParts: matrix
%   The (xDim x Np) resampled particle matrix.
% newWeights: row vector
%   The (1 x Np) updated weights matrix matrix.
% x: column vector, optional
%   The (xDim x 1) state mean estimate vector
%   (Optional, only computed if requested)
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
   % Update particle weights
    newWeights = ParticleFilterX_UpdateWeights(lik,y,parts,weights);
    
    if nargout==2
        % Compute estimated state
        x = sum(newWeights.*newParts,2);      
    end
end
