function [newWeights] = ParticleFilterX_UpdateWeights(lik,y,parts,weights)
% PARTICLEFILTERX_UPDATEWEIGHTS Perform the discrete-time PF weight update
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
% newWeights: row vector
%   The (1 x Np) updated weights matrix matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Update and normalise weights
    newWeights = weights.*lik(y,parts);
    if all(~newWeights)
%         warning("All weights evaluated as 0! Setting equal to eps");
        newWeights = newWeights + eps;
    end
    newWeights = newWeights./sum(newWeights,2);       
end
