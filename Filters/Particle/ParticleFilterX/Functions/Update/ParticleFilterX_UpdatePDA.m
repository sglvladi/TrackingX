function [newWeights] = ParticleFilterX_UpdatePDA(lik,y,parts,weights,W,LikelihoodMatrix)
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
% W: row vector
%   The (1 x nY+1) measurement association/mixture weights 
%   vector. (dummy measurement assumed at index 1)
% LikelihoodMatrix: matrix
%   A (Nm x Np) measurement likelihood matrix
%   per particle. (Optional, If not provided then it gets computed internally)
%
% Returns
% -------
% newParts: matrix
%   The (xDim x Np) resampled particle matrix.
% newWeights: row vector
%   The (1 x Np) updated weights matrix matrix.
% x : column vector, optional
%   The (xDim x 1) state mean estimate vector
%  (Optional, only computed if requested)
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
    if(nargin<6)
        LikelihoodMatrix = lik(y, parts);
    end
    
    numMeasurements = size(y,2);
    
    if(numMeasurements==0)
        newWeights = weights;
        return;
    else
        Ck = sum(LikelihoodMatrix.*weights,2)+eps;
        weightsPerMeasurement = LikelihoodMatrix./Ck.*weights;
        newWeights = W(1)*weights + sum(W(2:end)'.*weightsPerMeasurement,1);

        % Normalize weight vector
        newWeights = newWeights./sum(newWeights,2);
    end
end
