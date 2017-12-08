function [newParts,newWeights] = ParticleFilterX_Resample(parts,weights,strategy)
% PARTICLEFILTERX_RESAMPLE Perform particle resampling
%
% INPUTS:   parts    - A (xDim x Np) particle matrix.
%           weights  - A (1 x Np) weights matrix (normalised).
%           strategy - A char vector indicating the desired resampling
%                      strategy. Options are:
%                       1) 'multinomial_resampling'
%                       2) 'systematic_resampling' 
%
% OUTPUTS:  newParts   - The (xDim x Np) resampled particle matrix.
%           newWeights - The (1 x Np) uniform weights matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    Np = length(weights);  % Np = number of particles
    switch strategy
       case 'multinomial_resampling'
          with_replacement = true;
          idx = randsample(1:Np, Np, with_replacement, weights);
       case 'systematic_resampling'
          % this is performing latin hypercube sampling on wk
          edges = min([0 cumsum(weights)],1); % protect against accumulated round-off
          edges(end) = 1;                 % get the upper edge exact
          u1 = rand/Np;
          % this works like the inverse of the empirical distribution and returns
          % the interval where the sample is to be found
          [~, ~, idx] = histcounts(u1:1/Np:1, edges);
       otherwise
          error('Resampling strategy not implemented\n')
    end
    newParts = parts(:,idx);           % extract new particles
    newWeights = repmat(1/Np, 1, Np);  % uniform weights      
end
