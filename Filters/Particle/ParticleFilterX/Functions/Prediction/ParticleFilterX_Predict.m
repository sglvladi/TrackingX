function [predictedParts] = ParticleFilterX_Predict(f,parts,wk)
% PARTICLEFILTERX_PREDICT Perform the discrete-time PF state prediction step,
% under the assumption of additive process noise.
% 
% NOTES: 
% * Control inputs are NOT currently supported.
%
% Parameters
% ----------
% parts: matrix
%   A (xDim x Np) particle matrix from the previous 
%   time-step, where xDim is the number of state dimentions and Np is the
%   number of particles.
% f: function handle
%   A (non-linear) state transition function of the form f(x).
% wk: matrix
%   A (xDim x Np) process noise matrix.
%
% Returns
% -------
% predictedParts: matrix
%   The (xDim x Np) predicted particle matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
   predictedParts = f(parts, wk); % Simply propagate all particles        
end