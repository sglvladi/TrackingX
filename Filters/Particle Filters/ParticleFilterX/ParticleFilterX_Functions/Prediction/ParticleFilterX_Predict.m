function [predictedParts] = ParticleFilterX_Predict(f,parts,wk)
% PARTICLEFILTERX_PREDICT Perform the discrete-time PF state prediction step,
% under the assumption of additive process noise.
% 
% NOTES: 
% * Control inputs are NOT currently supported.
%
% INPUTS:   parts  - A (xDim x Np) particle matrix from the previous 
%                    time-step.
%           f      - A (non-linear) state transition function.
%           wk     - A (xDim x Np) process noise matrix.
%
% OUTPUTS:  predictedParts - The (xDim x Np) predicted particle matrix.
%
% October 2017 Lyudmil Vladimirov, University of Liverpool.
    
   predictedParts = f(parts, wk); % Simply propagate all particles        
end