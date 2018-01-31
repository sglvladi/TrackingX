function [predictedParts] = EParticleFilterX_Predict(parts,wk,f,Q,h,R,u,b,Qu)
% EPARTICLEFILTERX_PREDICT Perform the discrete-time EPF state prediction step,
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

    % Compute EKF prior mean and covariance
    this.ekf.Params.x = sum(repmat(this.Params.w,size(this.Params.particles,1),1).*this.Params.particles,2);
    this.ekf.Params.P = weightedcov(this.Params.particles',this.Params.w');

    % Iterate EKF to obtain Optimal Proposal
    this.ekf.Predict();
    this.Params.xPred = this.ekf.Params.xPred;
    this.Params.PPred = this.ekf.Params.PPred;        
end