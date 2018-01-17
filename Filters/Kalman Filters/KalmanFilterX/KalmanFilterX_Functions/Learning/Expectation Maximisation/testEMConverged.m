function [converged, delta_loglik, positive] = testEMConverged(loglik_k, loglik_km1, c)
% TESTEMCONVERGED Perform EM convergence check as described by [1]
%
% INPUTS:   loglik_k    The loglikelihood computed for time k
%           loglik_km1  The loglikelihood computed for time k-1
%           c           The relative change cut-off threshold, as denoted by [1]
%
% OUTPUTS:  converged   A boolean stating if EM has converged
%
% [1] A. W. Blocker, An EM algorithm for the estimation of affine state-space systems with or without known inputs, (2008)
%
% January 2018 Lyudmil Vladimirov, University of Liverpool.

    converged = false;
    delta_loglik = loglik_k-loglik_km1;
    if(delta_loglik>=0)
            positive = true;
    else
            positive = false;
    end
    delta_loglik = 2*abs(delta_loglik)/abs(loglik_k+loglik_km1+eps);
    if(delta_loglik<c)
        converged = true;
    end
end