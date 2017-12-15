function [pSmooth] = ParticleFilterX_SmoothBS(pFilt,wFilt,f,f_eval,Ns)
% PARTICLEFILTERX_SMOOTHBS Perform Backward Simulation PF Smoothing (Alg. 1 of [1])
%
% INPUTS:   pFilt       The (xDim x Np x N) matrix of filtered particles states for each time-step
%           wFilt       The (1 x Np x N) matrix of filtered particle weights for each time-step
%           f           The (noiseless) state transition function x_kp1 = f(x_k), 
%                       where x_k and x_kp1 are both (xDim x Np) matrices.
%                       (e.g. f = @(x) sin(x))
%           f_eval      The state transition pdf evaluation function w_k+1 = p(x_k+1|x_k) -> w_k+1 = f_eval(x_k+1,x_k),
%                       where x_kp1,x_k are both (xDim x Np) matrices, while w_kp1 is a (1 x Np) vector.
%                       (e.g. f_eval = @(x_kp1,x_k) mvnpdf(x_kp1,x_k,Q) assuming Gaussian process noise with covariance Q
%           Ns          The number of backward simulations to be performed
%
% OUTPUTS:  pSmooth     The (xDim x Ns x N) matrix of smoothed particle states for each time-step
%
% [1] Godsill et al., Monte Carlo Smoothing for Nonlinear Time Series, (2004).
%
% December 2017 Lyudmil Vladimirov, University of Liverpool.

    % Allocate memory
    [xDim,Np,N] = size(pFilt); 
    pSmooth = cell(1,N); %zeros(xDim,Ns,N);
  
    % Perform backward simulation
    parfor i=1:Ns
        
        pSmooth{i} = zeros(xDim,N);
        
        % Select a particle with probability wFilt
        ind = categ_rnd(wFilt(:,:,N));        
        pSmooth{i}(:,N) = pFilt(:,ind,N);
        
        % Perform smoothing
        for k = N-1:-1:1
            % Get predicted particles
            pPred = f(pFilt(:,:,k));
            
            % Compute weights
            wPred  = wFilt(:,:,k) .* f_eval(pSmooth{i}(:,k+1),pPred);
            if(~sum(wPred)); wPred = wPred + eps; end
            wPred  = wPred ./ sum(wPred);   % Normalise

            % Select a new particle with probability wPred
            ind = categ_rnd(wPred);
            pSmooth{i}(:,k) = pFilt(:,ind,k);
        end
    end
    pSmooth = permute(cat(3,pSmooth{:}),[1 3 2]);
end