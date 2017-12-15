function [pSmooth] = ParticleFilterX_SmoothBS(pFilt,wFilt,f_eval,Ns,parallel)
% PARTICLEFILTERX_SMOOTHBS Perform Backward Simulation PF Smoothing (Alg. 1 of [1])
% 
% NOTE: Backward simulations are performed using a "parfor" loop to increase
% performance. If this is causing an issue, simply set "parallel=false".
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
%           parallel    Set (true|false) to (enable|disable) the usage of "parfor" for backward simulations
%                       (Default = true)
%
% OUTPUTS:  pSmooth     The (xDim x Ns x N) matrix of smoothed particle states for each time-step
%
% [1] Godsill et al., Monte Carlo Smoothing for Nonlinear Time Series, (2004).
%
% December 2017 Lyudmil Vladimirov, University of Liverpool.
    
    % Validate inputs 
    if(nargin<5); parallel = true; end
    
    % Allocate memory
    [xDim,Np,N] = size(pFilt); 
    pSmooth = cell(1,N); 
  
    % Perform backward simulation
    if(parallel)
        % Use "parfor"
        parfor i=1:Ns
            pSmooth{i} = zeros(xDim,N);

            % (Step 1) Select a particle with probability wFilt
            ind = categ_rnd(wFilt(:,:,N));        
            pSmooth{i}(:,N) = pFilt(:,ind,N);

            % (Step 2) Perform smoothing
            for k = N-1:-1:1

                % Compute backward predicted weights
                wPred  = wFilt(:,:,k) .* f_eval(pSmooth{i}(:,k+1),pFilt(:,:,k));
                % IMPORTANT!: Protect from failure due to precision errors
                if(~sum(wPred)); wPred = wPred + eps; end
                % Normalise weights
                wPred  = wPred ./ sum(wPred);

                % Select a new particle with probability wPred
                ind = categ_rnd(wPred);
                pSmooth{i}(:,k) = pFilt(:,ind,k);
            end
        end
    else
        % Do NOT use "parfor"
        for i=1:Ns
            pSmooth{i} = zeros(xDim,N);

            % (Step 1) Select a particle with probability wFilt
            ind = categ_rnd(wFilt(:,:,N));        
            pSmooth{i}(:,N) = pFilt(:,ind,N);

            % (Step 2) Perform smoothing
            for k = N-1:-1:1

                % Compute backward predicted weights
                wPred  = wFilt(:,:,k) .* f_eval(pSmooth{i}(:,k+1),pFilt(:,:,k));
                % IMPORTANT!: Protect from failure due to precision errors
                if(~sum(wPred)); wPred = wPred + eps; end
                % Normalise weights
                wPred  = wPred ./ sum(wPred);

                % Select a new particle with probability wPred
                ind = categ_rnd(wPred);
                pSmooth{i}(:,k) = pFilt(:,ind,k);
            end
        end
    end
    
    % (Step 3) 
    pSmooth = permute(cat(3,pSmooth{:}),[1 3 2]);
end