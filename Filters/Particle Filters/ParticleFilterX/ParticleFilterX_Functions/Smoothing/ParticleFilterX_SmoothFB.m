function [wSmooth] = ParticleFilterX_SmoothFB(pFilt,wFilt,f_eval)
% PARTICLEFILTERX_SMOOTHBS Perform Backward recursion of Forward-Backward PF Smoother 
%  (As described in Sec. 3.1 of [1])
%
% INPUTS:   pFilt       The (xDim x Np x N) matrix of filtered particles states for each time-step
%           wFilt       The (1 x Np x N) matrix of filtered particle weights for each time-step
%           f_eval      The state transition pdf evaluation function w_k+1 = p(x_k+1|x_k) -> w_k+1 = f_eval(x_k+1,x_k),
%                       where x_kp1,x_k can be (xDim x Np) matrices, producing the (Np x Np) w_kp1 matrix.
%                       (e.g. f_eval = @(x_kp1,x_k) mvnpdf(x_kp1,x_k,Q) assuming Gaussian process noise with covariance Q
%
% OUTPUTS:  wSmooth     The (1 x Np x N) matrix of smoothed particle weights for each time-step
%
% [1] M. Klaas, M. Briers, N. de Freitas, A. Doucet, S. Maskell, and D. Lang. 2006. Fast particle smoothing: if I had a million particles.
%
% December 2017 Lyudmil Vladimirov, University of Liverpool.

    % Allocate memory
    [xDim,Np,N] = size(pFilt); 
    wSmooth = zeros(1,Np,N);
    wSmooth(1,:,N) = wFilt(1,:,N);
    % Perform Backward Recursion
    for k = N-1:-1:1
        disp(k);
        lik = f_eval(pFilt(:,:,k+1), pFilt(:,:,k));
        denom = sum(wFilt(ones(Np,1),:,k).*lik,2)'; % denom(1,j)
        wSmooth(1,:,k) = wFilt(1,:,k) .* sum(wSmooth(ones(Np,1),:,k+1).*lik'./denom(ones(Np,1),:),2)';   
    end
    
    % The above is equivalent to (but much faster than):
    %  for k = N-1:-1:1
    %     disp(k);
    %     denom = zeros(1,this.Params.Np);
    %     lik = zeros(Np,Np);
    %     for j = 1:this.Params.Np
    %         denom_j = 0;
    %         for t = 1:this.Params.Np
    %             lik(j,t) = f_eval(pFilt(:,j,k+1), pFilt(:,t,k));
    %             denom_j = denom_j + wFilt(1,t,k)*lik(j,t);
    %         end
    %         denom(j) = denom_j;
    %     end
    %     for i = 1:Np
    %         sum_j = 0;
    %         for j = 1:Np
    %             enum = wSmooth(1,j,k+1)*lik(j,i);
    %             sum_j = sum_j + enum/denom(j);
    %         end
    %         wSmooth(1,i,k) = wFilt(1,i,k) * sum_j; 
    %     end
    % end
end