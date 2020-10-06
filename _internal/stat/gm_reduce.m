function [x,P] = gm_reduce(x_i, P_i, w_i)
%GM_REDUCE Reduce a Gaussian Mixture to a single Gaussian
%   [x,P] = GM_REDUCE(x_i, P_i, w_i) reduces the Gaussian Mixture described
%   by the set of means x_i, coriances P_i and weights w_i to a single
%   Gaussian with mean x and covariace P. The means x_i are expected to have
%   size(nx, nm), while the covariance matrices P_i size(nx, nx, nm). The
%   mixture weights w_i shoud be provided as a vector of size(nm,1) or 
%   size(1,nm). Reduction performed according to equations (18)-(20) of [1].
%   
%   [1] D. F. Crouse, P. Willett, K. Pattipati and L. Svensson, "A look at 
%       Gaussian mixture reduction algorithms," 14th International 
%       Conference on Information Fusion, Chicago, IL, 2011, pp. 1-8.
%   
%   Author: Lyudmil Vladimirov

    [nx, nm] = size(x_i);   % Compute the state size and number of models
    w_i = w_i(:)/sum(w_i);  % Ensure w_i is a normalised column vector
    x = x_i * w_i;          % Equivalent to (19) of [1]
    v = x_i - x;          % Mixed mean differences of (11.6.6-10)
    P = zeros(nx, nx); 
    % Compute sum of terms within brackets in (11.6.6-10)
    for i = 1:nm
        P = P + w_i(i)*(P_i(:,:,i) + v(:,i)*v(:,i)');     
    end
end