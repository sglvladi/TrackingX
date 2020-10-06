function [x_0j,P_0j] = imm_merge(x_i, P_i, m_ij)
%IMM_MERGE Perform the IMM merging/mixing step
%   [x_0j, P_0j] = IMM_merge(x_i, P_i, m_ij) computes and returns the 
%   merged means x_0j and covariance matrices P_0j, using the provided
%   means x_i, covariance matrices P_i and mixing weights m_ij. The returned
%   matrices x_0j and P_0j have the same size as the provided inputs x_i
%   and P_i, respectively, where the means x_i are expected to have
%   size(nx, nm), while the covariance matrices P_i size(nx, nx, nm). The
%   mixing weights matrix m_ij shoud be square with size(nm, nm). Merging
%   is performed according to equations (11.6.6-9) and (11.6.6-10) of [1].
%   
%   [1] Yaakov Bar-Shalom, Thiagalingam Kirubarajan, and X.-Rong Li. 2002. 
%       Estimation with Applications to Tracking and Navigation. 
%       John Wiley & Sons, Inc., New York, NY, USA.
%   
%   Author: Lyudmil Vladimirov

    [nx, nm] = size(x_i);    % Compute the state size and number of models
    x_0j = x_i * m_ij;       % Equivalent to (11.6.6-9) of [1]
    v = x_i - x_0j;          % Mixed mean differences of (11.6.6-10)
    P_t = zeros(nx, nx, nm); 
    % Compute sum of terms within brackets in (11.6.6-10)
    for i = 1:nm
        P_t(:,:,i) = P_i(:,:,i)+ v(:,i)*v(:,i)';     
    end
    P_0 = reshape(P_t,[nx^2 nm]); % This is equivalent to iterating through
    P_0 = P_0 * m_ij;             % all i's, multiplying by the mixing 
    P_0j = reshape(P_0,nx,nx,nm); % weights and performing the summation
end