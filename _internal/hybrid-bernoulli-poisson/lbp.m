function [pupd,pnew] = lbp(wupd,wnew)
%LBP: LOOPY BELIEF PROPAGATION APPROXIMATION OF MARGINAL ASSOCIATION PROBABILITIES
%Syntax: [pupd,pnew] = lbp(wupd,wnew)
%Input:
% wupd(i,j+1) is PDA likelihood for track i/measurment j 
%  e.g., P_d N(z_j;H*mu_i,H*P_i*H'+R)
% wupd(i,1) is miss likelihood for target i
%  e.g., (1-P_d)
% wnew(j) is false alarm/new target intensity for measurement j
%  e.g., lambda_fa(z_j) + lambda^u(z_j)
%Output:
% Estimates of marginal association probabilities in similar format.
%Software implements pseudocode described in http://arxiv.org/abs/1203.2995 
% (accepted for publication, IEEE Transactions on Aerospace and Electronic 
% Systems)
%Copyright 2012 Defence Science and Technology Organisation
%Note:
% As this software is provided free charge and to the full extent permitted
% by applicable law the software is provided "as is" and DSTO and any
% copyright holder in material included in the software make no
% representation and gives no warranty of any kind, either expressed or
% implied, including, but not limited to, implied warranties of
% merchantability, fitness for a particular purpose, non-infringement, or
% the presence or absence of defects or errors, whether discoverable or
% not. The entire risk as to the quality and performance of the software is
% with you.  
%Note: this implementation is matched directly to the pseudocode in the
% journal paper, structured for readability. Large speed improvements can
% (and have) be made by vectorising, clustering, checking for convergence
% less often, etc.

[n,mp1] = size(wupd);
m = mp1-1;
% wupd dimensions n x m+1
% wnew dimensions m x 1
% pupd, pnew dimensions same

eps_conv_threshold = 1e-4;

mu = ones(n,m); % mu_ba
mu_old = zeros(n,m);
nu = zeros(n,m); % mu_ab

pupd = zeros(n,mp1);
pnew = zeros(m,1);

% Run LBP iteration
while (max(abs(mu(:)-mu_old(:))) > eps_conv_threshold),
  mu_old = mu;
  
  for (i = 1:n),
    prd = wupd(i,2:end).*mu(i,:);
    s = wupd(i,1) + sum(prd);
    nu(i,:) = wupd(i,2:end) ./ (s - prd);
  end;
  
  for (j = 1:m),
    s = wnew(j) + sum(nu(:,j));
    mu(:,j) = 1./(s - nu(:,j));
  end;
end;

% Calculate outputs--for existing tracks then for new tracks
for (i = 1:n),
  s = wupd(i,1) + sum(wupd(i,2:end).*mu(i,:));
  pupd(i,1) = wupd(i,1)/s;
  pupd(i,2:end) = wupd(i,2:end).*mu(i,:)/s;
end;

for (j = 1:m),
  s = wnew(j) + sum(nu(:,j));
  pnew(j) = wnew(j)/s;
end;