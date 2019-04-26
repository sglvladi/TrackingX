function [r,x,P,lambdau,xu,Pu] = predict(r,x,P,lambdau,xu,Pu,model)
%PREDICT: PREDICT MULTI-BERNOULLI AND POISSON COMPONENTS
%Syntax: [r,x,P] = tomb(pupd,rupd,xupd,Pupd,pnew,rnew,xnew,Pnew)
%Input:
% r(i), x(:,i) and P(:,:,i) give the probability of existence, state 
%  estimate and covariance for the i-th multi-Bernoulli component (track)
% lambdau(k), xu(:,k) and Pu(:,:,k) give the intensity, state estimate and 
%  covariance for the k-th mixture component of the unknown target Poisson
%  Point Process (PPP)
%model is structure describing target birth and transition models
%Output:
% Predicted components in same format as input
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

% Get multi-Bernoulli prediction parameters from model
F = model.F;
Q = model.Q;
Ps = model.Ps;

% Get birth parameters from model
lambdab = model.lambdab;
nb = length(lambdab);
xb = model.xb;
Pb = model.Pb;
lambdab_threshold = 1e-4;

% Interpret length of inputs
n = length(r);
nu = length(lambdau);

% Implement prediction algorithm

% Predict existing tracks
for (i = 1:n),
  r(i) = Ps*r(i);
  x(:,i) = F*x(:,i);
  P(:,:,i) = F*P(:,:,i)*F' + Q;
end;

% Predict existing PPP intensity
for (k = 1:nu),
  lambdau(k) = Ps*lambdau(k);
  xu(:,k) = F*xu(:,k);
  Pu(:,:,k) = F*Pu(:,:,k)*F' + Q;
end;

% Incorporate birth intensity into PPP

% Allocate memory
lambdau(end+nb) = 0;
xu(:,end+nb) = 0;
Pu(:,:,end+nb) = 0;
for (k = 1:nb),
  lambdau(nu+k) = lambdab(k);
  xu(:,nu+k) = xb(:,k);
  Pu(:,:,nu+k) = Pb(:,:,k);
end;

% Not shown in paper--truncate low weight components
ss = lambdau > lambdab_threshold;
lambdau = lambdau(ss);
xu = xu(:,ss);
Pu = Pu(:,:,ss);