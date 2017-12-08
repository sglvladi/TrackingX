function [z,J] = jaccsd(fun,x)
% JACCSD Compute Jacobian through complex step differentiation
% 
% INPUTS:   x    - A state vector
%           fun  - A (non-linear) transition function
%
% OUTPUTS:  z    - The transformed vector
%           J    - The computed Jacobian
%
% USAGE:    [z,J] = jaccsd(fun,x)         
    z   = fun(x);
    n   = numel(x);
    h   = n*eps;
    J   = imag(fun(repmat(x,1,n)+eye(n)*h*1i))./h;
end