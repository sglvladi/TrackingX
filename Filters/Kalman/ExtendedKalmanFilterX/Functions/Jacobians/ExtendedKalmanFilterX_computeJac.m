function [xT,Jac] = ExtendedKalmanFilterX_computeJac(fun,x)
% EXTENDEDKALMANFILTERX_COMPUTEJAC Compute Jacobian through complex step 
% differentiation
% 
% Parameters
% ----------
% x: column vector
%   A state vector
% fun: function handle
%   A (non-linear) transition function
%   Must be of the form "fun = @(x)...." 
%
% Returns
% -------
% xT: column vector
%   The transformed vector
% Jac: matrix
%   The computed Jacobian
    
    xT   = fun(x);
    Jac     = jacobian(fun,x);
%     xDim    = numel(x);
%     h       = xDim*eps;
%     Jac     = imag(fun(x(:,ones(1,xDim))+eye(xDim)*h*1i))./h;
end