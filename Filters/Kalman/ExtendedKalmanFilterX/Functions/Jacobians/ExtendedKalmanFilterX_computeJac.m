function [xT,Jac] = ExtendedKalmanFilterX_computeJac(fun,x)
% EXTENDEDKALMANFILTERX_COMPUTEJAC Compute Jacobian through complex step 
% differentiation
% 
% INPUTS:   x    - A state vector
%           fun  - A (non-linear) transition function
%                  Must be of the form "fun = @(x)...." 
%
% OUTPUTS:  xT   - The transformed vector
%           Jac  - The computed Jacobian
    
    xT   = fun(x);
    xDim    = numel(x);
    h       = xDim*eps;
    Jac     = imag(fun(x(:,ones(1,xDim))+eye(xDim)*h*1i))./h;
end