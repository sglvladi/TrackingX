function x = anglewrap(x, angleidx)
% ANGLEWRAP Ensure that angle indices are within (-pi,pi)
%
% Usage
% =====
% x = anglewrap(x, angleidx) correct the indices angleidx of x.
%
% Author: Paul Horridge

x(angleidx,:) = mod(x(angleidx,:) + pi, 2*pi) - pi;
