function [x_new,P_new,w_new]= gauss_prune(x,P,w,elim_threshold)
% gauss_prune Prune mixture components whose weights fall below a given
%   threshold
%
% TODO:
% ----
% * Generate documentation
%
% Credit: 
% -------
% * This code has been taken from http://ba-tuong.vo-au.com/codes.html.
% * Minor modifications may have been applied to the original version.

idx= find( w > elim_threshold );
w_new= w(idx);
x_new= x(:,idx);
P_new= P(:,:,idx);