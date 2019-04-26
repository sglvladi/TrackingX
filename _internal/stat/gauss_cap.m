function [x_new,P_new,w_new]= gauss_cap(x,P,w,max_number)
% gauss_cap Restrict the number of components in a Gaussian Mixture. 
%
% TODO:
% ----
% * Generate documentation
%
% Credit: 
% -------
% * This code has been taken from http://ba-tuong.vo-au.com/codes.html.
% * Minor modifications may have been applied to the original version.

if length(w) > max_number
    [notused,idx]= sort(w,2,'descend');
    w_new= w(idx(1:max_number)); w_new = w_new * (sum(w)/sum(w_new));
    x_new= x(:,idx(1:max_number));
    P_new= P(:,:,idx(1:max_number));
else
    x_new = x;
    P_new = P;
    w_new = w;
end
