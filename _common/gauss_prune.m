function [x_new,P_new,w_new]= gauss_prune(x,P,w,elim_threshold)

idx= find( w > elim_threshold );
w_new= w(idx); %w_new = w_new*(sum(w)/sum(w_new));
x_new= x(:,idx);
P_new= P(:,:,idx);