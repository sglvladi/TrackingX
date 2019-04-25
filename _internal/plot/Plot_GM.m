function Plot_GM(W,M,V,flag,X)
% Plot_GM(W,M,V,X)
%
% Plots a 1D or 2D Gaussian mixture. 
% 
% Inputs:
%   W(1,k)   - weights of each GM component, k is the number of component
%   M(d,k)   - mean of each GM component, d is the dimension of the random variable 
%   V(d,d,k) - covariance of each GM component
%   flag(1,1)- optional for 2D only,
%               'e' for standard deviation ellipse plot only (default),
%               'l' for landscape plot only using function MESHC,
%               'b' for both ellipse and landscape plots 
%   X(n,d)   - optional for 2D only, data used to learn the GM 
% 
% Written by
%   Patrick P. C. Tsui,
%   PAMI research group
%   Department of Electrical and Computer Engineering
%   University of Waterloo, 
%   October, 2005
%
if nargin < 3,
    disp('W, M and V must be provided!');
    return
elseif nargin == 3,
    [WMVerr,d] = VerifyWMV(W,M,V);
    flag = 'e'; X = [];
    if WMVerr,
        return
    end
elseif nargin == 4,
    [WMVerr,d] = VerifyWMV(W,M,V);
    Ferr = VerifyF(flag);
    X = [];
    if WMVerr | Ferr,
        return
    end
elseif nargin == 5,
    [WMVerr,d] = VerifyWMV(W,M,V);
    Ferr = VerifyF(flag);
    Xerr = VerifyX(X,d);
    if WMVerr | Ferr | Xerr,
        return
    end
elseif nargin > 5,
    disp('Too many inputs!');
    return
end
[d,k] = size(M);
S = zeros(d,k);
R1 = zeros(d,k);
R2 = zeros(d,k);
for i=1:k,  % Determine plot range as 4 x standard deviations
    S(:,i) = sqrt(diag(V(:,:,i)));
    R1(:,i) = M(:,i)-4*S(:,i);
    R2(:,i) = M(:,i)+4*S(:,i);
end
Rmin = min(min(R1));
Rmax = max(max(R2));
R = [Rmin:0.001*(Rmax-Rmin):Rmax];

if d==1,
    clf, hold on
    Q = zeros(size(R));
    for i=1:k,
        P = W(i)*normpdf(R,M(:,i),sqrt(V(:,:,i)));
        Q = Q + P;
        plot(R,P,'r-'); grid on,
    end
    plot(R,Q,'k-');
    xlabel('X');
    ylabel('Probability density');
    title('Gaussian Mixture estimated by EM');
else % d==2
    if flag=='e',
        Plot_Std_Ellipse(M,V,R,X);
    elseif flag=='l',    
        Plot_GM_Landscape(W,M,V,R);
    else
        Plot_Std_Ellipse(M,V,R,X);        
        Plot_GM_Landscape(W,M,V,R);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Plot_GM %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_Std_Ellipse(M,V,R,X)
clf, hold on
if ~isempty(X),
    plot(X(:,1),X(:,2),'r.');
end
[d,k] = size(M);
for i=1:k,
    if V(:,:,i)==zeros(d,d),
        V(:,:,i) = eyes(d,d)*eps;
    end
    [Ev,D] = eig(V(:,:,i));
    iV = inv(V(:,:,i));
    % Find the larger projection
    P = [1,0;0,0];  % X-axis projection operator
    P1 = P * 2*sqrt(D(1,1)) * Ev(:,1);
    P2 = P * 2*sqrt(D(2,2)) * Ev(:,2);
    if abs(P1(1)) >= abs(P2(1)),
        Plen = P1(1);
    else
        Plen = P2(1);
    end
    count = 1;
    step = 0.001*Plen;
    Contour1 = zeros(2001,2);
    Contour2 = zeros(2001,2);
    for x = -Plen:step:Plen,
        a = iV(2,2);
        b = x * (iV(1,2)+iV(2,1));
        c = (x^2) * iV(1,1) - 1;
        Root1 = (-b + sqrt(b^2 - 4*a*c))/(2*a);
        Root2 = (-b - sqrt(b^2 - 4*a*c))/(2*a);
        if isreal(Root1),
            Contour1(count,:) = [x,Root1] + M(:,i)';
            Contour2(count,:) = [x,Root2] + M(:,i)';
            count = count + 1;
        end
    end
    Contour1 = Contour1(1:count-1,:);
    Contour2 = [Contour1(1,:);Contour2(1:count-1,:);Contour1(count-1,:)];
    plot(M(1,i),M(2,i),'k+');
    plot(Contour1(:,1),Contour1(:,2),'k-');
    plot(Contour2(:,1),Contour2(:,2),'k-');
end
xlabel('1^{st} dimension');
ylabel('2^{nd} dimension');
title('Gaussian Mixture estimated by EM');
Rmin = R(1);
Rmax = R(length(R));
axis([Rmin Rmax Rmin Rmax])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Plot_Std_Ellipse %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Plot_GM_Landscape(W,M,V,R)
figure
[d,k] = size(M);
Rmin = min(R);
Rmax = max(R);
R = [Rmin:0.01*(Rmax-Rmin):Rmax];
Rlen= length(R);
P = zeros(Rlen,Rlen);
Q = zeros(Rlen,Rlen);
C = (2*pi)^(d/2);
for i=1:k,
    iV = inv(V(:,:,i));
    dV = sqrt(det(V(:,:,i)));
    for r=1:Rlen,
        for c=1:Rlen,
            dXM = [R(r);R(c)] - M(:,i);
            P(r,c) = exp(-0.5*dXM'*iV*dXM)/(C*dV);        
        end
    end    
    Q = Q + W(i)*P;
end
meshc(R,R,Q); colorbar('vert');
xlabel('1^{st} dimension');
ylabel('2^{nd} dimension');
title('Gaussian Mixture estimated by EM');
axis([Rmin Rmax Rmin Rmax])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Plot_GM_Landscape %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err,d] = VerifyWMV(W,M,V)
err = 0;
[Wd,Wk] = size(W);
[Md,Mk] = size(M);
[Vd1,Vd2,Vk] = size(V);
if Wd~=1,
    disp('W must be a row vector!');
    err = 1;
elseif Wk~=Mk | Wk~=Vk | Mk~=Vk,
    disp('The highest dimension of W, M and V must match!');
    err = 1;
elseif Md~=Vd1 | Md~=Vd2 | Vd1~=Vd2,
    disp('The random variable dimensions in M and V must match!');
    err = 1;
end
d = Md;
if d<1 | d>2,
    disp('Can only plot 1 or 2 dimensional applications!/n');
    err = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of VerifyWMV %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%

function err = VerifyF(flag)
err = 0;
if ~(flag=='e' | flag=='l' | flag=='b'),
    disp('flag must be ''e'', ''l'' or ''b''!');
    err = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of VerifyF %%%%
%%%%%%%%%%%%%%%%%%%%%%%%

function err = VerifyX(X,d)
err = 0;
[n,Xd] = size(X);
if Xd~=d,
    disp('Dimension of X must match the corresponding dimension of W, M and V');
    err = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of VerifyX %%%%
%%%%%%%%%%%%%%%%%%%%%%%%
