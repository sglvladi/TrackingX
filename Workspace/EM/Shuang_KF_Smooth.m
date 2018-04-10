function [xhat, VtT, J, Pt, Vtt1T, Ptt1] = Shuang_KF_Smooth(xtt,Vtt,Vtt1,F,H,B,u,KT)

T = length(xtt);
%% Backward Rauch Recursion %Equation 2
xtT(T) = xtt(T);
VtT(T) = Vtt(T);
for t = T:-1:2
    J(t-1) = Vtt(t-1)*F'/Vtt1(t);
    xtT(t-1) = xtt(t-1) + J(t-1)*(xtT(t)-(F*xtt(t-1)+B*u(t-1)));
    VtT(t-1) = Vtt(t-1) + J(t-1)*(VtT(t)-Vtt1(t))*J(t-1);
end

%% Another Recursion
xhat = xtT;
Pt = VtT + xtT.^2; %size 1*T
Vtt1T(T) = (1-H*KT(T))*F*Vtt(T-1);
for t = T:-1:3
   Vtt1T(t-1) = Vtt(t-1)*J(t-2) + J(t-1)*(Vtt1T(t)-Vtt(t-1))*J(t-2);  %Equation 3
end
%Ptt1=[NaN Vtt1T(2:T)+xtT(2:T).*xtT(1:(T-1))];
Ptt1 = [0 (Vtt1T(2:T)+xtT(2:T).*xtT(1:T-1))]; %Equation 4

end