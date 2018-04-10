function [Fnew,Qnew,Hnew,Rnew,Bnew,Params] = KalmanFilterX_LearnEM_Mstep_new(est_f,est_s,F,H,B,Q,R)
    
    if(nargin<5)
        Bnew = 0;
    end
    
    y       = cell2mat(cellfun(@(x)x.y,est_f,'un',0));    % Observations {N} x (yDim x 1)
    u       = cell2mat(cellfun(@(x)x.u,est_f,'un',0));    % Control inputs {N{ x (uDim x 1) 
    V_f     = cellfun(@(x)x.P,est_f,'un',0);    % Filtered co-variences {N} x (xDim x xDim)
    x_s     = cell2mat(cellfun(@(x)x.x,est_s,'un',0));    % Smoothed estimates {N} x (xDim x 1)
    V_s     = cellfun(@(x)x.P,est_s,'un',0);    % Smoothed co-variences {N} x (xDim x xDim)
    K       = cellfun(@(x)x.K,est_f,'un',0);    % Filtered Kalman gains {N} x (xDim x yDim)
    [xDim,N] = size(x_s);
    est_s{N}.C = est_s{N-1}.C;
    C       = cellfun(@(x)x.C,est_s,'un',0);    % Smoothed Kalman gains {N} x (xDim x yDim)
   
    
    Vt_tm1  = zeros(xDim,xDim,N);     % V_{t,t-1}
    Pt      = zeros(xDim,xDim,N);     % P_{t}
    Pt_tm1  = zeros(xDim,xDim,N);     % P_{t,t-1}
    
    
    Vt_tm1(:,:,N) = F*V_f{N-1} - K{N}*H*F*V_f{N-1} ;
    Pt_tm1(:,:,N) = Vt_tm1(:,:,N) + x_s(:,N)*x_s(:,N-1)';
    Pt(:,:,N)     = V_s{N} + x_s(:,N)*x_s(:,N)';
    Pt(:,:,N-1)     = V_s{N-1} + x_s(:,N-1)*x_s(:,N-1)';
    
    
    for k = N-1:-1:2
        Vt_tm1(:,:,k) = V_f{k}*C{k-1}' + C{k}*(Vt_tm1(:,:,k+1) - F*V_f{k})*C{k-1}';
        Pt_tm1(:,:,k) = Vt_tm1(:,:,k) + x_s(:,k)*x_s(:,k-1)';
        Pt(:,:,k-1)   = V_s{k-1} + x_s(:,k-1)*x_s(:,k-1)';
    end
    
%     %% Backward Rauch Recursion %Equation 2
%     [xhat, VtT, J, Pt2, Vtt1T, Ptt1] = Shuang_KF_Smooth(cell2mat(cellfun(@(x)x.x,est_f,'un',0)),cell2mat(cellfun(@(x)x.P,est_f,'un',0)),...
%             cell2mat(cellfun(@(x)x.PPred,est_f,'un',0)), F, H, B, [u(2:end),0], K{N});
%     sum(xhat == x_s)
%     test = permute(Pt_tm1,[2 3 1]);
%     sum(Ptt1 == test)
%     test = permute(Vt_tm1,[2 3 1]);
%     %sum(Vtt1T == test)
%     sum(Vtt1T - test<eps)
%     sum(cell2mat(C(1,1:end-1)) == J)
%     %cell2mat(C(1,1:end-1)) == J
%     test = permute(Pt,[2 3 1]);
%     sum(Pt2 == test)
    
%     Hnew = sum(y(:,1:N).*xhat(:,1:N))/sum(Pt2(1:N));
%     
%     %Update noise covariance
%     Rnew = 1/N*sum(y.*y - Hnew*xhat.*y);
% 
%     %Update dynamics matrix
%     Fnew = sum(Ptt1(2:N)-B*u(2:N).*xhat(1:N-1))/sum(Pt2(1:N-1));
% 
%     %new version of B and Q
%     Bnew = sum((xhat(2:N)-F*xhat(1:N-1)).*u(2:N))/sum((u(2:N).^2));
%     
%     %Update state noise covariance
%     Qnew = 1/(N-1)*sum(Pt2(2:N) + xhat(2:N).*u(2:N)*B - F*Ptt1(2:N) + F*xhat(1:N-1).*u(2:N)*B - B*u(2:N).*xhat(2:N) + B*(u(2:N).^2)*B);
     

    % In the calcualtions below we have used the following simplification:
    %   summ = y*x';
    % is the same as sum_{t=1}^{N}(y_{t}*x_{t}'), i.e:
    %   summ = 0;
    %   for t=1:N
    %       summ = summ + y(:,t)*x(:,t)';
    %   end
      
    %Update output matrix
    Hnew = y*x_s'/sum(Pt,3); 

    %Update noise covariance
    Rnew = 1/N*( y*y' - y*(x_s'*H') - (H*x_s)*y' + H*sum(Pt,3)*H'); %sum( - H*xhat.*y);
    
    %Update dynamics matrix
    Fnew = (sum(Pt_tm1(:,:,2:N),3)-(B*u(:,2:N))*x_s(:,1:N-1)') / sum(Pt(:,:,1:N-1),3);

    %update Control matrix
    Bnew = (x_s(:,2:N)*u(:,2:end)' - F*x_s(:,1:N-1)*u(:,2:N)') / (u(:,2:N)*u(:,2:N)');
    
    %Update state noise covariance
    Qnew = 1/(N-1)*(sum(Pt(:,:,2:end),3) - sum(Pt_tm1(:,:,2:end),3)*F' - F*sum(Pt_tm1(:,:,2:end),3)' + F*sum(Pt(:,:,1:N-1),3)*F' ...
         + F*(x_s(:,1:N-1)*u(:,2:N)')*B' + B*(u(:,2:end)*x_s(:,1:N-1)')*F'- x_s(:,2:N)*u(:,2:N)'*B' - B*(u(:,2:N)*x_s(:,2:N)')+ B*(u(:,2:N)*u(:,2:N)')*B'); 
    
     
    Params.loglik = computeEMLoglik(est_f{1}.x, est_f{1}.P, x_s,y,Fnew,Hnew,Qnew,Rnew,Bnew,u,Pt,Pt_tm1);
     
    Params.Vt_tm1 = Vt_tm1;
    Params.Pt_tm1 = Pt_tm1;
    Params.Pt = Pt;

end