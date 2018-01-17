function [Fnew,Qnew,Hnew,Rnew,Bnew,Params] = KalmanFilterX_LearnEM_Mstep(est_f,est_s,F,H,B)
    
    if(nargin<5)
        Bnew = 0;
    end
    
    
    N       = numel(est_f);                     % Number of timesteps
    y       = cellfun(@(x)x.y,est_f,'un',0);    % Observations {N} x (yDim x 1)
    u       = cellfun(@(x)x.u,est_f,'un',0);    % Control inputs {N{ x (uDim x 1) 
    V_f     = cellfun(@(x)x.P,est_f,'un',0);    % Filtered co-variences {N} x (xDim x xDim)
    x_s     = cellfun(@(x)x.x,est_s,'un',0);    % Smoothed estimates {N} x (xDim x 1)
    V_s     = cellfun(@(x)x.P,est_s,'un',0);    % Smoothed co-variences {N} x (xDim x xDim)
    K       = cellfun(@(x)x.K,est_f,'un',0);    % Filtered Kalman gains {N} x (xDim x yDim)
    est_s{N}.C = est_s{N-1}.C;
    C       = cellfun(@(x)x.C,est_s,'un',0);    % Smoothed Kalman gains {N} x (xDim x yDim)
    Vt_tm1  = cell(1,N);                        % V_{t,t-1}
    Pt      = cell(1,N);                        % P_{t}
    Pt_tm1  = cell(1,N);                        % P_{t,t-1}
    
    
    Vt_tm1{N} = F*V_f{N-1} - K{N}*H*F*V_f{N-1} ;
    Pt_tm1{N} = Vt_tm1{N} + x_s{N}*x_s{N-1}';
    Pt{N}     = V_s{N} + x_s{N}*x_s{N}';
    Pt{N-1}   = V_s{N-1} + x_s{N-1}*x_s{N-1}';

    % Sum components required to compute optimal F
    %  - sumF{1} = sum(P_{t,t-1})
    %  - sumF{2} = sum(P_{t-1})
    sumF = cell(1,2);
    sumF{1} = Pt_tm1{N} - B*u{N}*x_s{N-1}';
    sumF{2} = Pt{N-1}; 
    
    % Sum components required to compute optimal Q
    sumQ = Pt{N} - Pt_tm1{N}*F' - F*Pt_tm1{N}' + F*Pt{N-1}*F' ...
           + F*x_s{N-1}*u{N}'*B' + B*u{N}*x_s{N-1}'*F'- x_s{N}*u{N}'*B' - B*u{N}*x_s{N}'+ B*u{N}*u{N}'*B';           
    
    % Sum components required to compute optimal H
    %  - sumH{1} = sum(y_t * x_{t}^{s})
    %  - sumH{2} = sum(P_t)
    sumH = cell(1,2);
    sumH{1} = y{N}*x_s{N}';
    sumH{2} = Pt{N};
    
    % Sum components required to compute optimal R
    sumR = y{N}*y{N}' - y{N}*x_s{N}'*H' - H*x_s{N}*y{N}' + H*Pt{N}*H';
    
    % Sum components required to compute optimal B
    sumB = cell(1,2);
    sumB{1} = x_s{N}*u{N}' - F*x_s{N-1}*u{N}';
    sumB{2} = u{N}*u{N}';
    
    for k = N-1:-1:1
        if(k>1)
            Vt_tm1{k} = V_f{k}*C{k-1}' + C{k}*(Vt_tm1{k+1} - F*V_f{k})*C{k-1}';
            Pt_tm1{k} = Vt_tm1{k} + x_s{k}*x_s{k-1}';
            Pt{k-1}   = V_s{k-1} + x_s{k-1}*x_s{k-1}';
        
            sumF{1} = sumF{1} + Pt_tm1{k} - B*u{k}*x_s{k-1}';
            sumF{2} = sumF{2} + Pt{k-1};
        
            sumQ = sumQ + Pt{k} - Pt_tm1{k}*F' - F*Pt_tm1{k}' + F*Pt{k-1}*F' ...
                   + F*x_s{k-1}*u{k}'*B' + B*u{k}*x_s{k-1}'*F'- x_s{k}*u{k}'*B' - B*u{k}*x_s{k}'+ B*u{k}*u{k}'*B'; 
            
            sumB{1} = sumB{1} + x_s{k}*u{k}' - F*x_s{k-1}*u{k}';
            sumB{2} = sumB{2} + u{k}*u{k}';
        end
        
        sumH{1} = sumH{1} + y{k}*x_s{k}';
        sumH{2} = sumH{2} + Pt{k};
        
        sumR = sumR + y{k}*y{k}' - y{k}*x_s{k}'*H' - H*x_s{k}*y{k}' + H*Pt{k}*H';
    end
    
    Fnew = sumF{1}/sumF{2};%sumPt_tm1/sumPtm1;
    Hnew = sumH{1}/sumH{2};
    Bnew = sumB{1}/sumB{2};
    Rnew = sumR/N; Rnew = (Rnew+Rnew')/2;
    Qnew = sumQ/(N-1); Qnew = (Qnew+Qnew')/2;
    Params.Vt_tm1 = Vt_tm1;
    Params.Pt_tm1 = Pt_tm1;
    Params.Pt = Pt;
end