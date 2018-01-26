% Plot settings
ShowPlots = 1;              % Enable|Disable plot outputs
SkipSimFrames = 10;          % Number of EM output plot frames to skip 
ShowPredict = 0;            % Show KF Prediction Step plot Output
ShowUpdate = 0;             % Show KF Update Step plot Output
nEMIter = 50000;            % Number of EM iterations

% Desired parameters
F_pre = @(~) 1;
H_pre = @(~) 1;
Q_pre = @(~) .1;
R_pre = @(~) 1;
B_pre = @(~) -.002;
c     = 10^(-6);

% Recording Settings
Record = 0;                 % Enable|Disable recording
clear Frames                % Clear stored frames from previous simulations

% Simulation Settings
N = 2733;                   % Number of timesteps

% Results Log Container
Log.xFilt = zeros(1,N);     % Log to store filtered state vectors (for each EM iteration)             
Log.exec_time = 0;          % Log to store total execution time
Log.estFilt = cell(1,N);    % Log to store all filtered estimates (for each EM iteration)
    

% Create figure windows
if(ShowPlots)
    figure('units','normalized','outerposition',[0 0 1 1])
    ax(1) = gca;
end

% Instantiate a generic dynamic model
Params_dyn.xDim   = 1;
Params_dyn.q      = 1;                          
DynModel          = GenericDynamicModelX(Params_dyn);
DynModel.Params.F = F_pre;     % Set Transition matrix F = 1;
DynModel.Params.Q = Q_pre;    % Set Process noise covariance Q = q^2

% Instatiate a generic observation model
% (H = 1, R = r^2)
Params_obs.xDim = 1;
Params_obs.yDim = 1;
Params_obs.r = 1;
ObsModel = GenericObservationModelX(Params_obs);
ObsModel.Params.H = H_pre;
ObsModel.Params.R = R_pre;

Params_ctr.xDim = 1;
Params_ctr.uDim = 1;
Params_ctr.B    = B_pre;
CtrModel = LinearNoiselessCtrModelX(Params_ctr); 


Q_old = DynModel.Params.Q(1);
F_old = DynModel.Params.F(1);
R_old = ObsModel.Params.R(1);

% % Generate ground truth and measurements
% sV = 5;
% zV = ObsModel.sample(0, sV(1),1);
% %uV = 2*rand()-1;
% clear pErr mErr;
% mErr = zV - sV;
% for k = 2:N
%     % Generate new measurement from ground truth
%     %uV(:,k) = 2*rand()-1;
%     pErr(k) = DynModel.sys_noise(1,1);
%     sV(:,k) = DynModel.sys(1,sV(:,k-1),pErr(k)) + CtrModel.ctr(1,uV(:,k));     % save ground truth
%     zV(:,k) = ObsModel.sample(0, sV(:,k),1);     % generate noisy measurment
%     mErr(k) = zV(k) - ObsModel.obs(0,sV(k));
%     
% end

for i = 1:size(uV,2)
    if uV(i) > 100
        uV(i) = 100;
    elseif uV(i) < -100
        uV(i) = -100;
    end
end

% Comput initial estimates
B_init = (tilde(zV,N-3)-tilde(zV,1))/(N-4);
z_tilde = zeros(1,N-3);
for k = 1:N-3
    z_tilde(k) = tilde(zV,k);
end

F_init = 1;
H_init = H_pre();
Q_init = 1/3 *(var(z_tilde(5:end) - z_tilde(1:end-4)) - var(z_tilde(2:end)-z_tilde(1:end-1)));
R_init = 1/2*(var(z_tilde(2:end)-z_tilde(1:end-1))-Q_init);
if(Q_init<0)
    Q_init = 0.1;
end
if(R_init<0)
    R_init = 0.1;
end
%F_init = mean( (H_init*B_init*(uV(2:end)-uV(1:end-1))-(zV(2:end)-zV(1:end)))*(H_init*
%F_init = mean((zV(:,2:end)-B_init*uV(:,2:end)-(zV(:,1:end-1)-B_init*uV(:,1:end-1)))/H);


% Calculate and store the true process and measurement noise covariances
% Q_true = var(pErr);
% R_true = std(mErr)^2;

% Corrupt the model parameters
DynModel.Params.F = @(~) F_init;
%ObsModel.Params.H = @(~) 1;
ObsModel.Params.R = @(~) R_init;
DynModel.Params.Q = @(~) Q_init;
CtrModel.Params.B = @() B_init;

% Initiate Kalman Filter
Params_kf.k        = 1;
Params_kf.x_init   = zV(1);
Params_kf.P_init   = 0.1;
Params_kf.DynModel = DynModel;
Params_kf.ObsModel = ObsModel;
Params_kf.CtrModel = CtrModel;
KFilter            = KalmanFilterX(Params_kf);

% For all EM iterations
for EMIter = 1:nEMIter
    
    fprintf('\nEMIter: %d/%d\n', EMIter, nEMIter);

    % FILTERING
    % ===================>
    
    % For all timesteps
    for k = 1:N
        
        % Update KF measurement vector
        KFilter.Params.y = zV(:,k);
        KFilter.Params.u = uV(:,k);
        
        if(k == 1)
            KFilter.Params.xPred = KFilter.Params.x;
            KFilter.Params.PPred = KFilter.Params.P;
            [KFilter.Params.yPred,KFilter.Params.S,KFilter.Params.Pxy] = ...
                KalmanFilterX_PredictObs(KFilter.Params.xPred,KFilter.Params.PPred,...
                            KFilter.ObsModel.obs(), KFilter.ObsModel.obs_cov());
            KFilter.Update();
        else
            % Iterate Kalman Filter
            KFilter.Iterate();
        end
        % Store Logs
        Log.xV(:,k)     = KFilter.Params.x;
        Log.estFilt{k}  = KFilter.Params;

        % Plot update step results
        if(ShowPlots && ShowUpdate)
            % Plot data
            cla(ax(1));
            hold on;
            h2 = plot(ax(1), k, zV(k),'k*','MarkerSize', 10);
            %h3 = plot(ax(1), 1:k, sV(1:k),'b.-','LineWidth',1);
            %plot(ax(1), k, sV(k),'bo','MarkerSize', 10);
            plot(k, Log.xV(k), 'o', 'MarkerSize', 10);
            h4 = plot(1:k, Log.xV(1:k), '.-', 'MarkerSize', 10);
            title(sp1,'\textbf{State evolution}','Interpreter','latex')
            xlabel(sp1,'Time (s)','Interpreter','latex')
            ylabel(sp1,'x (m)','Interpreter','latex')
            legend(sp1,[h2 h4], 'Measurements', 'Filtered state', 'Interpreter','latex');
            pause(0.01);
        end
    end
    estFilt = Log.estFilt;
    xV_filt = cell2mat(cellfun(@(x)x.x,estFilt,'un',0));
    
    % SMOOTHING
    % ===================>
    estSmooth = cell(1,N);
    estSmooth{N}.x = estFilt{N}.x;
    estSmooth{N}.P = estFilt{N}.P;
    for k = N-1:-1:1
        [estSmooth{k}.x, estSmooth{k}.P, estSmooth{k}.C] = KalmanFilterX_SmoothRTS_Single(estFilt{k}.x,estFilt{k}.P,estFilt{k+1}.xPred,estFilt{k+1}.PPred, estSmooth{k+1}.x, estSmooth{k+1}.P, KFilter.DynModel.sys());%, KFilter.CtrModel.ctr(), estFilt{k}.u);
    end
    %smoothed_estimates = FilterList{i}.Filter.Smooth(filtered_estimates);
    xV_smooth = zeros(1,N);
    PV_smooth = zeros(1,N);
    for k=1:N
        xV_smooth(:,k) = estSmooth{k}.x;          %estmate        % allocate memory
        PV_smooth(:,k) = estSmooth{k}.P;
    end
    
    % MAXIMIZATION
    % ===================>
    [F,Q,H,R,B,EMParams] = KalmanFilterX_LearnEM_Mstep(estFilt, estSmooth,KFilter.DynModel.sys(),KFilter.ObsModel.obs(),KFilter.CtrModel.ctr());
    
    % Update KF instance with optimised parameters
    KFilter = KalmanFilterX(Params_kf);
    %KFilter.DynModel.Params.F = @(~)F;
    KFilter.DynModel.Params.Q = @(~)Q;
    %KFilter.ObsModel.Params.H = @(~)H;
    KFilter.ObsModel.Params.R = @(~)R; %diag(diag(R));
    KFilter.CtrModel.Params.B= @(~)B;
    
    % Compute log-likelihood
    loglik_k = computeEMLoglik(Params_kf.x_init, Params_kf.P_init, xV_smooth,zV,KFilter.DynModel.Params.F(),KFilter.ObsModel.Params.H(),...
                                KFilter.DynModel.Params.Q(),KFilter.ObsModel.Params.R(),KFilter.CtrModel.Params.B(),uV);
    
    % Print debugging
    disp('F H B Q R :');
    disp([F H B Q R]);
    
    % Evaluate termination condition
    terminate = false;
    if(EMIter>1)
        [terminate,delta_loglik,positive] = testEMConverged(loglik_k,loglik_km1,c);
        if(~positive)
            warning('Loglik decreasing!!');
        end
        disp('delta_loglik:');
        disp(delta_loglik);
    end  
    loglik_km1 = loglik_k;
   
    
    
    % Plot results
    if(Record || (ShowPlots && (EMIter==1 || rem(EMIter,SkipSimFrames)==0)) || terminate)
        
        clf;
        sp1 = subplot(2,4,[1 2 3 4]);

        % Plot State evolution
        hold on;
        h6 = errorbar(sp1,1:k,xV_smooth,PV_smooth);
        h2 = plot(sp1,1:k,zV(1:k),'k*','MarkerSize', 10);
        %h3 = plot(sp1,1:k,sV(1:k),'b.-','LineWidth',1);
        %plot(sp1,k,sV(k),'bo','MarkerSize', 10);
        plot(sp1,k, Log.xV(k), 'ro', 'MarkerSize', 10);
        h4 = plot(sp1,1:k, Log.xV(1:k), 'r.-', 'MarkerSize', 10);
        plot(sp1,k, xV_smooth(k), 'go', 'MarkerSize', 10);
        h5 = plot(sp1,1:k, xV_smooth(1:k), 'g.-', 'MarkerSize', 10);
        title(sp1,'\textbf{State evolution}','Interpreter','latex')
        xlabel(sp1,'Time (s)','Interpreter','latex')
        ylabel(sp1,'x (m)','Interpreter','latex')
        legend(sp1,[h2 h4 h5 h6], 'Measurements', 'Filtered state', 'Smoothed state', 'Smoothed variance');                           
        
        sp3 = subplot(2,4,[5 6]);
        x = -5*sqrt(Q):10*sqrt(Q)/1000:5*sqrt(Q);
        y = mvnpdf(x',0,Q);
        plot(sp3,x,y);
        hold on;
        y = mvnpdf(x',0,KFilter.DynModel.Params.Q());
        plot(sp3,x,y);
        xlabel(sp3,'\textbf{Process noise $w_k \sim \mathcal{N}(0,Q)$}','Interpreter','latex');
        ylabel(sp3,'pdf($w_k$)','Interpreter','latex');
        title(sp3,'\textbf{True vs Estimated process noise pdf}','Interpreter','latex');
        
        sp4 = subplot(2,4,[7 8]);
        x = -5*sqrt(R):10*sqrt(R)/1000:5*sqrt(R);
        y = mvnpdf(x',0,R);
        plot(sp4,x,y);
        hold on;    
        y = mvnpdf(x',0,KFilter.ObsModel.Params.R());
        plot(sp4,x,y);
        xlabel(sp4,'\textbf{Measurement noise $\epsilon_k \sim \mathcal{N}(0,R)$}','Interpreter','latex');
        ylabel(sp4,'pdf($v_k$)','Interpreter','latex');
        title(sp4,'\textbf{True vs Estimated measurement noise pdf}','Interpreter','latex');
        
        pause(.01);
    end
    
    if(terminate)
        set(gcf,'color','g')
        title(sp1,'\textbf{State evolution (DONE)}','Interpreter','latex')
        break;
    end
end

% END OF SIMULATION
% ===================>

if(Record)
    Frames = Frames(2:end);
    vidObj = VideoWriter(sprintf('em_test.avi'));
    vidObj.Quality = 100;
    vidObj.FrameRate = 100;
    open(vidObj);
    writeVideo(vidObj, Frames);
    close(vidObj);
end
