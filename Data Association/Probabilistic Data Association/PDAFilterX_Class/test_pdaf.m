%% Plot settings
ShowPlots = 1;
ShowArena = 1;
SkipFrames = 0;

%% Constant Velocity model
Params_cv.q = 0.0001;
Params_cv.xDim = 4;
CVmodel = ConstantVelocityModelX(Params_cv);

%% Constant Heading Model
Params_ch.q_vel = 0.01;
Params_ch.q_head = 0.3;
CHmodel = ConstantHeadingModelX(Params_ch);

%% 2D linear-Gaussian Observation Model
Params_pom.r = 0.2;
Params_pom.xDim = 4;
Params_pom.yDim = 2;
POmodel = PositionalObsModelX(Params_pom);

%% Assign PF parameter values
Params_pf.k               = 1;                   % initial iteration number
Params_pf.Np              = 5000;                 % number of particles
Params_pf.DynModel = CVmodel;
Params_pf.ObsModel = POmodel;
Params_pf.resampling_strategy = 'systematic_resampling';

%% Initiate Kalman Filter
Params_kf.k = 1;


%% Set TrackNum
TrackNum = 2;

%% Generate DataList
%[DataList,x1,y1] = gen_obs_cluttered_multi2(TrackNum, x_true, y_true, Params_pom.r, 2, 20, 1);

%% Get GroundTruth
for i=1:TrackNum
    GroundTruth{i} = [x_true(:,i), y_true(:,i)]; % ith target's GroundTruth
end

%% Initiate TrackList
for i=1:TrackNum
    Params_kf.x_init = [GroundTruth{i}(1,1);GroundTruth{i}(1,2); 0; 0];
    Params_kf.P_init = CVmodel.Params.Q(1); %blkdiag(POmodel.Params.R(1)/2, 2^2, 2*pi);%CVmodel.Params.Q(1);
    %FilterList{1}.Filter = KalmanFilterX(Params_kf, CVmodel, obs_model);
    Params_pf.gen_x0 =  @(Np) [mvnrnd(repmat([GroundTruth{i}(1,1),GroundTruth{i}(1,2), 0, 0],Np,1),CVmodel.Params.Q(1))]; %2*pi*rand(Np,1)];
    Params_kf.DynModel = CVmodel;
    Params_kf.ObsModel = POmodel;
    %Params_pf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan((y_true(2,1)-y_true(1,1))/(x_true(2,1)-x_true(1,1)))]', Np,1), CVmodel.Params.Q(1));
    TrackList{i}.TrackObj = EKalmanFilterX(Params_kf);
    %TrackList{i}.TrackObj = ParticleFilterX(Params_pf);%UKalmanFilterX(Params_kf, CHmodel, POmodel);% 
end

%% Initiate PDAF parameters
Params_pdaf.DataList = DataList{1}(:,:);  
Params_pdaf.TrackList = TrackList;
Params_pdaf.pDetect = 0.8;
Params_pdaf.pGate = 0.998; %0.999;
Params_pdaf.gateLevel = chi2inv(Params_pdaf.pGate,2);%14;

%% Instantiate JPDAF
pdaf = JPDAFilterX(Params_pdaf);
pdaf.Params.TrackList = TrackList;
pdaf.Params.DataList = DataList{1}(:,:); 

%% Instantiate Log to store output
N=size(DataList,2);
Logs = cell(1, 5); % 4 tracks
N = size(x_true,1)-2;
for i=1:TrackNum
    Logs{i}.sv = zeros(4,N);
    Logs{i}.xV = zeros(4,N);          %estmate        % allocate memory
    Logs{i}.err = zeros(2,N);
    Logs{i}.pos_err = zeros(1,N);
    Logs{i}.exec_time = 0;
    Logs{i}.filtered_estimates = cell(1,N);
end

% Create figure windows
if(ShowPlots)
    if(ShowArena)
        img = imread('maze.png');

        % set the range of the axes
        % The image will be stretched to this.
        min_x = 0;
        max_x = 10;
        min_y = 0;
        max_y = 10;

        % make data to plot - just a line.
        x = min_x:max_x;
        y = (6/8)*x;
    end
    
    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
end

exec_time = 0;
for i = 1:N
    % Remove null measurements   
    pdaf.Params.k = 1;
    pdaf.Params.DataList = DataList{i}(:,:);
    pdaf.Params.DataList( :, ~any(pdaf.Params.DataList,1) ) = [];
    
    for j=1:TrackNum
        %Logs{j}.xV_ekf(:,i) = pdaf.Params.TrackList{j}.TrackObj.Params.x;
        st = [x1(i,j); y1(i,j)];
        Logs{j}.sV_ekf(:,i)= st;
        %Logs{j}.eV_ekf(:,i) = (pdaf.Params.TrackList{j}.TrackObj.Params.x(1:2,1) - st).*(pdaf.Params.TrackList{j}.TrackObj.Params.x(1:2,1) - st);
    end
    tic;
    pdaf.Predict();
    %exec_time = 
    if (ShowPlots && i>1)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:TrackNum
                h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum
                colour = 'r';
                if(j==2)
                   colour = 'c';
                elseif (j==3)
                   colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,:),Logs{j}.xV_ekf(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                if(isa(pdaf.Params.TrackList{j}.TrackObj,'ParticleFilterX'))
                    c_mean = sum(repmat(pdaf.Params.TrackList{j}.TrackObj.Params.w,size(pdaf.Params.TrackList{j}.TrackObj.Params.particles,1),1).*pdaf.Params.TrackList{j}.TrackObj.Params.particles,2);
                    c_cov = [std(pdaf.Params.TrackList{j}.TrackObj.Params.particles(1,:),pdaf.Params.TrackList{j}.TrackObj.Params.w)^2,0;0,std(pdaf.Params.TrackList{j}.TrackObj.Params.particles(2,:),pdaf.Params.TrackList{j}.TrackObj.Params.w)^2];
                else
                    c_mean = pdaf.Params.TrackList{j}.TrackObj.Params.xPred;
                    c_cov = pdaf.Params.TrackList{j}.TrackObj.Params.PPred(1:2,1:2);
                end
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov);
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(pdaf.Params.TrackList{j}.TrackObj.pf.particles(1,:),pdaf.Params.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                % set the y-axis back to normal.
            set(gca,'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis([0 10 0 10])
            pause(0.01)
        end
    end
    
    pdaf.Update();
        
    %store Logs
    for j=1:TrackNum
        Logs{j}.xV_ekf(:,i) = pdaf.Params.TrackList{j}.TrackObj.Params.x;
        st = [x1(i,j); y1(i,j)];
        Logs{j}.sV_ekf(:,i)= st;
        Logs{j}.eV_ekf(:,i) = (pdaf.Params.TrackList{j}.TrackObj.Params.x(1:2,1) - st).*(pdaf.Params.TrackList{j}.TrackObj.Params.x(1:2,1) - st);
    end

    if (ShowPlots)
        if(i==1 || rem(i,SkipFrames+1)==0)
            % Plot data
            clf;
             % Flip the image upside down before showing it
            imagesc([min_x max_x], [min_y max_y], flipud(img));

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.

            hold on;
            for j=1:TrackNum
                h2 = plot(Logs{j}.sV_ekf(1,1:i),Logs{j}.sV_ekf(2,1:i),'b.-','LineWidth',1);
                if j==2
                    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
                end
                h2 = plot(Logs{j}.sV_ekf(1,i),Logs{j}.sV_ekf(2,i),'bo','MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            h2 = plot(DataList{i}(1,:),DataList{i}(2,:),'k*','MarkerSize', 10);
            for j=1:TrackNum
                colour = 'r';
                if(j==2)
                   colour = 'c';
                elseif (j==3)
                   colour = 'm';
                end
                h4 = plot(Logs{j}.xV_ekf(1,:),Logs{j}.xV_ekf(2,:),strcat(colour,'.-'),'LineWidth',1);
                %h4 = plot(Logs{j}.xV_ekf(1,i),Logs{j}.xV_ekf(2,i),strcat(colour,'o'),'MarkerSize', 10);
                if(isa(pdaf.Params.TrackList{j}.TrackObj,'ParticleFilterX'))
                    c_mean = sum(repmat(pdaf.Params.TrackList{j}.TrackObj.Params.w,size(pdaf.Params.TrackList{j}.TrackObj.Params.particles,1),1).*pdaf.Params.TrackList{j}.TrackObj.Params.particles,2);
                    c_cov = [std(pdaf.Params.TrackList{j}.TrackObj.Params.particles(1,:),pdaf.Params.TrackList{j}.TrackObj.Params.w)^2,0;0,std(pdaf.Params.TrackList{j}.TrackObj.Params.particles(2,:),pdaf.Params.TrackList{j}.TrackObj.Params.w)^2];
                else
                    c_mean = pdaf.Params.TrackList{j}.TrackObj.Params.x;
                    c_cov = pdaf.Params.TrackList{j}.TrackObj.Params.P(1:2,1:2);
                end
                h2=plot_gaussian_ellipsoid(c_mean(1:2), c_cov);
                set(h2,'color',colour);
                set(h2,'LineWidth',1);
                %plot(pdaf.Params.TrackList{j}.TrackObj.pf.particles(1,:),pdaf.Params.TrackList{j}.TrackObj.pf.particles(2,:),strcat(colour,'.'),'MarkerSize', 3);
                set(get(get(h4,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
                % set the y-axis back to normal.
            set(gca,'ydir','normal');
            str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
            title(str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
%            h_legend = legend('Real', 'Meas', 'Target 1', 'Target 2');
%            set(h_legend,'FontSize',9, 'Orientation', 'horizontal', 'Location', 'north');
            axis([0 10 0 10])
            pause(0.01)
        end
    end
end

    
ospa_vals= zeros(N,3);
ospa_c= 10;
ospa_p= 1;
for k=1:N
    trueX = [x_true(k,:);y_true(k,:)];
    estX = zeros(2,TrackNum);
    for i=1:TrackNum
        estX(:,i) = Logs{i}.xV(1:2,k);
    end
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= ospa_dist(trueX,estX,ospa_c,ospa_p);
end