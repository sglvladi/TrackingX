% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps

% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update
TrackNum = size(x_true,2);

% Simulation parameters
% Constant Velocity Model
ParamsCV.q = 0.0001;
ParamsCV.xDim = 4;
CVmodel = ConstantVelocityModelX(ParamsCV);

% % Constant Heading Model
Params.q_vel = 0.01;
Params.q_head = 0.3;
CHmodel = ConstantHeadingModelX(Params);
% 
% % Constant Heading Model
% Params.q_vel = 0.01;
% Params.q_head = 0.3;
% CHmodel2 = ConstantHeadingModel2X('q_vel',0.01,'q_head',0.16,'smartArgs',false);

% Positional Observation Model
Params_meas.xDim = 4;
Params_meas.yDim = 2;
Params_meas.r = .2;
obs_model = PositionalObsModelX(Params_meas);

% n_x = 4;      % state dimensions
% q = 0.01;     % std of process noise 
% n_y = 2;      % measurement dimensions
% r = 0.1;      % std of measurement noise
lambdaV = 50; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
[DataList,x1,y1] = gen_obs_cluttered_multi3(TrackNum, x_true, y_true, 0.2, lambdaV, 1); 
N=size(DataList,2); % timesteps 

% % Constant Heading (CH) dynamic model
% F_ch = @(k,xkm1) [xkm1(1,:)+k*xkm1(3,:).*cos(xkm1(4,:)); xkm1(2,:)+k*xkm1(3,:).*sin(xkm1(4,:)); xkm1(3,:); xkm1(4,:)];
% Q_ch = diag([0,0,q^2,0.16^2]);
% gen_sys_noise_ch = @(Np) mvnrnd(zeros(Np, n_x), Q_ch)';
% sys_ch = @(k, xkm1, uk) F_ch(k,xkm1) + uk; 
% 
% % Constant Velocity (CV) dynamic model
% F_cv = @(k,xkm1) [xkm1(1,:)+k*xkm1(3,:); xkm1(2,:)+k*xkm1(4,:); xkm1(3,:); xkm1(4,:)];
% Q_cv = [1^3/3, 0, 1^2/2, 0;  0, 1^3/3, 0, 1^2/2; 1^2/2, 0, 1, 0; 0, 1^2/2, 0, 1]*q^2;
% sys_cv = @(k, xkm1, uk) F_cv(k,xkm1) + uk;
% gen_sys_noise_cv = @(Np) mvnrnd(zeros(Np, n_x), Q_cv)';
% 
% % Linear Observation model
% H = @(xk) [xk(1,:); xk(2,:)];
% R = [r^2 0; 0 r^2];
% likelihood = @(k, yk, xk) mvnpdf(yk', H(xk)', R); % p(y_k|x_k)

% Assign PHD parameter values
config.k               = 1;                   % initial iteration number
config.Np              = 50000;              % number of particles
config.resampling_strategy = 'systematic_resampling'; % resampling strategy
config.DynModel = CVmodel;
config.ObsModel = obs_model;
%config.sys = sys_ch; 
%config.sys_noise = gen_sys_noise_ch; 
config.gen_x0 = @(Np) [(V_bounds(2)-V_bounds(1))*rand(Np,1),(V_bounds(4)-V_bounds(3))*rand(Np,1), mvnrnd(zeros(Np,2), eye(2)*CVmodel.Params.q^2')]; % Uniform position and heading, Gaussian speed
%config.gen_x0 = @(Np) [(V_bounds(2)-V_bounds(1))*rand(Np,1),(V_bounds(4)-V_bounds(3))*rand(Np,1), mvnrnd(zeros(Np,1), CVmodel.Params.q^2), 2*pi*rand(Np,1)]; % Uniform position and heading, Gaussian speed
config.particles_init = config.gen_x0(config.Np)'; % Generate inital particles as per gen_x0
config.w_init = repmat(1/config.Np, config.Np, 1)'; % Uniform weights
%config.likelihood = likelihood;
%config.obs_model = H;
config.pBirth = 0.1;
config.pDeath = 0.005;
config.Jk = 500;
config.pDetect = 0.9;
config.lambda = lambdaV/V;
config.pConf = 0.9;
config.NpConf = 1000;
config.type = 'search';
config.birth_strategy = 'mixture';

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config);

% Create figure windows
if(ShowPlots)
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

    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;
end

% START OF SIMULATION
% ===================>
for k=1:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    
    % Extract DataList at time k
    tempDataList = DataList{k}(:,:);
    tempDataList( :, ~any(tempDataList,1) ) = [];       
    
    % Change PHD filter parameters
    myphd.Params.k = 1; % Time index
    myphd.Params.y = tempDataList; % New observations
    %myphd.config.lambda = size(tempDataList,2)/V;
    
    % Predict PHD filter
    myphd.Predict();
    
    % Plot prediction step results
    if(ShowPlots && ShowPrediction)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        for j=1:TrackNum
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end

        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Prediction)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.Params.particles(1:2,:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);
        shading interp
        colormap(ax(2), jet(3000))
        set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.Params.particles(1,:), myphd.Params.particles(2,:), '.')
        hold on;
        plot(ax(2), myphd.Params.y(1,:), myphd.Params.y(2,:), 'y*');
        axis(ax(2), [V_bounds 0 10]);
        str = sprintf('PHD intensity (Prediction)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
        
    % Update PHD filter
    myphd.Update();
    fprintf("Estimated number of targets: %f\n\n", myphd.Params.Nk);
    
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        for j=1:TrackNum
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
            
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.Params.particles(1:2,:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.Params.particles(1,:), myphd.Params.particles(2,:), '.')
        hold on;
        plot(ax(2), myphd.Params.y(1,:), myphd.Params.y(2,:), 'y*');
        axis(ax(2), [V_bounds 0 20]);
        str = sprintf('PHD intensity (Update)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
end