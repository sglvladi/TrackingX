% Plot settings
ShowPlots = 1;
ShowUpdate = 1;
ShowArena = 0;
ShowPredict = 0;
SimNum = 100;
V_bounds = [0 25 0 15];

% Instantiate a Tracklist to store each filter
FilterList = [];
FilterNum = 3;

% Containers
Logs = cell(1, 3); % 4 tracks
N = size(x_true,1)-2;
for i=1:FilterNum
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



for SimIter = 1:SimNum
    
    FilterList = [];
    % Constant Velocity Model
    Params_cv.q = .0001;
    Params_cv.xDim = 4;
    CVmodel = ConstantVelocityModelX(Params_cv);

    % Constant Heading Model
    Params_ch1.q_vel = 0.01;
    Params_ch1.q_head = 0.3;
    CHmodel = ConstantHeadingModelX(Params_ch1);
    
%     % Constant Heading Model
%     Params_ch2.q_vel = 0.01;
%     Params_ch2.q_head = 0.3;
%     CHmodel2 = ConstantHeadingModel2X(Params_ch2);

    % Positional Observation Model
%     Params_meas.xDim = 4;
%     Params_meas.yDim = 2;
%     Params_meas.r = 0.3;
%     obs_model = PositionalObsModelX(Params_meas);
    
    % Polar2Cart Observation Model
    Params_meas.xDim = 4;
    Params_meas.R = [0.0001 0 ; 0 0.0001];
    obs_model = Polar2CartGaussModelX(Params_meas);

    % Initiate Kalman Filter
%     Params_kf.k = 1;
%     Params_kf.x_init = [x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)];
%     Params_kf.P_init = CVmodel.Params.Q(1);
%     Params_kf.DynModel = CVmodel;
%     Params_kf.ObsModel = obs_model;
%     FilterList{1}.Filter = KalmanFilterX(Params_kf);

    % Initiate Extended Kalman Filter
%     Params_ekf.k = 1;
%     %Params_ekf.x_init = [x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)];
%     Params_ekf.x_init = [x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan2((y_true(2,1)-y_true(1,1)),(x_true(2,1)-x_true(1,1)))];
%     Params_ekf.P_init = CHmodel.Params.Q(1);
%     Params_ekf.DynModel = CHmodel;
%     Params_ekf.ObsModel = obs_model;
%     FilterList{end+1}.Filter = EKalmanFilterX(Params_ekf);

    % Initiate Unscented Kalman Filter
    Params_ukf.k = 1;
    %Params_ukf.x_init = [x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)];
    Params_ukf.x_init = [x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan2((y_true(2,1)-y_true(1,1)),(x_true(2,1)-x_true(1,1)))];
    Params_ukf.P_init = CHmodel.Params.Q(1);
    Params_ukf.DynModel = CHmodel;
    Params_ukf.ObsModel = obs_model;
    FilterList{end+1}.Filter = UKalmanFilterX(Params_ukf);

    % Initiate Particle Filter
    Params_pf.k = 1;
    Params_pf.Np = 5000;
    %Params_pf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)]', Np,1), CVmodel.Params.Q(1));
    Params_pf.DynModel = CHmodel;
    Params_pf.ObsModel = obs_model;
    Params_pf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan2((y_true(2,1)-y_true(1,1)),(x_true(2,1)-x_true(1,1)))]', Np,1), CVmodel.Params.Q(1));
    FilterList{end+1}.Filter = ParticleFilterX(Params_pf);
    
    % Initiate Particle Filter
%     Params_pf.k = 1;
%     Params_pf.Np = 5000;
%     Params_pf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)]', Np,1), CVmodel.Params.Q(1));
%     %Params_pf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan((y_true(2,1)-y_true(1,1))/(x_true(2,1)-x_true(1,1)))-pi]', Np,1), CVmodel.Params.Q(1));
%     FilterList{5}.Filter = ParticleFilterX(CHmodel2, obs_model, 'Np', 5000, 'gen_x0', Params_pf.gen_x0);
 
    % Initiate EParticle Filter
%     Params_epf.k = 1;
%     Params_epf.Np = 5000;
%     %Params_epf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)]', Np,1), CVmodel.Params.Q(1));
%     Params_epf.DynModel = CHmodel;
%     Params_epf.ObsModel = obs_model;
%     Params_epf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan((y_true(2,1)-y_true(1,1))/(x_true(2,1)-x_true(1,1)))]', Np,1), CVmodel.Params.Q(1));
%     FilterList{end+1}.Filter = EParticleFilterX(Params_epf);
     
    % Initiate UParticle Filter
    Params_upf.k = 1;
    Params_upf.Np = 5000;
    %Params_upf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(1,1)]', Np,1), CVmodel.Params.Q(1));
    Params_upf.DynModel = CHmodel;
    Params_upf.ObsModel = obs_model;
    Params_upf.gen_x0 = @(Np) mvnrnd(repmat([x_true(2,1); y_true(2,1); sqrt((x_true(2,1)-x_true(1,1))^2+(y_true(2,1)-y_true(1,1))^2); atan2((y_true(2,1)-y_true(1,1)),(x_true(2,1)-x_true(1,1)))]', Np,1), CVmodel.Params.Q(1));
    FilterList{end+1}.Filter = UParticleFilterX(Params_upf);


    DataList = dataGen(obs_model,GroundTruth,1);
    % Generate ground truth and measurements
    for k = 1:N
        % Generate new measurement from ground truth
        sV(:,k) = [x_true(k+2,1); y_true(k+2,1)];     % save ground truth
        zV(:,k) = DataList{k}; %obs_model.sample(0, sV(:,k),1);     % generate noisy measurment
    end
    % START OF SIMULATION
    % ===================>
    tic;
    for k = 1:N

        % Update measurements
        for i=1:FilterNum
            FilterList{i}.Filter.Params.y = zV(:,k);
        end

        % Iterate all filters
        for i=1:FilterNum
            tic;
            FilterList{i}.Filter.Iterate();
            Logs{i}.exec_time = Logs{i}.exec_time + toc;
        end

        % Store Logs
        for i=1:FilterNum
            Logs{i}.err(:,k) = Logs{i}.err(:,k) + (sV(:,k) - FilterList{i}.Filter.Params.x(1:2))/SimNum;
            Logs{i}.pos_err(1,k) = Logs{i}.pos_err(1,k) + (sV(1,k) - FilterList{i}.Filter.Params.x(1))^2 + (sV(2,k) - FilterList{i}.Filter.Params.x(2))^2;
            Logs{i}.xV(:,k) = FilterList{i}.Filter.Params.x;
            Logs{i}.filtered_estimates{k} = FilterList{i}.Filter.Params;
        end

      % Plot update step results
        if(ShowPlots && ShowUpdate)
            % Plot data
            cla(ax(1));

            if(ShowArena)
                 % Flip the image upside down before showing it
                imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));
            end

            % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
            hold on;
            h2 = plot(ax(1),zV(2,k)*sin(zV(1,k)),zV(2,k)*cos(zV(1,k)),'k*','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),sV(1,1:k),sV(2,1:k),'b.-','LineWidth',1);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            h2 = plot(ax(1),sV(1,k),sV(2,k),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend

            for i=1:FilterNum
                h2 = plot(Logs{i}.xV(1,k), Logs{i}.xV(2,k), 'o', 'MarkerSize', 10);
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
                %plot(pf.Params.particles(1,:), pf.Params.particles(2,:), 'b.', 'MarkerSize', 10);
                plot(Logs{i}.xV(1,1:k), Logs{i}.xV(2,1:k), '.-', 'MarkerSize', 10);
            end
            legend('UKF', 'PF', 'UPF')
            %legend('KF','EKF', 'UKF', 'PF', 'EPF', 'UPF')

            if(ShowArena)
                % set the y-axis back to normal.
                set(ax(1),'ydir','normal');
            end

            str = sprintf('Robot positions (Update)');
            title(ax(1),str)
            xlabel('X position (m)')
            ylabel('Y position (m)')
            axis(ax(1),V_bounds)
            pause(0.01);
        end
      %s = f(s) + q*randn(3,1);                % update process 
    end
end

figure
for i=1:FilterNum
    hold on;
    plot(sqrt(Logs{i}.pos_err(1,:)/SimNum), '.-');
end
legend('UKF', 'PF', 'UPF')
%legend('KF','EKF', 'UKF', 'PF', 'EPF', 'UPF');%, 'EPF', 'UPF')

figure
bars = zeros(1, FilterNum);
c = {'UKF', 'PF', 'UPF'};
%c = {'KF','EKF', 'UKF', 'PF', 'EPF', 'UPF'};
%c = categorical(c, {'KF','EKF', 'UKF', 'PF', 'EPF', 'UPF'},'Ordinal',true); %, 'EPF', 'UPF'
c = categorical(c, {'UKF', 'PF', 'UPF'},'Ordinal',true); %, 'EPF', 'UPF'
for i=1:FilterNum
    bars(i) =  Logs{i}.exec_time;
end
bar(c, bars);
%smoothed_estimates = pf.Smooth(filtered_estimates);
toc;
% END OF SIMULATION
% ===================>

