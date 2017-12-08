% Plot settings
ShowPlots = 1;
ShowUpdate = 1;

V_bounds = [0 10 0 10];

% Constant Velocity Model
Params_cv.xDim = 2;
Params_cv.q = 0.0001;
CVmodel = ConstantVelocityModelX(Params_cv);

% Constant Heading Model
Params.q_vel = 0.01;
Params.q_head = 0.3;
CHmodel = ConstantHeadingModelX(Params);

% Positional Observation Model
Params_meas.xDim = 4;
Params_meas.yDim = 2;
Params_meas.r = 0.5;
obs_model = PositionalObsModelX(Params_meas);

% Initiate Kalman Filter
Init.k = 1;
Init.DynModel = CHModel;
Init.ObsModel = ObsModel;
Init.x_init = [x_true(2,1); y_true(2,1); x_true(2,1)-x_true(1,1); y_true(2,1)-y_true(2,1)];
Init.P_init = CHModel.Params.Q(1);
ukf = UKalmanFilterX(Init);

% Containers
N = size(x_true,1)-2;
xV = zeros(4,N);          %estmate        % allocate memory
sV = zeros(2,N);          %actual
zV = zeros(2,N);
filtered_estimates = {};

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
end

% START OF SIMULATION
% ===================>
tic;
for k = 1:N
  
  % Generate new measurement from ground truth
  sV(:,k) = [x_true(k+2,1); y_true(k+2,1)];     % save ground truth
  zV(:,k) = obs_model.sample(0, sV(:,k),1);     % generate noisy measurment
  
  % Iterate Kalman Filter
  ukf.Params.y = zV(:,k);
  ukf.Predict(); 
  ukf.Update();
  
  xV(:,k) = ukf.Params.x;                    % save estimate
  filtered_estimates{k} = ukf.Params;
  
  % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1),zV(1,k),zV(2,k),'k*','MarkerSize', 10);
        h2 = plot(ax(1),sV(1,1:k),sV(2,1:k),'b.-','LineWidth',1);
        if j==2
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
        end
        h2 = plot(ax(1),sV(1,k),sV(2,k),'bo','MarkerSize', 10);
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        plot(xV(1,k), xV(2,k), 'ro', 'MarkerSize', 10);
        plot(xV(1,1:k), xV(2,1:k), 'r.-', 'MarkerSize', 10);
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
        pause(0.01);
    end
  %s = f(s) + q*randn(3,1);                % update process 
end
smoothed_estimates = ukf.Smooth(filtered_estimates);
toc;
% END OF SIMULATION
% ===================>

for i=1:N
    xV_smooth(:,i) = smoothed_estimates{i}.x;          %estmate        % allocate memory
end
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
    % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
    % Flip the image upside down before showing it
    % Plot data
    cla(ax(1));
     % Flip the image upside down before showing it
    imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

    % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
    hold on;
    h2 = plot(ax(1),zV(1,k),zV(2,k),'k*','MarkerSize', 10);
    h2 = plot(ax(1),sV(1,1:k),sV(2,1:k),'b.-','LineWidth',1);
    if j==2
        set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
    end
    h2 = plot(ax(1),sV(1,k),sV(2,k),'bo','MarkerSize', 10);
    set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
    plot(xV(1,k), xV(2,k), 'ro', 'MarkerSize', 10);
    plot(xV(1,1:k), xV(2,1:k), 'r.-', 'MarkerSize', 10);
    plot(xV_smooth(1,k), xV_smooth(2,k), 'go', 'MarkerSize', 10);
    plot(xV_smooth(1,1:k), xV_smooth(2,1:k), 'g.-', 'MarkerSize', 10);
    % set the y-axis back to normal.
    set(ax(1),'ydir','normal');
    str = sprintf('Robot positions (Update)');
    title(ax(1),str)
    xlabel('X position (m)')
    ylabel('Y position (m)')
    axis(ax(1),V_bounds)
    pause(0.01);
end

err(3,:) = sqrt((xV(1,:)-sV(1,:)).^2 + (xV(2,:)-sV(2,:)).^2);

% for k=1:2                                 % plot results
%   subplot(2,1,k)
%   plot(1:N, sV(k,:), '-', 1:N, xV(k,:), '--')
% end