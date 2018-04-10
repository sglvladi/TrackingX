% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps

load('example.mat');
tot_elapsed = 0;
% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update
TrackNum = size(x_true,2);

tot_ellapsed = 0;

% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.04,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% n_x = 4;      % state dimensions
% q = 0.01;     % std of process noise 
% n_y = 2;      % measurement dimensions
% r = 0.1;      % std of measurement noise
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
[DataList,x1,y1] = gen_obs_cluttered_multi3(TrackNum, x_true, y_true, 0.2, lambdaV, 1); 
N=size(DataList,2); % timesteps 

% Assign PHD parameter values
config.Model = ssm;
q = dyn.covariance();
BirthComponents.Means = [(V_bounds(2)-V_bounds(1))*rand(50,1), mvnrnd(zeros(50,1), 0.05'),(V_bounds(4)-V_bounds(3))*rand(50,1),mvnrnd(zeros(50,1), 0.05')]';
BirthComponents.Covars = repmat(10*q,1,1,50);
BirthComponents.Weights = 0.001*ones(1,50);
config.PriorComponents.Means = [];
config.PriorComponents.Covars = [];
config.PriorComponents.Weights = [];
config.BirthIntFcn = BirthComponents; % Uniform position and heading, Gaussian speed
config.ProbOfDeath = 0.005;
config.ProbOfDetection = 0.9;
config.ProbOfGating = 0.9;
config.ClutterRate = lambdaV/V;
config.Filter = UnscentedKalmanFilterX(ssm);

% Instantiate PHD filter
myphd = GM_PHDFilterX(config);

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
    myphd.MeasurementList = tempDataList; % New observations
    BirthComponents.Means = [(V_bounds(2)-V_bounds(1))*rand(50,1), mvnrnd(zeros(50,1), 0.05'),(V_bounds(4)-V_bounds(3))*rand(50,1),mvnrnd(zeros(50,1), 0.05')]';
    BirthComponents.Covars = repmat(10*q,1,1,50);
    BirthComponents.Weights = 0.001*ones(1,50);
    myphd.BirthIntFcn = BirthComponents; % New birth components

    tic;
    % Predict PHD filter
    myphd.predict();
            
    % Update PHD filter
    myphd.update();
    
    ellapsed = toc;
    tot_ellapsed = tot_ellapsed + ellapsed;
    fprintf("Estimated number of targets: %f\n", myphd.NumTargets);
    fprintf("Ellapsed time: %f\n\n", ellapsed);
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
        particles = zeros(2,0);
        numParticles = 10000;
        for i=1:myphd.NumComponents
            numParts = round(myphd.Components.Weights(i)/sum(myphd.Components.Weights)*numParticles);
            particles = [particles, mvnrnd(myphd.Components.Means([1,3],i)',myphd.Components.Covars([1,3],[1,3],i)',numParts)'];
        end
        %gm = gmdistribution(myphd.Components.Means([1,3],:)',myphd.Components.Covars([1,3],[1,3],:),myphd.Components.Weights/sum(myphd.Components.Weights));
        cla(ax(2), 'reset');
        %fsurf(ax(2),@(x,y)pdf(gm,[x y]),[0 10 0 10])
        [bandwidth,density,X,Y]=kde2d(particles');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        %camlight
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), particles(1,:), particles(2,:), '.')
        hold on;
        plot(ax(2), myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds]);
        str = sprintf('PHD intensity (Update)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
end