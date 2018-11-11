% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps

% Load dataset
load('3_robots.mat');

tot_ellapsed = 0;
% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update
NumTracks = 3;

% Recording settings
clear F;
Record = 1;                 % Set to (0|1) to turn video recording (off|on)
FrameRate = 10;            % Number of frames per second
VideoQuality = 100;         % Set to desired quality percentage
VideoPathName = 'tomb_tracks_only.avi'; % Set to the desired path and name of produced recording

% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.01,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
lambda = lambdaV/V;

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
[DataList,x1,y1] = gen_obs_cluttered_multi2(NumTracks, x_true, y_true, sqrt(obs.ObsErrVariance), V_bounds, lambdaV, 1, 1); 
N=size(DataList,2); % timesteps 

% Assign PHD parameter values
config.Model = ssm;
config.BirthModel.ProbOfBirth = 0.005;
xb = [5;0;5;0];
Pb = diag([10,0.1,10,0.1]);
lambdab = 0.05;
config.BirthModel.BirthIntFcn = @ () deal(xb,Pb,lambdab);
config.ProbOfSurvive = 0.9;
config.ProbOfDetection = 0.9;
config.ClutterIntFcn = @(z) lambda;

% Instantiate PHD filter
filter = TrackOrientedMeMBerPoissonGMFilterX(config);
filter.Posterior.Poisson{1}.Mean = xb;
filter.Posterior.Poisson{1}.Covar = Pb;
filter.Posterior.Poisson{1}.Weight = 3;

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
    axis(ax(1), 'manual');
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;
    
    axis
end

% START OF SIMULATION
% ===================>
for k=1:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    
    % Extract DataList at time k
    tempDataList = DataList{k}(:,:);
    tempDataList( :, ~any(tempDataList,1) ) = [];       
    
    % Change PHD filter parameters
    filter.MeasurementList = tempDataList; % New observations
    
    tic;
    % Predict PHD filter
    filter.predict();
        
    % Update PHD filter
    filter.update();
    ellapsed = toc;
    tot_ellapsed = tot_ellapsed + ellapsed;
    %fprintf("Probability of existence: %f\n", filter.ProbOfExistence);
    fprintf("Ellapsed time: %f\n\n", ellapsed);
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        cla(ax(1));
         % Flip the image upside down before showing it

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        axis(ax(1),V_bounds)
        hold(ax(1),'on');

        hold(ax(1),'on');
%         for j=1:NumTracks
%             h2 = plot(GroundTruth{k}(1,1:i),Logs{j}.Groundtruth.State(2,1:i),'b.-','LineWidth',1);
%             if j==2
%                 set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
%             end
%             h2 = plot(Logs{j}.Groundtruth.State(1,i),Logs{j}.Groundtruth.State(2,i),'bo','MarkerSize', 10);
%             set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % 
%         end
        for i=1:numel(filter.Posterior.Bernoulli)
            %if(numel(filter.Posterior.Bernoulli{i}.Trajectory.ProbOfExistence)>10&&mean(filter.Posterior.Bernoulli{i}.Trajectory.ProbOfExistence)>0.8)
            if(filter.Posterior.Bernoulli{i}.ProbOfExistence>0.8)
                plot(ax(1),filter.Posterior.Bernoulli{i}.Trajectory.Mean(1,:),filter.Posterior.Bernoulli{i}.Trajectory.Mean(3,:),'.-');
                %plot_gaussian_ellipsoid(filter.Posterior.Bernoulli{i}.Trajectory.Mean([1,3],end),filter.Posterior.Bernoulli{i}.Trajectory.Covariance([1,3],[1,3],end),'r',1,50,ax(1));
            end
                %hold on;
        end
        
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
            
        % Plot PHD
%         cla(ax(2), 'reset');
%         [bandwidth,density,X,Y]=kde2d(filter.Poisson.Updated.Particles([1,3],:)');
%         %contour3(X,Y,density,50);
%         h = surf(ax(2),X,Y,density);        
%         shading interp
%         colormap(ax(2), jet(3000))
%         %set(h, 'edgecolor','none')
%         hold on;
%         plot(ax(2), filter.MeasurementList(1,:), filter.MeasurementList(2,:), 'y*');
%         axis(ax(2), [V_bounds]);
%         str = sprintf('PHD intensity (Update)');
%         xlabel(ax(2),'X position (m)')
%         ylabel(ax(2),'Y position (m)')
%         zlabel(ax(2),'Intensity')
%         title(ax(2),str)
%         pause(0.01)
        
        % Store video frame
        if(Record)
            F(k) = getframe(ax(1));
        end
    end
end

% Create video file and write to it
if(Record)
    %F = F(2:end);
    vidObj = VideoWriter(char(VideoPathName));
    vidObj.Quality = VideoQuality;
    vidObj.FrameRate = FrameRate;
    open(vidObj);
    writeVideo(vidObj, F);
    close(vidObj);
end
