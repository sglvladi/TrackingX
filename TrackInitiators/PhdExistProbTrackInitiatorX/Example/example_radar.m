% phd_test.m
% ====================================================>
% This is a test script which demonstrates the usage of the "SMC_PHD" class.
% 
% SETUP:
%  * Before running the simulation, open "2_crossing_targets.mat" or "3_roaming_targets.mat" datasets, from the "datasets" folder
%  * The datasets have been extracted by simulating the motion of differential robots in a 2D-plane (x,y)
%  * The "gen_obs_cluttered_multi3" function takes as an input the ground truth data, including information about the measurement noise and clutter rate
%     and then produces 1xNk cell array of corrupted and cluttered measurements, Nk being the total number of timesteps
clear F;
load('example.mat');
tot_elapsed = 0;
% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
ShowPrediction = 0;         % Set to 0 to skip showing prediction
ShowUpdate = 1;             % Set to 0 to skip showing update
TrackNum = size(x_true,2);

tot_ellapsed = 0;

% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',1);
dyn.TimeVariant = 2;

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',50,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
%V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [-8154.72944624983 -212.289393440959 -7548.44272179096 4355.32645897434];%[-7548.44272179096,4361.40291082772,-8154.72944624983,5210.48002192510];%[-2500 200 -3000 3000]; % [x_min x_max y_min y_max]
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3)));

% Generate observations (Poisson number with rate of lambdaV, positions are uniform over surveillance region)
%[DataList,x1,y1] = gen_obs_cluttered_multi2(numTrueTracks, x_true, y_true, sqrt(obs.covariance()), V_bounds, lambdaV, 1); 
N=size(DataList,2); % timesteps 

% Assign PHD parameter values
config.NumParticles = 50000;              % number of particles
config.Model = ssm;
q = dyn.covariance();
%config.BirthIntFcn = @(Np) [(V_bounds(2)-V_bounds(1))*rand(Np,1), mvnrnd(zeros(Np,1), q(3,3)'),(V_bounds(4)-V_bounds(3))*rand(Np,1),mvnrnd(zeros(Np,1), q(4,4)')]'; % Uniform position and heading, Gaussian speed
config.BirthIntFcn = @(Np)[abs(V_bounds(2)-V_bounds(1))*rand(Np,1)+V_bounds(1), 3^2*rand(Np,1), abs(V_bounds(4)-V_bounds(3))*rand(Np,1)+V_bounds(3), 3^2*rand(Np,1)]'; % Uniform position and heading, Gaussian speed

config.PriorDistFcn = @ (Np) deal(config.BirthIntFcn(Np), repmat(1/Np, Np, 1)');
config.BirthScheme = {'Mixture', 0.01};
%config.BirthScheme = {'Expansion', 500};
config.ProbOfDeath = 0.005;
config.ProbOfDetection = 0.7;
config.ClutterRate = lambdaV/V;
config.ResamplingScheme = 'Multinomial';

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config);

% Instantiate PF filter
mypf = ParticleFilterX(ssm);

% Initiate PDAF parameters
Params_pdaf.Clusterer = NaiveClustererX();
Params_pdaf.Gater = EllipsoidalGaterX(2,'ProbOfGating',0.99)';
Params_pdaf.ProbOfDetect = 0.7;
mypdaf = JointProbabilisticDataAssocX(Params_pdaf);
mypdaf.MeasurementList = DataList{1}(:,:); 

% Initiate Track Initiator
config_ti.Filter = mypf;
config_ti.PHDFilter = myphd;
config_ti.ProbOfGating = 0.99;
config_ti.ProbOfConfirm = 0.8;
myti = PhdExistProbTrackInitiatorX(config_ti);

TrackList = [];

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
    plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2);
    axis([-4.19499354141717 -4.10940387037447 50.3123450958395 50.3735587491040])
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;
    
    plots = [];
end

% START OF SIMULATION
% ===================>
for k=1:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    if(k==92)
    
    end
    % Extract DataList at time k
    tempDataList = DataList{k}(:,:);
    tempDataList( :, ~any(tempDataList,1) ) = [];
    
    % Process JPDAF
    mypdaf.MeasurementList = tempDataList;
    mypdaf.TrackList = TrackList;
    for j=1:numel(TrackList)
        mypdaf.TrackList{j}.Filter.predict();
    end
    mypdaf.associate();    
    mypdaf.updateTracks();
    
    % Perform Track initiation
    myti.MeasurementList = tempDataList; % New observations
    myti.PHDFilter.ClutterRate = size(myti.MeasurementList,2)/V;
    myti.TrackList = TrackList;
    myti.AssocWeightsMatrix = mypdaf.AssocWeightsMatrix;
    TrackList = myti.initiateTracks();
    tic;
    
    % Append state to target trajectories
    for j=1:numel(TrackList)
        if(isempty(TrackList{j}.Trajectory))
            TrackList{j}.Trajectory = TrackList{j}.Filter.StateMean;
        else
            TrackList{j}.Trajectory = [TrackList{j}.Trajectory, TrackList{j}.Filter.StateMean];
        end
    end
    
    % Plot prediction step results
    if(ShowPlots && ShowPrediction)
        % Plot data
        cla(ax(1));
%         for i = 1:numel(plots)
%             delete(plots(i))
%         end
         plots = [];
         % Flip the image upside down before showing it
        imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        plots(end+1) = h2;
        for j=1:TrackNum
            h2 = plot(ax(1), x_true(1:k,j),y_true(1:k,j),'b.-','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            plots(end+1) = h2;
            h2 = plot(ax(1), x_true(k,j),y_true(k,j),'bo','MarkerSize', 10);
            set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            plots(end+1) = h2;
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
        [bandwidth,density,X,Y]=kde2d(myphd.PredParticles([1,3],:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);
        shading interp
        colormap(ax(2), jet(3000))
        set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.PredParticles(1,:), myphd.PredParticles(3,:), '.')
        hold on;
        plot(ax(2), myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds 0 10]);
        str = sprintf('PHD intensity (Prediction)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
    end
        
    % Update PHD filter
    myphd.update();
    ellapsed = toc;
    tot_ellapsed = tot_ellapsed + ellapsed;
    fprintf("Estimated number of targets: %f\n", myphd.NumTargets);
    fprintf("Ellapsed time: %f\n", ellapsed);
    for j=1:numel(TrackList)
        fprintf("Track %d: %f\n", TrackList{j}.TrackID, TrackList{j}.ProbOfExist);
    end
    disp("");
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        % Plot data
        %cla(ax(1));
        for i = 1:numel(plots)
            delete(plots(i))
        end
        plots = [];

        % NOTE: if your image is RGB, you should use flipdim(img, 1) instead of flipud.
        hold on;
        [lat,lon,h] = ned2geodetic(DataList{k}(2,:),DataList{k}(1,:),0,50.346069,-4.113670,0,referenceEllipsoid('wgs84'));
        h2 = plot(ax(1), lon,lat,'k*','MarkerSize', 10);
        %h2 = plot(ax(1), DataList{k}(1,:),DataList{k}(2,:),'k*','MarkerSize', 10);
        plots(end+1) = h2;
        for j=1:numel(TrackList)
            %h2 = plot(ax(1), TrackList{j}.Filter.StateMean(1),TrackList{j}.Filter.StateMean(3),'.','LineWidth',1);
            [lat,lon,h] = ned2geodetic(TrackList{j}.Trajectory(3,:),TrackList{j}.Trajectory(1,:),0,50.346069,-4.113670,0,referenceEllipsoid('wgs84'));
            h2 = plot(ax(1), lon,lat,'r','LineWidth',3);
            %h2 = plot(ax(1), TrackList{j}.Trajectory(1,:),TrackList{j}.Trajectory(3,:),'r','LineWidth',1);
            plots(end+1) = h2;
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            %[lat,lon,h] = ned2geodetic(North,East,0,50.346069,-4.113670,0,referenceEllipsoid('wgs84'));
            %h2 = plot_gaussian_ellipsoid(TrackList{j}.Filter.StateMean([1 3]), TrackList{j}.Filter.StateCovar([1 3],[1 3]),1,20,ax(1));
            %plots(end+1) = h2;
%             plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-5,int2str(TrackList{j}.TrackID));
%             plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-130,num2str(TrackList{j}.ProbOfExist,2));
%             plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-5,int2str(TrackList{j}.TrackID));
%             plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-130,num2str(TrackList{j}.ProbOfExist,2));
            %set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        %axis(ax(1),V_bounds)
        F(k) = getframe(ax(1));
        % Plot PHD
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.Particles([1,3],:)');
        %contour3(X,Y,density,50);
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        %set(h, 'edgecolor','none')
        hold on;
        plot(ax(2), myphd.Particles(1,:), myphd.Particles(3,:), '.')
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