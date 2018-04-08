% example_radar.m
% ====================================================>
% The following seem to be descent parameters:
%   - ProbOfDetection = 0.8
%   - ProbOfDeath = 0.005
%   - VelocityErrVariance = 2
%   - ObsErrVariance = 50
%   - ProbOfConfirm = 0.9
%   - PHD.BirthScheme = {'Expansion', 5000}
%
%% Clear all data and import dataset
clear all;
load('radar_ned.mat');
N = size(DataList,2); % Simulation length

%% Plot & Recording settings
% Plot settings
ShowPlots = 1;              % Set to 0 to prevent showing any plots
ShowUpdate = 1;             % Set to 0 to skip showing update plots

% Recording settings
Record = 1;                 % Set to (0|1) to turn video recording (off|on)
FrameRate = 0.5;            % Number of frames per second
VideoQuality = 100;         % Set to desired quality percentage
VideoPathName = "./example_radar.avi"; % Set to the desired path and name of produced recording

%% Instantiation of necessary components

% Instantiate a Dynamic model (CV model with q = 1 m/s^2)
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',2);
dyn.TimeVariant = 2;

% Instantiate an Observation model (Variance of 50m^2 on each coordinate)
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',50,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Surveillance region parameters
V_bounds = [-8154.72944624983;... % X-min | 
            -212.289393440959;... % X-max | Surveillance region bounding
            -7548.44272179096;... % Y-min | box coordinates (m)
            4355.32645897434]';   % Y-max |
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3))); % Total area of surveillance region

% Assign PHD parameter values
config.NumParticles = 50000;
config.Model = ssm;
config.BirthIntFcn = @(Np)[abs(V_bounds(2)-V_bounds(1))*rand(1,Np)+V_bounds(1);... % Uniform position across the surveillance 
                           9*rand(1,Np);...                                      % region.                         
                           abs(V_bounds(4)-V_bounds(3))*rand(1,Np)+V_bounds(3);... % Uniform speed between 0-9 m/s in both
                           9*rand(1,Np)];                                        % X and Y.
config.PriorDistFcn = @ (Np) deal(config.BirthIntFcn(Np), repmat(1/Np, Np, 1)');   % Uniform position and weights.
config.BirthScheme = {'Expansion', 5000};
config.ProbOfDeath = 0.005;                                                        % Probability of death = 0.5%
config.ProbOfDetection = 0.8;                                                      % Probability of detection = 70%
config.ResamplingScheme = 'Multinomial';                                           % Use Multinomial Resampling

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config);

% Instantiate PF filter
mypf = ParticleFilterX(ssm);

% Initiate PDAF parameters
Params_pdaf.Clusterer = NaiveClustererX();
Params_pdaf.Gater = EllipsoidalGaterX(2,'ProbOfGating',0.99)';
Params_pdaf.ProbOfDetect = 0.8;
mypdaf = JointProbabilisticDataAssocX(Params_pdaf);

% Initiate Track Initiator
config_ti.Filter = mypf;
config_ti.PHDFilter = myphd;
config_ti.ProbOfGating = 0.99;
config_ti.ProbOfConfirm = 0.9;
myti = PhdExistProbTrackInitiatorX(config_ti);

TrackList = [];

%% Create plot windows
if(ShowPlots)
    
    % Map plot
    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
    plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2);
    axis([-4.19499354141717 -4.10940387037447 50.3123450958395 50.3735587491040])
    
    % PHD Intensity plot
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;
    
    plots = [];
end

%% START OF SIMULATION
% ===================>
for k=1:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);

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
    
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        
        % Delete all plots (other than map)
        for i = 1:numel(plots)
            delete(plots(i))
        end
        plots = [];
        hold on;
        
        % Convert measurements to LLA and plot them
        [lat,lon,~] = ned2geodetic(DataList{k}(2,:),...
                                   DataList{k}(1,:),...
                                   0,50.346069,...
                                   -4.113670,...
                                   0,...
                                   referenceEllipsoid('wgs84'));
        plots(end+1) = plot(ax(1), lon,lat,'k*','MarkerSize', 10);
        
        % Plot all existing tracks
        for j=1:numel(TrackList)
            
            % Convert track trajectory to LLA and plot it
            [lat,lon,~] = ned2geodetic(TrackList{j}.Trajectory(3,:),...
                                       TrackList{j}.Trajectory(1,:),...
                                       0,50.346069,...
                                       -4.113670,...
                                       0,...
                                       referenceEllipsoid('wgs84'));
            plots(end+1) = plot(ax(1), lon,lat,'r','LineWidth',2);
            
            % TODO: Convert heading and state covariance to LLA and plot them
            % [lat,lon,h] = ned2geodetic(North,East,0,50.346069,-4.113670,0,referenceEllipsoid('wgs84'));
            % h2 = plot_gaussian_ellipsoid(TrackList{j}.Filter.StateMean([1 3]), TrackList{j}.Filter.StateCovar([1 3],[1 3]),1,20,ax(1));
            % plots(end+1) = h2;
            % plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-5,int2str(TrackList{j}.TrackID));
            % plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-130,num2str(TrackList{j}.ProbOfExist,2));
        end
        
        % Add axis labels
        xlabel('X position (m)')
        ylabel('Y position (m)')
        
        % Plot PHD Intensity
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.Particles([1,3],:)');
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        hold on;
        plot(ax(2), myphd.Particles(1,:), myphd.Particles(3,:), '.',...
                    myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds]);
        str = sprintf('PHD intensity (Update)');
        xlabel(ax(2),'X position (m)')
        ylabel(ax(2),'Y position (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
        
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