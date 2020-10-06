% example_radar.m
% ====================================================>
% The following seem to be descent parameters:
%   - ProbOfDetection = 0.8
%   - ProbOfDeath = 0.2
%   - VelocityErrVariance = 2
%   - ObsErrVariance = 50
%   - ProbOfConfirm = 0.9
%   - PHD.BirthScheme = {'Expansion', 5000}
%

%% Radar specific settings
RadarName = 'staddon';
switch(RadarName)
    case('staddon')
        % Surveillance region parameters
        V_bounds = [-8154.72944624983;... % X-min | 
                    -212.289393440959;... % X-max | Surveillance region bounding
                    -7548.44272179096;... % Y-min | box coordinates (m)
                    4355.32645897434]';   % Y-max |
%         RadarCoords.lat = 50.33933333;
%         RadarCoords.lon = -4.12527778;
        RadarCoords.lat = 50.346069;
        RadarCoords.lon = -4.113670;
    case('longroom')
        % Surveillance region parameters
        V_bounds = [-8154.72944624983;... % X-min | 
                    3000;...              % X-max | Surveillance region bounding
                    -7548.44272179096;... % Y-min | box coordinates (m)
                    4355.32645897434]';   % Y-max |
        RadarCoords.lat = 50.36286670714617;
        RadarCoords.lon = -4.156833300366998;
end

% Load dataset
load(strcat(RadarName,'.mat'));
N = size(DataList,2); % Simulation length

%% Focus on region of high clutter
lon_lim = [-4.195 -4.11];
lat_lim = [50.31 50.375];
% lon_lim = [-4.1538   -4.1152];
% lat_lim = [50.3097   50.3373];
% lon_lim = [-4.15 -4.13];
% lat_lim = [50.335 50.35];
% Convert LLA limits to NED
[y_lim, x_lim] = geodetic2ned(lat_lim,...
                              lon_lim, ...
                              0,...
                              RadarCoords.lat,...
                              RadarCoords.lon,...
                              0,...
                              referenceEllipsoid('wgs84'));
V_bounds = [x_lim(1);... % X-min | 
            x_lim(2);... % X-max | Surveillance region bounding
            y_lim(1);... % Y-min | box coordinates (m)
            y_lim(2)]';  % Y-max |                         
for i = 1:N
    j = 1;
    while j <= size(DataList{i},2)
        if(DataList{i}(3,j)>lat_lim(2) || DataList{i}(3,j)<lat_lim(1) || DataList{i}(4,j)<lon_lim(1) || DataList{i}(4,j)>lon_lim(2))
            DataList{i}(:,j) = [];
            j = j -1;
        end
        j = j + 1;
    end
end

radar_meas_convert;
%% Plot & Recording settings
% Plot settings
ShowPlots = 1;              % Set to 0 to prevent showing any plots
ShowUpdate = 1;             % Set to 0 to skip showing update plots
ShowTrackInfo = 1;
NumPersistFrames = 50;           

% Recording settings
Record = 1;                 % Set to (0|1) to turn video recording (off|on)
FrameRate = 0.5;            % Number of frames per second
VideoQuality = 100;         % Set to desired quality percentage
VideoPathName = strcat(RadarName,'_heavy_clutter_staddon.avi'); % Set to the desired path and name of produced recording

% Model parameter shortcuts
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3))); % Total area of surveillance region
P_D = 0.7;    % Probability of detection
timestep_duration = duration(0,0,2);

%% Instantiation of necessary components

% Instantiate a Dynamic model (CV model with q = 1 m/s^2)
transition_model = ConstantVelocityX('VelocityErrVariance', 0.1,...
                                     'NumDims', 2,...
                                     'TimestepDuration', timestep_duration);
% Instantiate an Observation model (Variance of 50m^2 on each coordinate)
measurement_model = LinearGaussianX('NumMeasDims', 2,...
                                    'NumStateDims', 4,...
                                    'MeasurementErrVariance', 10^2,...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[V_bounds(1:2);...
                                                      V_bounds(3:4)]);
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-5 5 ];...
                                                                                V_bounds(3:4);...
                                                                                [-5 5 ]]),...
                                           'BirthIntensity', 0.0001);                                                  
% Compile the State-Space model
model = StateSpaceModelX(transition_model, measurement_model,...
                        'Clutter',clutter_model,...
                        'Detection', detection_model,...
                        'Birth', birth_model);

%% Base Filter
obs_covar= measurement_model.covar();
PriorState = GaussianStateX(zeros(4,1), transition_model.covar() + blkdiag(obs_covar(1,1), 0, obs_covar(2,2),0));
base_filter = KalmanFilterX('Model', model, 'StatePrior', PriorState);

%% Data Associator
config.ClutterModel = clutter_model;
config.Clusterer = NaiveClustererX();
config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
config.DetectionModel = detection_model;
jpdaf = JointIntegratedProbabilisticDataAssocX(config);

%% Track Initiator

% Initiate Data Associator
config_phd.Model = model;
[priorParticles, priorWeights] = model.Birth.random(50000);
config_phd.StatePrior = ParticleStateX(priorParticles,10*priorWeights);
config_phd.BirthScheme = {'Expansion', 5000};
config_phd.SurvivalProbability = 0.995;

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config_phd);

% Initiate Tag Generator
tag_gen = UuidTagGeneratorX();

% Prepare initiator parameters
config_ti.TagGenerator = tag_gen;
config_ti.InitFilter = base_filter;
config_ti.PhdFilter = myphd;
config_ti.ConfirmThreshold = 0.8;

% Create the track initiator
myti = PhdTrackInitiatorX(config_ti);

%% Track Deleter
config_td = struct('Fieldname', 'ExistenceProbability',...
                   'ReferenceValue', 0.1,...
                   'ReferenceOperand', 'lt');
mytd = FieldBasedDeleterX(config_td);

TrackList = [];

%% Create plot windows
if(ShowPlots)
    
    % Map plot
    figure('units','normalized','outerposition',[0 0 .5 1])
    ax(1) = gca;
    plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
    axis(ax(1),[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
%     axis(ax(1),[-4.195 -4.11 50.31 50.375])
    
%     % PHD Intensity plot
%     figure('units','normalized','outerposition',[.5 0 .5 1])
%     ax(2) = gca;
    
    plots = [];
end

%% START OF SIMULATION
% ===================>
for k=2:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);

    %% Extract DataList at time k
    MeasurementList = MeasurementScans(k);
    timestamp_km1 = MeasurementScans(k-1).Timestamp;
    timestamp_k = MeasurementList.Timestamp;
    dt = timestamp_k - timestamp_km1;
    transition_model.TimestepDuration = dt;
    fprintf('Timestamp = %s\n================>\n',timestamp_k);
    
    %% Process JPDAF
    jpdaf.MeasurementList = MeasurementList;
    jpdaf.TrackList = TrackList;
    jpdaf.predictTracks();
    jpdaf.associate(TrackList, MeasurementList);    
    jpdaf.updateTracks();
    
    % Delete all plots (other than map)
    for i = 1:numel(plots)
        delete(plots(i))
    end
    plots = [];
    hold on;

%     if(MeasurementList.NumMeasurements>0)
%         % Convert measurements to LLA and plot them
%         [lat,lon,~] = ned2geodetic(MeasurementList.Vectors(2,:),...
%                                    MeasurementList.Vectors(1,:),...
%                                    0,...
%                                    RadarCoords.lat,...
%                                    RadarCoords.lon,...
%                                    0,...
%                                    referenceEllipsoid('wgs84'));
%         plots(end+1) = plot(ax(1), lon,lat,'y*','MarkerSize', 10);
%         
%         for i=1:numel(lat)
%             plots(end+1) = text(ax(1), lon(:,i)+0.0001,lat(:,i)+0.00027,num2str(i),'FontSize',8,'Color','w');
%         end
% %             plot(ax(1), RadarCoords.lon,RadarCoords.lat,...
% %                                '-s','MarkerSize',20,...
% %                                'MarkerEdgeColor','red',...
% %                                'MarkerFaceColor',[1 .6 .6]);
%     end
%     pause(0.01);
     %% Perform Track initiation
    [TrackList] = myti.initiateTracks(jpdaf.TrackList, MeasurementList, jpdaf.AssocWeightsMatrix);
        
    %% Perform Track deletion
    TrackList = mytd.deleteTracks(TrackList);
    
    % Plot update step results
    if(ShowPlots && ShowUpdate)
        
        % Delete all plots (other than map)
        for i = 1:numel(plots)
            delete(plots(i))
        end
        plots = [];
        hold on;
        
        if(MeasurementList.NumMeasurements>0)
            % Convert measurements to LLA and plot them
            [lat,lon,~] = ned2geodetic(MeasurementList.Vectors(2,:),...
                                       MeasurementList.Vectors(1,:),...
                                       0,...
                                       RadarCoords.lat,...
                                       RadarCoords.lon,...
                                       0,...
                                       referenceEllipsoid('wgs84'));
            plots(end+1) = plot(ax(1), lon,lat,'y*','MarkerSize', 10);
%             plot(ax(1), RadarCoords.lon,RadarCoords.lat,...
%                                '-s','MarkerSize',20,...
%                                'MarkerEdgeColor','red',...
%                                'MarkerFaceColor',[1 .6 .6]);

        end
        
        % Plot all existing tracks
        for j=1:numel(TrackList)
            track = TrackList{j};
            % Convert track trajectory to LLA and plot it
            means = [track.Trajectory.Mean];
            [lat,lon,~] = ned2geodetic(means(3,:),...
                                       means(1,:),...
                                       0,...
                                       RadarCoords.lat,...
                                       RadarCoords.lon,...
                                       0,...
                                       referenceEllipsoid('wgs84'));
            traj_length = size(lon,2);
            if(traj_length>NumPersistFrames)
                start = traj_length-NumPersistFrames;
            else
                start = 1;
            end
            plots(end+1) = plot(ax(1), lon(:,start:end),lat(:,start:end),'-.w','LineWidth',2);
            plots(end+1) = plot(ax(1), lon(:,end),lat(:,end),'ws','MarkerSize',15);
            
            % Convert track velocity to LLA and plot it
            [lat_vel,lon_vel,~] = ned2geodetic(track.State.Mean(4,end),...
                                               track.State.Mean(2,end),...
                                               0,...
                                               lat(:,end),...
                                               lon(:,end),...
                                               0,...
                                               referenceEllipsoid('wgs84'));
            lat_vel = lat_vel-lat(:,end);
            lon_vel = lon_vel-lon(:,end);
            
            plots(end+1) = quiver(ax(1), lon(:,end),lat(:,end),20*lon_vel,20*lat_vel,'r','LineWidth',1.5);
            speed_kmph = sqrt(track.State.Mean(4,end)^2+track.State.Mean(2,end)^2)*3.6;
            speed_knot = speed_kmph/1.852;
            
            if(ShowTrackInfo)
                plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)+0.00027,strcat("Sog:",num2str(speed_knot,2)," kt"),'FontSize',8,'Color','w');
                plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)-0.00027,strcat("PoE:",num2str(track.ExistenceProbability*100,3)," %"),'FontSize',8,'Color','w');
            end
            % TODO: Convert heading and state covariance to LLA and plot them
            % [lat,lon,h] = ned2geodetic(North,East,0,50.346069,-4.113670,0,referenceEllipsoid('wgs84'));
            % h2 = plot_gaussian_ellipsoid(TrackList{j}.Filter.StateMean([1 3]), TrackList{j}.Filter.StateCovar([1 3],[1 3]),1,20,ax(1));
            % plots(end+1) = h2;
            % plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-5,int2str(TrackList{j}.TrackID));
            % plots(end+1) = text(ax(1),TrackList{j}.Filter.StateMean(1)+20,TrackList{j}.Filter.StateMean(3)-130,num2str(TrackList{j}.ProbOfExist,2));
        end
        
        % Add axis labels
        xlabel('Longitude')
        ylabel('Latitude')
                
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