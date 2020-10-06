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
RadarName = 'staddon_clutter';

RadarCoords.lat = 50.346069;
RadarCoords.lon = -4.113670;

% Load dataset
load(strcat(RadarName,'.mat'));
N = size(MeasurementScans,2); % Simulation length

%% Focus on region of high clutter

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

%% Plot & Recording settings
% Plot settings
ShowPlots = 0;              % Set to 0 to prevent showing any plots
ShowUpdate = 1;             % Set to 0 to skip showing update plots
ShowTrackInfo = 0;
NumPersistFrames = 10000;
PlotLatLon = 1;

% Recording settings
Record = 0;                 % Set to (0|1) to turn video recording (off|on)
FrameRate = 0.5;            % Number of frames per second
VideoQuality = 100;         % Set to desired quality percentage
VideoPathName = strcat(RadarName,'_smcphd_pf.avi'); % Set to the desired path and name of produced recording

% Model parameter shortcuts
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = (abs(V_bounds(2)-V_bounds(1))*abs(V_bounds(4)-V_bounds(3))); % Total area of surveillance region
P_D = 0.5;    % Probability of detection
timestep_duration = duration(0,0,2);

%% Instantiation of necessary components

% Instantiate a Dynamic model (CV model with q = 1 m/s^2)
transition_model = OrnsteinUhlenbeckModelX('NumDims',2,...
                             'VelocityErrVariance',0.1,...
                             'DampingCoefficient',0.01,...
                             'TimestepDuration',timestep_duration);
% Instantiate an Observation model (Variance of 50m^2 on each coordinate)
measurement_model = LinearGaussianX('NumMeasDims', 2,...
                                    'NumStateDims', 4,...
                                    'MeasurementErrVariance', 10^2,...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPosition2X('ClutterRate',lambdaV,...
                                            'Limits',[V_bounds(1:2);...
                                                      V_bounds(3:4)]);
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-5 5 ];...
                                                                                V_bounds(3:4);...
                                                                                [-5 5 ]]),...
                                           'BirthIntensity', 0.1);
                                       
% Compile the State-Space model
model = StateSpaceModelX(transition_model, measurement_model,...
                        'Clutter',clutter_model,...
                        'Detection', detection_model,...
                        'Birth', birth_model);

%% Base Filter
obs_covar= measurement_model.covar();
PriorState = GaussianStateX(zeros(4,1), transition_model.covar() + blkdiag(obs_covar(1,1), 0, obs_covar(2,2), 0));
base_filter = ExtendedKalmanFilterX('Model', model, 'StatePrior', PriorState);
% base_filter = ParticleFilterX('Model', model, 'StatePrior', PriorState);

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
config_phd.StatePrior = ParticleStateX(priorParticles, priorWeights);
config_phd.BirthScheme = {'Expansion', 5000};
config_phd.SurvivalProbability = 0.99;

% Instantiate PHD filter
myphd = SMC_PHDFilterX(config_phd);
% Initiate Tag Generator
% tag_gen = UuidTagGeneratorX();
tag_gen = RandSampleTagGeneratorX(1:100000);

% Prepare initiator parameters
config_ti.TagGenerator = tag_gen;
config_ti.InitFilter = base_filter;
config_ti.PhdFilter = myphd;
config_ti.ConfirmThreshold = 0.8;

% Create the track initiator
myti = PhdTrackInitiatorX(config_ti);

%% Track Deleter
config_td = struct('Fieldname', 'ExistenceProbability',...
                   'ReferenceValue', 0.05,...
                   'ReferenceOperand', 'lt');
mytd = FieldBasedDeleterX(config_td);

TrackList = {};
DeletedTracks = {};
%% Create plot windows
if(ShowPlots)
    
    % Map plot
    figure('units','normalized','outerposition',[0 0 .5 0.9])
    ax(1) = gca;
    if PlotLatLon
        plot_google_map('Axis',ax(1),'APIKey','AIzaSyBXKujdtXRZiqya1soVS9pxBzYR4g7aGvM','Resize',3,'Scale',2,'MapType','satellite');
        axis(ax(1),[lon_lim(1) lon_lim(2) lat_lim(1) lat_lim(2)])
    end
    
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
%     for i = 1:numel(plots)
%         delete(plots(i))
%     end
%     plots = [];
%     hold on;

     %% Perform Track initiation
    [TrackList] = myti.initiateTracks(jpdaf.TrackList, MeasurementList, jpdaf.AssocWeightsMatrix);
        
    %% Perform Track deletion
    [TrackList, del_tracks] = mytd.deleteTracks(TrackList);
    DeletedTracks = [DeletedTracks, del_tracks];
    
    del_tracks = [];
    for t=1:numel(TrackList)
        track = TrackList{t};
        xk = track.State.Mean;
        [ylim, xlim] = geodetic2ned(lat_lim,...
                                lon_lim, ...
                                0,...
                                RadarCoords.lat,...
                                RadarCoords.lon,...
                                0,...
                                referenceEllipsoid('wgs84'));
        if xk(1,:)>x_lim(2) || xk(1,:)<x_lim(1) || xk(3,:)>y_lim(2) || xk(3,:)<y_lim(1)
            del_tracks{end+1} = track;
            TrackList{t} = [];
        end
    end
    TrackList = TrackList(~cellfun('isempty',TrackList));
    DeletedTracks = [DeletedTracks, del_tracks];
    
    %% Plot update step results
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
        end
        
        % Plot region of low PD
        x = [-3422.85261310741,-3191.58441229017,-2993.55608122595,-3228.08558459284, -3422.85261310741];
        y = [1.472966334316646e+03,1.123217741935625e+03,1.243952901078573e+03,1.579667558371468e+03, 1.472966334316646e+03];
        [y,x,~] = ned2geodetic(y,...
                               x,...
                               0,...
                               RadarCoords.lat,...
                               RadarCoords.lon,...
                               0,...
                               referenceEllipsoid('wgs84'));
        plots(end+1) = plot(ax(1), x, y, 'm-');

        % Plot region of high clutter
        x = [-469.710602736631,-2572.92295458513,-2516.31234127506,-1320.74163107899,-835.307363807233, -469.710602736631];
        y = [-4915.41679013513,-4065.50122858976,-2858.60918927224,-1380.19032924249,-1261.08679763585, -4915.41679013513];
%         x = [-2858.62556879160,-2856.97056906902,-108.925126527319,-108.988225201227,-2858.62556879160];
%         y = [-4044.74840882411,-974.655039469575,-975.424226884065,-4045.51804181679,-4044.74840882411];
        [y,x,~] = ned2geodetic(y,...
                               x,...
                               0,...
                               RadarCoords.lat,...
                               RadarCoords.lon,...
                               0,...
                               referenceEllipsoid('wgs84'));
        plots(end+1) = plot(ax(1), x, y, 'r-');
            
        % Plot all existing tracks
        for j=1:numel(TrackList)
            track = TrackList{j};
            % Convert track trajectory to LLA and plot it
            means = [track.Trajectory.Mean];
            if PlotLatLon
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
                if ShowTrackInfo
                    speed_kmph = sqrt(track.State.Mean(4,end)^2+track.State.Mean(2,end)^2)*3.6;
                    speed_knot = speed_kmph/1.852;
                    plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)+0.00027,strcat("Sog:",num2str(speed_knot,2)," kt"),'FontSize',8,'Color','w');
                    plots(end+1) = text(ax(1), lon(:,end)+0.001,lat(:,end)-0.00027,strcat("PoE:",num2str(track.ExistenceProbability*100,3)," %"),'FontSize',8,'Color','w');
                end
            else
                traj_length = size(lon,2);
                if(traj_length>NumPersistFrames)
                    start = traj_length-NumPersistFrames;
                else
                    start = 1;
                end
                plots(end+1) = plot(ax(1), means(1,start:end),means(3,start:end),'-.k','LineWidth',2);
                plots(end+1) = plot(ax(1), means(1,end),means(3,end),'ks','MarkerSize',15);
                
                x_vel = track.State.Mean(3,end)*cos(track.State.Mean(4,end));
                y_vel = track.State.Mean(3,end)*sin(track.State.Mean(4,end));
                plots(end+1) = quiver(ax(1), means(1,end),means(2,end),20*x_vel,20*y_vel,'r','LineWidth',1.5);
                
                if ShowTrackInfo
                    speed_kmph = track.State.Mean(3,end)*3.6;
                    speed_knot = speed_kmph/1.852;
                    plots(end+1) = text(ax(1), means(1,end)+60,means(2,end)+50,strcat("Sog:",num2str(speed_knot,2)," kt"),'FontSize',8,'Color','k');
                    plots(end+1) = text(ax(1), means(1,end)+60,means(2,end)-50,strcat("PoE:",num2str(track.ExistenceProbability*100,3)," %"),'FontSize',8,'Color','k');
                end
            end
        end
        
        % Add axis labels
%         xlabel('Longitude')
%         ylabel('Latitude')
        title(datestr(timestamp_k));
        a = xlim;
        plots(end+1) = text(ax(1), a(1), 50.3522, datestr(timestamp_k), 'FontSize',15,'Color','w');
%         set(gca,'visible','off')
        pause(0.0001)
        % Store video frame
        if(Record)
            F(k) = getframe(ax(1));
        end
    end

end

% Create video file and write to it
if(Record)
    F = F(2:end);
    vidObj = VideoWriter(char(VideoPathName));
    vidObj.Quality = VideoQuality;
    vidObj.FrameRate = FrameRate;
    open(vidObj);
    writeVideo(vidObj, F);
    close(vidObj);
end