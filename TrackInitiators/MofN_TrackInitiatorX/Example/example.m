% This is a test script which demonstrates the usage of the "MofN_TrackInitiatorX" class.
% =========================================================================>

% Load the ground truth data
load('multiple-robot-tracking.mat');

% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
numTrueTracks = 3;

% Model parameter shortcuts
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]
P_D = 1;    % Probability of detection
timestep_duration = duration(0,0,1);

%% Models
transition_model = ConstantVelocityX('VelocityErrVariance', 0.0001,...
                                     'NumDims', 2,...
                                     'TimestepDuration', timestep_duration);
measurement_model = LinearGaussianX('NumMeasDims', 2,...
                                    'NumStateDims', 4,...
                                    'MeasurementErrVariance', 0.02,...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[V_bounds(1:2);...
                                                      V_bounds(3:4)]);
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model, 'Detection', detection_model);


%% Generate DataList
meas_simulator = MultiTargetMeasurementSimulatorX('Model',model);
% meas_simulator.DetectionProbability = 1;
DataList = meas_simulator.simulate(GroundTruthStateSequence);
N = numel(DataList);

%% Base Filter
obs_covar= measurement_model.covar();
PriorState = GaussianStateX(zeros(4,1), transition_model.covar() + blkdiag(obs_covar(1,1), 0, obs_covar(2,2),0));
base_filter = KalmanFilterX('Model', model, 'StatePrior', PriorState);

%% Data Associator
config.ClutterModel = clutter_model;
config.Clusterer = NaiveClustererX();
config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
config.DetectionProbability = P_D;
jpdaf = JointProbabilisticDataAssocX(config);

%% Track Initiator

% Initiate Data Associator
config_pdaf.ClutterModel = clutter_model;
% config.Clusterer = NaiveClustererX();
config_pdaf.Gater = EllipsoidalGaterX(2,'GateLevel',5);
config_pdaf.DetectionProbability = P_D;
pdaf = ProbabilisticDataAssocX(config_pdaf);

% Initiate Tag Generator
tag_gen = RandSampleTagGeneratorX(1:10000);

% Prepare initiator parameters
config_ti.TagGenerator = tag_gen;
config_ti.InitFilter = base_filter;
config_ti.DataAssociator = pdaf;
config_ti.ConfirmThreshold = [8,10];
config_ti.DeleteThreshold = [5,10];
CovarThreshold = 4*PriorState.Covar;
config_ti.CustomDeleteConditionFcn = ...
    @(x,t) t.Filter.StatePosterior.Covar(1,1)>CovarThreshold(1,1) ...
           || t.Filter.StatePosterior.Covar(3,3)>CovarThreshold(3,3);

% Create the track initiator
myti = MofN_TrackInitiatorX(config_ti);

%% START OF SIMULATION
%  ===================>

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

TrackList = [];
for k=2:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    
    %% Extract DataList at time k
    MeasurementList = DataList(k);
    timestamp_km1 = DataList(k-1).Timestamp;
    timestamp_k = MeasurementList.Timestamp;
    dt = timestamp_k - timestamp_km1;
    transition_model.TimestepDuration = dt;
    fprintf('Timestamp = %s\n================>\n',timestamp_k);
    
    %% Process JPDAF
    jpdaf.MeasurementList = MeasurementList;
    jpdaf.TrackList = TrackList;
    jpdaf.predictTracks();
    jpdaf.associate();    
    jpdaf.updateTracks();
    
    %% Perform Track initiation
    [TrackList, TentativeTrackList] = myti.initiateTracks(jpdaf.TrackList, MeasurementList, jpdaf.AssocWeightsMatrix);
        
    %% Plot update step results
    if(ShowPlots)
            
        cla(ax(1));
        %imagesc(ax(1),[min_x max_x], [min_y max_y], flipud(img));
        hold on;
        if(exist('data_plot','var'))
            delete(data_plot);
        end
        data_inv = measurement_model.finv(MeasurementList.Vectors);
        data_plot = plot(ax(1), data_inv(1,:), data_inv(3,:),'k*','MarkerSize', 10);

        % Plot confirmed tracks
        for j=1:numel(TrackList)
            means = [TrackList{j}.Trajectory.Mean];
            h2 = plot(ax(1), means(1,:),means(3,:),'-','LineWidth',1);
            h2 = plotgaussellipse(TrackList{j}.Filter.StatePosterior.Mean([1 3]),...
                                  TrackList{j}.Filter.StatePosterior.Covar([1 3],[1 3]),...
                                  'Color','r',...
                                  'Axis',ax(1)); 
        end

        % Plot tentative tracks
        for j=1:numel(TentativeTrackList)
            h2 = plot(ax(1), TentativeTrackList{j}.Filter.StatePosterior.Mean(1),TentativeTrackList{j}.Filter.StatePosterior.Mean(3),'.','LineWidth',1);
            if j==2
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            end
            h2 = plotgaussellipse(TentativeTrackList{j}.Filter.StatePosterior.Mean([1 3]),...
                                  TentativeTrackList{j}.Filter.StatePosterior.Covar([1 3],[1 3]),...
                                  'Color','g',...
                                  'Axis',ax(1));
        end
        
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Robot positions (Update)');
        title(ax(1),str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1),V_bounds)
        pause(0.01)
    end
end