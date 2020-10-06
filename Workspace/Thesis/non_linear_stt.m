% This is a test/example script which demonstrates the usage of the 
% ProbabilisticDataAssociationX class.
% =========================================================================>

%% Load the ground truth data
load('maneuvering_robot.mat');
NumIter = size(TrueTrack.Trajectory,2);

%% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots

%% Model parameter shortcuts
lambdaV = 0; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 25 0 15]; % [x_min x_max y_min y_max]
P_D = 1;    % Probability of detection
timestep_duration = duration(0,0,1);

%% Models
transition_model = ConstantVelocityX('VelocityErrVariance', 0.01^2,...
                                     'NumDims', 2,...
                                     'TimestepDuration', timestep_duration);
measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',...
                                            [(pi/180)^2,0.1^2],'Mapping',[1 3]);
% measurement_model = LinearGaussianX('NumMeasDims', 2,...
%                                     'NumStateDims', 4,...
%                                     'MeasurementErrVariance', 0.2^2,...
%                                     'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[V_bounds(1:2);...
                                                      V_bounds(3:4)]);   
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model, 'Detection', detection_model);

%% Simulation
% Data Simulator
dataSim = SingleTargetMeasurementSimulatorX(model);

% Simulate some measurements from ground-truth data
MeasurementScans = dataSim.simulate(TrueTrack);
measurements = [MeasurementScans.Vectors];

%% Base Filter
% base_filter = ExtendedKalmanFilterX('Model', model);
base_filter = UnscentedKalmanFilterX('Model', model);
% base_filter = ParticleFilterX('Model', model);

%% Initiation
% Use the first measurement scan to perform single-point initiation
measurement = MeasurementScans(1).Measurements;
timestamp = measurement.Timestamp;
NumTracks = 1;
TrackList = cell(1,NumTracks);

% Setup prior
% xPrior = measurement.Model.Measurement.finv(measurement.Vector);
true_mean = TrueTrack.Trajectory(1).Vector;
xPrior = [true_mean(1);0.2;true_mean(3); 0.2];
PPrior = 30*transition_model.covar();
dist = GaussianDistributionX(xPrior,PPrior);
StatePrior = ParticleStateX(dist,5000,timestamp);

% Initiate a track using the generated prior
track = TrackX(StatePrior, TagX(1));
TrackList{1} = track;
TrackList{1}.addprop('Filter');
TrackList{1}.Filter = copy(base_filter);
TrackList{1}.Filter.initialise('Model',model,'StatePrior',StatePrior);

%% Data Associator
config.ClutterModel = clutter_model;
config.Clusterer = NaiveClustererX();
config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
config.DetectionModel = detection_model;
assocFilter = ProbabilisticDataAssocX(config);
assocFilter.TrackList = TrackList;

%% Metric Generator
rmse_gen = RMSEX();

%% START OF SIMULATION
%  ===================>

% Create figure windows
if(ShowPlots)
    figure('units','normalized','outerposition',[0 0 .5 1])
end

for k=2:NumIter
    fprintf('Iteration = %d/%d\n================>\n',k,NumIter);
    
    %% Extract DataList at time k
    MeasurementList = MeasurementScans(k);
    
    %% Process JPDAF
    assocFilter.MeasurementList = MeasurementList;
    assocFilter.TrackList = TrackList;
    assocFilter.predictTracks();
    assocFilter.associate();    
    assocFilter.updateTracks();
        
     %% Plot update step results
    if(ShowPlots)
            
        ax(1) = gca;
        cla(ax(1));
        hold on;
        if(exist('data_plot','var'))
            delete(data_plot);
        end
        if MeasurementList.NumMeasurements ~= 0
            data_inv = measurement_model.finv(MeasurementList.Vectors);
            data_plot = plot(ax(1), data_inv(1,:), data_inv(3,:),'k*','MarkerSize', 10);
        end
        
        % Plot tracks
        t_means = [TrueTrack.Trajectory.Vector];
        plot(ax(1),t_means(1,1:k),t_means(3,1:k),'k.-','LineWidth',1);
        for j=1:numel(assocFilter.TrackList)
            means = [assocFilter.TrackList{j}.Trajectory.Mean];
            h2 = plot(ax(1), means(1,:),means(3,:),'.-','LineWidth',1);
            h2 = plotgaussellipse(assocFilter.TrackList{j}.Filter.StatePosterior.Mean([1 3]),...
                                  assocFilter.TrackList{j}.Filter.StatePosterior.Covar([1 3],[1 3]),...
                                  'Color','r',...
                                  'Axis',ax(1)); 
        end      
        % set the y-axis back to normal.
        str = sprintf('Target positions (Update)');
        set(ax(1), 'ydir','normal');
        title(str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis(ax(1), V_bounds)
        pause(0.01);
        
    end
end
rmse = rmse_gen.evaluate(TrueTrack, track)
