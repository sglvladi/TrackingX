% This is a test script which demonstrates the usage of the "PhdTrackInitiatorX" class.
% =========================================================================>

% Load the ground truth data
load('multi-target-tracking-1.mat');
GroundTruthStateSequence = GroundTruthStateSequence(~cellfun('isempty',GroundTruthStateSequence));

% Plot settings
ShowPlots = 0;              % Set to 0 to hide plots
numTrueTracks = 3;

% Model parameter shortcuts
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = 5000^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [-2500 2500 -2500 2500]; % [x_min x_max y_min y_max]
P_D = 0.9;    % Probability of detection
timestep_duration = duration(0,0,1);

%% Models
transition_model = ConstantVelocityX('VelocityErrVariance', 5^2,...
                                     'NumDims', 2,...
                                     'TimestepDuration', timestep_duration);
% measurement_model = LinearGaussianX('NumMeasDims', 2,...
%                                     'NumStateDims', 4,...
%                                     'MeasurementErrVariance', 10^2,...
%                                     'Mapping', [1 3]);
% clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
%                                             'Limits',[V_bounds(1:2);...
%                                                       V_bounds(3:4)]);
measurement_model = RangeBearing2CartesianX('NumStateDims', 4,...
                                    'MeasurementErrVariance', [(pi/180)^2, 10^2],...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[[-pi, pi];...
                                                      [0, 2500]]);
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-10 10 ];...
                                                                                V_bounds(3:4);...
                                                                                [-10 10 ]]),...
                                           'BirthIntensity', 0.0001);

% Compile the State-Space model
model = StateSpaceModelX(transition_model, measurement_model,...
                        'Clutter',clutter_model,...
                        'Detection', detection_model,...
                        'Birth', birth_model);


%% Simulate measurements, based on ground-truth
meas_simulator = MultiTargetMeasurementSimulatorX('Model',model);
DataList = meas_simulator.simulate(GroundTruthStateSequence);
N = numel(DataList);

%% Base Filter
obs_covar= measurement_model.covar();
PriorState = GaussianStateX(zeros(4,1), transition_model.covar() + blkdiag(obs_covar(1,1), 0, obs_covar(2,2),0));
base_filter = UnscentedKalmanFilterX('Model', model, 'StatePrior', PriorState);
% base_filter = ParticleFilterX('Model', model, 'StatePrior', PriorState);

%% Data Associator
config.ClutterModel = clutter_model;
config.DetectionModel = detection_model;
config.Clusterer = NaiveClustererX();
config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
config.DetectionProbability = P_D;
jpdaf = JointIntegratedProbabilisticDataAssocX(config);

%% Track Initiator

% Initiate Data Associator
config_phd.Model = model;
[priorParticles, priorWeights] = model.Birth.random(50000);
config_phd.StatePrior = ParticleStateX(priorParticles,10*priorWeights);
config_phd.BirthScheme = {'Expansion', 5000};
config_phd.SurvivalProbability = 0.99;

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

%% Metric Generators
ospa = OSPAX('CutOffThreshold',100,'Order',2);
gospa = GOSPAX('CutOffThreshold',100,'Order',2);

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
gospa_vals= zeros(N,4);
timestamp_km1 = DataList(1).Timestamp;
for k=2:N
    fprintf('Iteration = %d/%d\n================>\n',k,N);
    
    %% Extract DataList at time k
    MeasurementList = DataList(k);
    if(MeasurementList.NumMeasurements == 0)
        continue;
    end
    timestamp_km1 = DataList(k-1).Timestamp;
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
    
    %% Perform Track initiation
    [TrackList] = myti.initiateTracks(jpdaf.TrackList, MeasurementList, jpdaf.AssocWeightsMatrix);
        
    %% Perform Track deletion
    TrackList = mytd.deleteTracks(TrackList);
    
    %% Evaluate metric
    [gospa_vals(k,1), gospa_vals(k,2), gospa_vals(k,3), gospa_vals(k,4)]= ...
        gospa.evaluate(GroundTruthStateSequence{k},jpdaf.TrackList);
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
figure
title("GOSPA");
subplot(2,3,[1 2 3]), plot(1:k,gospa_vals(1:k,1)), title("GOSPA Metric");
subplot(2,3,4), plot(1:k,gospa_vals(1:k,2)), title("GOSPA Localisation");
subplot(2,3,5), plot(1:k,gospa_vals(1:k,3)), title("GOSPA Missed");
subplot(2,3,6), plot(1:k,gospa_vals(1:k,4)), title("GOSPA False");