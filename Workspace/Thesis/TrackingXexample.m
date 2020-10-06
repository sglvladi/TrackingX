% This is a test script which demonstrates the usage of the "PhdTrackInitiatorX" class.
% =========================================================================>

% Load the ground truth data
load('multi-target-tracking-1.mat');
GroundTruthStateSequence = GroundTruthStateSequence(~cellfun('isempty',GroundTruthStateSequence));

% Plot settings
ShowPlots = 0;              % Set to 0 to hide plots
numTrueTracks = 3;

% Model parameter shortcuts
lambdaV = 50; % Expected number of clutter measurements over entire surveillance region
V = 5000^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [-2500 2500 -2500 2500]; % [x_min x_max y_min y_max]
P_D = 0.9;    % Probability of detection
timestep_duration = duration(0,0,1);

%% Models
transition_model = ConstantVelocityX('VelocityErrVariance', 5^2,...
                                     'NumDims', 2,...
                                     'TimestepDuration', timestep_duration);
measurement_model = RangeBearing2CartesianX('NumStateDims', 4,...
                                    'MeasurementErrVariance', [(pi/270)^2, 10^2],...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[[-pi, pi];...
                                                      [0, 2500]]);
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);

birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-10 10 ];...
                                                                                V_bounds(3:4);...
                                                                                [-10 10 ]]),...
                                           'BirthIntensity', 0.1);

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
base_filter = ParticleFilterX('Model', model, 'StatePrior', PriorState);

%% Data Associator
config.ClutterModel = clutter_model;
config.Clusterer = NaiveClustererX();
config.Gater = EllipsoidalGaterX(2,'GateLevel',10)';
config.DetectionModel = detection_model;
jpdaf = JointIntegratedProbabilisticDataAssocX(config);

%% Track Initiator
% Initiate SMC PHD Filter
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
config_ti.ConfirmThreshold = 0.6;

% Create the track initiator
track_initiator = PhdTrackInitiatorX(config_ti);

%% Track Deleter
config_td = struct('Fieldname', 'ExistenceProbability',...
                   'ReferenceValue', 0.1,...
                   'ReferenceOperand', 'lt');
track_deleter = FieldBasedDeleterX(config_td);

%% Metric Generators
gospa = GOSPAX('CutOffThreshold',200,'Order',1);
gospa_vals= zeros(N,4);

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
    
    % PHD Intensity plot
    figure('units','normalized','outerposition',[.5 0 .5 1])
    ax(2) = gca;
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
    jpdaf.associate(TrackList, MeasurementList);    
    jpdaf.updateTracks();
    
    %% Perform Track initiation
    [TrackList] = track_initiator.initiateTracks(jpdaf.TrackList, MeasurementList, jpdaf.AssocWeightsMatrix);
        
    %% Perform Track deletion
    TrackList = track_deleter.deleteTracks(TrackList);
    
    %% Evaluate performance metrics
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
        for j=1:numel(GroundTruthTracks)
            track = GroundTruthTracks(j);
            if track.TimeOfInitiation>timestamp_k || track.TimeOfLastUpdate<timestamp_k
                continue
            end
            means = [];
            for i = 1: numel(track.Trajectory)
                state=track.Trajectory(i);
                if state.Timestamp > timestamp_k
                    break
                end
                means(:,end+1) = state.Vector;
            end
            h1 = plot(ax(1), means(1,:),means(3,:),'b--','LineWidth',1);
            if j~=numel(GroundTruthTracks)
                set(get(get(h1,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
        end
        
        % Plot confirmed tracks
        for j=1:numel(TrackList)
            means = [TrackList{j}.Trajectory.Mean];
            h2 = plot(ax(1), means(1,:),means(3,:),'r-','LineWidth',1);
            if j~=numel(TrackList)
                set(get(get(h2,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
            end
            h3 = plotgaussellipse(TrackList{j}.Filter.StatePosterior.Mean([1 3]),... 
                                  TrackList{j}.Filter.StatePosterior.Covar([1 3],[1 3]),...
                                  'Color','r',...
                                  'Axis',ax(1));
            set(get(get(h3,'Annotation'),'LegendInformation'),'IconDisplayStyle','off'); % Exclude line from legend
        end
        
        % set the y-axis back to normal.
        set(ax(1),'ydir','normal');
        str = sprintf('Track positions');
        title(ax(1),str)
        xlabel(ax(1),'X (m)')
        ylabel(ax(1),'Y (m)')
        %h_legend = legend([h1 h4 h3 h2],'Ground Truth', 'measurements', 'PDAF-EKF', 'PDAF-PF');%, 'PF-PCH', 'PF-PCHR');   
        legend(ax(1), 'Measurements', 'GroundTruth Tracks','Estimated Tracks');
%         set(h_legend, 'FontSize', 18);
        axis(ax(1),V_bounds)
        
        % Plot PHD Intensity
        cla(ax(2), 'reset');
        [bandwidth,density,X,Y]=kde2d(myphd.StatePosterior.Particles([1,3],:)');
        h = surf(ax(2),X,Y,density);        
        shading interp
        colormap(ax(2), jet(3000))
        hold on;
%         plot(ax(2), myphd.StatePosterior.Particles(1,:), myphd.StatePosterior.Particles(3,:), '.',...
%                     myphd.MeasurementList(1,:), myphd.MeasurementList(2,:), 'y*');
        axis(ax(2), [V_bounds]);
        str = sprintf('PHD intensity');
        xlabel(ax(2),'X (m)')
        ylabel(ax(2),'Y (m)')
        zlabel(ax(2),'Intensity')
        title(ax(2),str)
        pause(0.01)
        
        pause(0.01)
    end
end

figure
title("GOSPA");
subplot(2,3,[1 2 3]), plot(1:k,gospa_vals(1:k,1)), title("GOSPA Metric");
subplot(2,3,4), plot(1:k,gospa_vals(1:k,2)), title("GOSPA Localisation");
subplot(2,3,5), plot(1:k,gospa_vals(1:k,3)), title("GOSPA Missed");
subplot(2,3,6), plot(1:k,gospa_vals(1:k,4)), title("GOSPA False");