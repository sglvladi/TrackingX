% This is a test script which demonstrates the usage of the 
% JointProbabilisticDataAssociation class.
% =========================================================================>

% Load the ground truth data
load('multiple-robot-tracking.mat');

% Plot settings
ShowPlots = 1;              % Set to 0 to hide plots
NumTracks = 3;

% Model parameter shortcuts
lambdaV = 10; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]
P_D = 0.8;    % Probability of detection

%% Models
transition_model = ConstantVelocityX('VelocityErrVariance', 0.01^2,...
                                     'NumDims', 2,...
                                     'TimestepDuration', duration(0,0,1));
measurement_model = LinearGaussianX('NumMeasDims', 2,...
                                    'NumStateDims', 4,...
                                    'MeasurementErrVariance', 0.2^2,...
                                    'Mapping', [1 3]);
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,...
                                            'Limits',[V_bounds(1:2);...
                                                      V_bounds(3:4)]);
% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model);

%% Generate DataList
meas_simulator = MultiTargetMeasurementSimulatorX('Model',model);
% meas_simulator.DetectionProbability = 1;
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
config.DetectionProbability = P_D;
pdaf = JointIntegratedProbabilisticDataAssocX(config);

%% Initiate TrackList
TrackList = cell(1,NumTracks);
for i=1:NumTracks
    xPrior = [GroundTruth{1}(1,i); 0; GroundTruth{1}(2,i); 0];
    PPrior = 10*transition_model.covar();
    StatePrior = GaussianStateX(xPrior,PPrior);
    TrackList{i} = TrackX(StatePrior);
    TrackList{i}.addprop('Filter');
    TrackList{i}.Filter = copy(base_filter);
    TrackList{i}.Filter.initialise('Model',model,'StatePrior',StatePrior);
    TrackList{i}.addprop('ExistenceProbability');
    TrackList{i}.ExistenceProbability = 0.5;
end

pdaf.TrackList = TrackList;

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
    pdaf.MeasurementList = MeasurementList;
    pdaf.TrackList = TrackList;
    pdaf.predictTracks();
    pdaf.associate();    
    pdaf.updateTracks();
    
    for j=1:NumTracks
        fprintf("%f, ",pdaf.TrackList{j}.ExistenceProbability);
    end
    
    %% Update target trajectories
    for t = 1: numel(pdaf.TrackList)
        pdaf.TrackList{t}.Trajectory(end+1) = pdaf.TrackList{t}.Filter.StatePosterior;
    end

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

        % Plot tracks
        for j=1:numel(pdaf.TrackList)
            means = [pdaf.TrackList{j}.Trajectory.Mean];
            h2 = plot(ax(1), means(1,:),means(3,:),'-','LineWidth',1);
            h2 = plot_gaussian_ellipsoid(pdaf.TrackList{j}.Filter.StatePosterior.Mean([1 3]), pdaf.TrackList{j}.Filter.StatePosterior.Covar([1 3],[1 3]),'r',1,20,ax(1)); 
        end        
        % set the y-axis back to normal.
        set(gca,'ydir','normal');
        str = sprintf('Estimated state x_{1,k} vs. x_{2,k}');
        title(str)
        xlabel('X position (m)')
        ylabel('Y position (m)')
        axis([0 10 0 10])
        pause(0.01)
    end
end

    
ospa_vals= zeros(N,3);
ospa_c= 1;
ospa_p= 1;
for k=1:N
    trueX = [x_true(k,:);y_true(k,:)];
    estX = zeros(2,0);
    for i=1:NumTracks
        if(Logs{i}.Estimates.ExistenceProbability(k)>0.2)
            estX = [estX, Logs{i}.Estimates.StateMean([1 3],k)];
        end
    end
    [ospa_vals(k,1), ospa_vals(k,2), ospa_vals(k,3)]= OSPAX.ospa_dist(trueX,estX,ospa_c,ospa_p);
end
ospa = mean(ospa_vals,1);
figure
subplot(2,2,[1 2]), plot(1:k,ospa_vals(1:k,1));
subplot(2,2,3), plot(1:k,ospa_vals(1:k,2));
subplot(2,2,4), plot(1:k,ospa_vals(1:k,3));