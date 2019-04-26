%% BernoulliParticleFilterX Demo 
% ----------------------
% * This script demonstrates the process of configuring an running a
%   BernoulliParticleFilterX object to perform single-target state estimation.
%
% * A toy single-target scenario is considered, for a target that generates
%   regular reports of it's position with the possibility of clutter and
%   missed detections
%
%% Extract the GroundTruth data from the example workspace
load('single-target-tracking.mat');
NumIter = size(TrueTrack.Trajectory,2);

%% Models

lambdaV = 5; % Expected number of clutter measurements over entire surveillance region
V = 10^2;     % Volume of surveillance region (10x10 2D-grid)
V_bounds = [0 10 0 10]; % [x_min x_max y_min y_max]
P_D = 0.8; 

% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate an Observation model
measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);

% Instantiate a clutter model
clutter_model = PoissonRateUniformPositionX('ClutterRate',lambdaV,'Limits',[V_bounds(1:2);V_bounds(3:4)]);

% Instantiate birth model
birth_model = DistributionBasedBirthModelX('Distribution', UniformDistributionX([V_bounds(1:2); ...
                                                                                [-0.1 0.1 ];...
                                                                                V_bounds(3:4);...
                                                                                [-0.1 0.1 ]]),...
                                           'BirthIntensity', 0.00005);
% Instantiate detection model                                       
detection_model = ConstantDetectionProbabilityX('DetectionProbability',P_D);
                                       

% Compile the State-Space model
model = StateSpaceModelX(transition_model,...
                       measurement_model,...
                       'Clutter',clutter_model,...
                       'Birth', birth_model,...
                       'Detection', detection_model);

%% Simulation
% Data Simulator
dataSim = SingleTargetMeasurementSimulatorX(model);

% Simulate some measurements from ground-truth data
MeasurementScans = dataSim.simulate(TrueTrack);
measurements = [MeasurementScans.Vectors];

%% Initiation

% Setup prior assuming we have some intuition of the target's initial
% position
priorParticles = TrueTrack.Trajectory(1).Vector ...
                 + model.Measurement.finv(model.Measurement.random(4000));
priorWeights = ones(1,4000)/4000;
config.Model = model;
config.StatePrior = ParticleStateX(priorParticles,priorWeights./sum(priorWeights));
config.StatePrior.Metadata.ExistenceProbability = 0.5;
config.BirthScheme = {'Expansion', 5000};
config.SurvivalProbability = 0.99;

% Initiate a track using the generated prior
track = TrackX(config.StatePrior, TagX(1));

%% Estimation           
% Instantiate a filter objects
filter = BernoulliParticleFilterX(config);                            
figure;
for t = 2:NumIter
    
    % Provide filter with the new measurement
    MeasurementList = MeasurementScans(t);
    
    % Perform filtering
    filter.MeasurementList = MeasurementList;
    filter.predict();
    filter.update();
    
    % Log the data
    track.Trajectory(end+1) = filter.StatePosterior;
    
    clf;
    hold on;
    measurements = [MeasurementScans(1:t).Vectors];
    meas = measurement_model.finv(measurements);
    true_means = [TrueTrack.Trajectory(1:t).Vector];
    track_means = [track.Trajectory(1:t).Mean];
    plot(true_means(1,1:t), true_means(3,1:t),'.-k', track_means(1,1:t), track_means(3,1:t), 'b-', meas(1,:), meas(3,:), 'rx');
    plot_gaussian_ellipsoid(track.Trajectory(t).Mean([1,3],1), track.Trajectory(t).Covar([1,3],[1,3]));
    legend('GroundTrouth','Estimated Mean','Measurements', 'Estimated Covariance');
    xlabel("x coordinate (m)");
    ylabel("y coordinate (m)");
    axis([2 9 1 9]);
    drawnow();
end    