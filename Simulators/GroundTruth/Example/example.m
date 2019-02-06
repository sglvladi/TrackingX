%% Models
% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate an Observation model
measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
%measurement_model = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

% Compile the State-Space model
model = StateSpaceModelX(transition_model,measurement_model);

%% Groud-Truth Simulator
num_timesteps = 100;
timestep_duration = duration(0,0,5);
initial_state = GroundTruthStateX([0; 0; 0; 0], datetime());
gnd_sim = SingleTargetGroundTruthSimulatorX('Model',model, ...
                                            'InitialState', initial_state,...
                                            'NumTimesteps', num_timesteps,...
                                            'TimestepDuration', timestep_duration);
                                        
%% Simulate Ground-Truth
track = gnd_sim.simulate();

%% Measurement simulator
det_sim = SingleTargetMeasurementSimulatorX(model);

%% Simulate measurements
measurements = det_sim.simulate(track);

%% Plot the output
plotter = BasicPlotterX(model);

plotter.figure();
plotter.plotObject(track);
plotter.plotObject(measurements);