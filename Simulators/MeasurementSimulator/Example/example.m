% Instantiate a Transitionamic model
transition_model = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate a Measurement model
measurement_model = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
%obs = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

% Instantiate a clutter model
clutter_model = PoissonRateUniformPositionX('ClutterRate',5,'Limits',[0 10; 0 10]);

% Compile the State-Space model
ssm = StateSpaceModelX(transition_model,measurement_model,'Clutter',clutter_model);

% Extract the ground truth data from the example workspace
load('example.mat');

% Simulate measurements
DataList = MeasurementSimulatorX(ssm,GroundTruth);

figure;
for i=1:numel(DataList)
    clf;
    hold on;
    measurements = ssm.Measurement.finv(DataList{i});
    plot(measurements(1,:),measurements(3,:),'r*');
    plot(GroundTruth{i}(1,:),GroundTruth{i}(2,:), 'ob');
    axis([0 10 0 10]);
    drawnow;
end