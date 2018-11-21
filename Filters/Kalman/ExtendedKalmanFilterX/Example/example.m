
% Instantiate a Transitionamic model
dyn = ConstantVelocityX('NumDims',2,'VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
%obs = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[0.001,0.02],'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Extract the ground truth data from the example workspace
load('example.mat');
NumIter = size(truth,2);

% Add faux velocity components to the data
truth = [truth(1,:);zeros(1,NumIter);truth(2,:);zeros(1,NumIter)];

% Simulate some measurements from ground-truth data
measurements = ssm.Measurement.feval(truth,true);

% Now let's estimate!!

% Use the first measurement as our prior mean and the measurement noise
% plus process noise covariance as our prior covariance
xPrior = truth(:,1);
PPrior = 10*dyn.covar();

% Setup prior
StatePrior = GaussianStateX(xPrior,PPrior);

% Instantiate a filter object
filter = ExtendedKalmanFilterX('Model',ssm, 'StatePrior', StatePrior);

Log.Estimates.StateMean = zeros(ssm.Transition.NumStateDims,NumIter);
Log.Estimates.StateCovar = zeros(ssm.Transition.NumStateDims,...
                                      ssm.Transition.NumStateDims,NumIter);
                        
                        
figure;
for t = 1:NumIter
    
    % Provide filter with the new measurement
    filter.MeasurementList = measurements(:,t);
    
    % Perform filtering
    filter.predict();
    filter.update();
    
    % Log the data
    Log.Estimates.StateMean(:,t) = filter.StatePosterior.Mean;
    Log.Estimates.StateCovar(:,:,t) = filter.StatePosterior.Covar;
    
    clf;
    hold on;
    meas = obs.finv(measurements(:,1:t));
    plot(truth(1,1:t),truth(3,1:t),'.-k', Log.Estimates.StateMean(1,1:t), Log.Estimates.StateMean(3,1:t), 'b-', meas(1,:), meas(3,:), 'rx');
    plot_gaussian_ellipsoid(Log.Estimates.StateMean([1,3],t), Log.Estimates.StateCovar([1,3],[1,3],t));
    legend('GroundTrouth','Estimated Mean','Measurements', 'Estimated Covariance');
    xlabel("x coordinate (m)");
    ylabel("y coordinate (m)");
    axis([2 9 1 9]);
    drawnow();
end    