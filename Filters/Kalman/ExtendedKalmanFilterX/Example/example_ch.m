
% Instantiate a Dynamic model
dyn = ConstantHeadingModelX('VelocityErrVariance',0.0001, 'HeadingErrVariance',0.09);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.04,'Mapping',[1 2]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Instantiate a Kalman Filter object
ekf = ExtendedKalmanFilterX(ssm);

% Extract the ground truth data from the example workspace
load('example.mat');
NumIter = size(truth,2);

% Add faux velocity components to the data
truth = [truth(1,:);truth(2,:);zeros(1,NumIter);zeros(1,NumIter)];

% Simulate some measurements from ground-truth data
%measurements = ssm.Obs.heval(truth,true);

% Now let's estimate!!

% Use the first measurement as our prior mean and the measurement noise
% plus process noise covariance as our prior covariance
xPrior = [measurements(1,1); measurements(2,1); 0; 0];
measErrCov = ssm.Obs.covariance();
stateErrCov = ssm.Dyn.covariance();
PPrior = stateErrCov + blkdiag(measErrCov(1,1),measErrCov(2,2),0,0);
ekf.initialise('PriorStateMean',xPrior,'PriorStateCovar',PPrior); 

Log.Estimates.StateMean = zeros(ssm.Dyn.NumStateDims,NumIter);
Log.Estimates.StateCovar = zeros(ssm.Dyn.NumStateDims,...
                                      ssm.Dyn.NumStateDims,NumIter);
                        
figure;
for t = 1:NumIter
    
    % Provide KalmanFilter with the new measurement
    ekf.Measurement = measurements(:,t);
    
    % Perform filtering
    ekf.predict();
    ekf.update();
    
    % Log the data
    Log.Estimates.StateMean(:,t) = ekf.StateMean;
    Log.Estimates.StateCovar(:,:,t) = ekf.StateCovar;
    
    clf;
    hold on;
    plot(truth(1,1:t),truth(2,1:t),'.-k', Log.Estimates.StateMean(1,1:t), Log.Estimates.StateMean(2,1:t), 'b-', measurements(1,1:t), measurements(2,1:t), 'rx');
    plot_gaussian_ellipsoid(Log.Estimates.StateMean([1,2],t), Log.Estimates.StateCovar([1,2],[1,2],t));
    legend('GroundTrouth','Estimated Mean','Measurements', 'Estimated Covariance');
    xlabel("x coordinate (m)");
    ylabel("y coordinate (m)");
    axis([2 9 1 9]);
    pause(0.1);
end


    