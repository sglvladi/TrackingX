
% Instantiate a Dynamic model
dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);

% Instantiate an Observation model
obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',1,'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Instantiate a Particle Filter object
pf = ParticleFilterX(ssm);
%ukf = UnscentedKalmanFilterX(ssm);

% Extract the ground truth data from the example workspace
load('example.mat');
NumIter = size(truth,2);

% Add faux velocity components to the data
truth = [truth(1,:);zeros(1,NumIter);truth(2,:);zeros(1,NumIter)];

% Simulate some measurements from ground-truth data
measurements = ssm.Obs.heval(truth,true);

% Now let's estimate!!

% Use the first measurement as out prior mean and the measurement noise
% plus process noise covariance as our prior covariance
xPrior = [measurements(1,1); 0; measurements(2,1); 0];
measErrCov = ssm.Obs.covariance();
stateErrCov = ssm.Dyn.covariance();
PPrior = stateErrCov + blkdiag(measErrCov(1,1),0,measErrCov(2,2),0);
pf.initialise('PriorDistFcn',@(N)deal(mvnrnd(xPrior(:,ones(1,N))',PPrior)',1/N*ones(1,N)));
%ukf.initialise('PriorStateMean',xPrior,'PriorStateCov',PPrior); 

Log.Estimates.StateMean = zeros(ssm.Dyn.NumStateDims,NumIter);
Log.Estimates.StateCovar = zeros(ssm.Dyn.NumStateDims,...
                                      ssm.Dyn.NumStateDims,NumIter);
                        
figure;
for t = 1:NumIter
    
    % Provide KalmanFilter with the new measurement
    pf.Measurement = measurements(:,t);
    
    % Perform filtering
    pf.predict();
    pf.update();
    
    % Log the data
    Log.Estimates.StateMean(:,t) = pf.StateMean;
    Log.Estimates.StateCovar(:,:,t) = pf.StateCovar;
    
    clf;
    hold on;
    plot(truth(1,1:t),truth(3,1:t),'.-k', Log.Estimates.StateMean(1,1:t), Log.Estimates.StateMean(3,1:t), 'b-o', measurements(1,1:t), measurements(2,1:t), 'rx');
    plot_gaussian_ellipsoid(Log.Estimates.StateMean([1,3],t), Log.Estimates.StateCovar([1,3],[1,3],t));
    legend('GroundTrouth','Estimated Mean','Measurements', 'Estimated Covariance');
    xlabel("x coordinate (m)");
    ylabel("y coordinate (m)");
    axis([2 9 1 9]);
    pause(0.1);
end


    