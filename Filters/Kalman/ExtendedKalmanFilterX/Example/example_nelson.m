
clear F
% Instantiate a Transitionamic model
dyn = ConstantVelocityX('NumDims',2,'VelocityErrVariance',2^2, 'TimestepDuration', duration(0,0,10));

% Instantiate an Observation model
%obs = LinearGaussianX('NumMeasDims',2,'NumStateDims',4,'MeasurementErrVariance',0.02,'Mapping',[1 3]);
obs = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[(pi/100)^2, 50^2],'Mapping',[1 3]);

% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Extract the ground truth data from the example workspace
load('nelson.mat');
NumIter = size(measurements,2);

% Now let's estimate!!

% Use the first measurement as our prior mean and the measurement noise
% plus process noise covariance as our prior covariance
xPrior = [-1.6741e+03, 0, 875.1847, 0]';
PPrior = diag([1000, 10, 1000, 10]);

% Setup prior
StatePrior = GaussianStateX(xPrior,PPrior);

% Instantiate a filter object
filter = ExtendedKalmanFilterX('Model',ssm, 'StatePrior', StatePrior);
%filter = UnscentedKalmanFilterX('Model',ssm, 'StatePrior', StatePrior);

Log.Estimates.StateMean = zeros(ssm.Transition.NumStateDims,NumIter);
Log.Estimates.StateCovar = zeros(ssm.Transition.NumStateDims,...
                                      ssm.Transition.NumStateDims,NumIter);
                        
                        
figure;
for t = 1:NumIter
    
    % Provide filter with the new measurement
    filter.MeasurementList = measurements(:,t);
    filter.Model.Transition.TimestepDuration = dt(t);
    
    % Perform filtering
    filter.predict();
    filter.update();
    
    % Log the data
    Log.Estimates.StateMean(:,t) = filter.StatePosterior.Mean;
    Log.Estimates.StateCovar(:,:,t) = filter.StatePosterior.Covar;
    
    clf;
    hold on;
    meas = obs.finv(measurements(:,1:t));
    plot(Log.Estimates.StateMean(1,1:t), Log.Estimates.StateMean(3,1:t), 'b-', meas(1,:), meas(3,:), 'rx');
    plot_gaussian_ellipsoid(Log.Estimates.StateMean([1,3],t), Log.Estimates.StateCovar([1,3],[1,3],t));
    legend('Estimated Mean','Measurements', 'Estimated Covariance');
    xlabel("x coordinate (m)");
    ylabel("y coordinate (m)");
    %axis([-1800 -400 -800 1000]);
    axis([-1673.31483408488,-427.506070065913,-103.492188459434,937.889546649686]);
    drawnow();
%     F(t) = getframe(gcf);
end    

%F = F(2:end);
% vidObj = VideoWriter(sprintf('nelson_radar_ukf.avi'));
% vidObj.Quality = 100;
% vidObj.FrameRate = 5;
% open(vidObj);
% writeVideo(vidObj, F);
% close(vidObj);