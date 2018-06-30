% Script used to generate an example of a multimodal distribution in the 
% case of (J)PDAF
% =========================================================================

%% Initialise a few things

% Instantiate a state prior
prior.mean = [2000;0;50;0];
prior.covar = diag([100; 100; 15; pi/45]);

% Instantiate a Dynamic model
%dyn = ConstantVelocityModelX_2D('VelocityErrVariance',0.0001);
dyn = ConstantHeadingModelX('VelocityErrVariance',(7)^2, 'HeadingErrVariance',0.20);


% Instantiate an Observation model
% obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.1,'Mapping',[1 3]);
obs = Polar2CartGaussModelX('NumStateDims',4,'RangeErrVariance',(10)^2,...
                          'ThetaErrVariance',(pi/180)^2,'Mapping',[1 2]);
obs2 = Polar2CartGaussModelX('NumStateDims',4,'RangeErrVariance',(50)^2,...
                          'ThetaErrVariance',(pi/90)^2,'Mapping',[1 2]);
                      
% Compile the State-Space model
ssm = StateSpaceModelX(dyn,obs);

% Initiate PDAF parameters
Params_jpdaf.Clusterer = NaiveClustererX();
Params_jpdaf.Gater = EllipsoidalGaterX(2,'GateLevel',20)';
Params_jpdaf.ProbOfDetect = 0.9;

% Instantiate a TrackList with a single track object
TrackList{1} = TrackX();
TrackList{1}.addprop('Filter');

% Generate DataList
NumMeasurements = 10;
pred.mean = dyn.feval(prior.mean);
%measurements = obs2.heval(repmat(pred.mean,1,NumMeasurements),true);
meas2state = obs.heval_inv(measurements);

%% EKF
% Allocate a Filter to the track
% Params_pf.Model = ssm;
% Params_pf.NumParticles = 5000000;
% Params_pf.PriorParticles = mvnrnd(prior.mean', prior.covar, Params_pf.NumParticles)';
% Params_pf.PriorWeights = ones(1,Params_pf.NumParticles)./Params_pf.NumParticles;
% TrackList{1}.Filter = ParticleFilterX(Params_pf);

% Allocate a Filter to the track
Params_kf.PriorStateMean = prior.mean;
Params_kf.PriorStateCovar = prior.covar; %blkdiag(POmodel.Params.R(1)/2, 2^2, 2*pi);%CVmodel.Params.Q(1);
Params_kf.Model = ssm;
TrackList{1}.Filter = ExtendedKalmanFilterX(Params_kf);
%TrackList{i}.Filter = UnscentedKalmanFilterX(Params_kf);

% Instantiate JPDAF
jpdaf = ProbabilisticDataAssocX(Params_jpdaf);
jpdaf.TrackList = TrackList;
jpdaf.MeasurementList = measurements; 

% Predict
jpdaf.TrackList{1}.Filter.predict();

% Associate and update
jpdaf.associate();    
jpdaf.updateTracks();

% Plot distribution
figure;
subplot(1,3,1)
%Xk = [Params_pf.PriorParticles(1:2,:)];
Xk = mvnrnd(jpdaf.TrackList{1}.Filter.StateMean',jpdaf.TrackList{1}.Filter.StateCovar,5000000)';
Xk = Xk(1:2,:);
[bandwidth,density1,X1,Y1]=kde2d(Xk([1,2],:)');
hold on;
%plot3(meas2state(1,:),meas2state(2,:),repmat(-10^6,1,NumMeasurements),'r+', 'MarkerSize', 2);
surf(X,Y,density);
xlim([1900 2200])
ylim([-150 150]);
xlabel('X (m)');
ylabel('Y (m)');
shading interp
contour(X1,Y1,density1,4,'r','LineWidth',2)
set(gca,'Ydir','reverse')
title('JPDA-EKF');
view(0,-90);

%% UKF
% Allocate a Filter to the track
Params_kf.PriorStateMean = prior.mean;
Params_kf.PriorStateCovar = prior.covar; %blkdiag(POmodel.Params.R(1)/2, 2^2, 2*pi);%CVmodel.Params.Q(1);
Params_kf.Model = ssm;
TrackList{1}.Filter = UnscentedKalmanFilterX(Params_kf);

% Instantiate JPDAF
jpdaf = ProbabilisticDataAssocX(Params_jpdaf);
jpdaf.TrackList = TrackList;
jpdaf.MeasurementList = measurements; 

% Predict
jpdaf.TrackList{1}.Filter.predict();

% Associate and update
jpdaf.associate();    
jpdaf.updateTracks();

% Plot distribution
subplot(1,3,2)
%Xk = [Params_pf.PriorParticles(1:2,:)];
Xk = mvnrnd(jpdaf.TrackList{1}.Filter.StateMean',jpdaf.TrackList{1}.Filter.StateCovar,5000000)';
Xk = Xk(1:2,:);
[bandwidth,density1,X1,Y1]=kde2d(Xk([1,2],:)');
hold on;
%plot3(meas2state(1,:),meas2state(2,:),repmat(-10^6,1,NumMeasurements),'r+', 'MarkerSize', 2);
surf(X,Y,density);
xlim([1900 2200])
ylim([-150 150]);
xlabel('X (m)');
ylabel('Y (m)');
shading interp
contour(X1,Y1,density1,4,'r','LineWidth',2)
set(gca,'Ydir','reverse')
title('JPDA-UKF');
view(0,-90);

%% PF

% Allocate a Filter to the track
Params_pf.Model = ssm;
Params_pf.NumParticles = 5000000;
Params_pf.PriorParticles = mvnrnd(prior.mean', prior.covar, Params_pf.NumParticles)';
Params_pf.PriorWeights = ones(1,Params_pf.NumParticles)./Params_pf.NumParticles;
TrackList{1}.Filter = ParticleFilterX(Params_pf);

% Instantiate JPDAF
jpdaf = ProbabilisticDataAssocX(Params_jpdaf);
jpdaf.TrackList = TrackList;
jpdaf.MeasurementList = measurements; 

% Predict
jpdaf.TrackList{1}.Filter.predict();

% Associate and update
jpdaf.associate();    
jpdaf.updateTracks();

subplot(1,3,3)
% clf;
Xk = [jpdaf.TrackList{1}.Filter.Particles(1:2,:)];
[bandwidth,density,X,Y]=kde2d(Xk([1,2],:)');
hold on;
%plot3(meas2state(1,:),meas2state(2,:),repmat(-10^6,1,NumMeasurements),'r+', 'MarkerSize', 2);
surf(X,Y,density);
shading interp
contour(X,Y,density,4,'r','LineWidth',2)
set(gca,'Ydir','reverse')
xlabel('X (m)');
ylabel('Y (m)');
xlim([1900 2200])
ylim([-150 150]);
title('JPDA-PF');
view(0,-90);
%
