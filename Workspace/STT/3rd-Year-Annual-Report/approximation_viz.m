% Script used to generate a visualisation of how well EKF, UKF and PF
% approximate a distribution
% =========================================================================

% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [2000; 0.1; 0; 0];

% Instantiate an Observation model
% obs = LinGaussObsModelX_2D('NumStateDims',4,'ObsErrVariance',0.2,'Mapping',[1 2]);
obs = RangeBearing2CartesianX('NumStateDims',4,'MeasurementErrVariance',[(pi/90)^2,(1)^2],...
                              'Mapping',[1 3]);

% View the transition matrix and process covariance matrices
R = obs.covar();

% Predict the target's position and velocity after the interval has passed
measurement  = obs.feval(state);

%% EKF
[ekf_mean, H_jac] = ExtendedKalmanFilterX.computeJac_(@(x)obs.finv(x),measurement);
ekf_covar = H_jac*R*H_jac';

%% UKF
alpha = 0.5;
kappa = 0;
beta  = 2;
% Calculate unscented transformation parameters
[c, Wmean, Wcov, OOM] = matlabshared.tracking.internal.calcUTParameters(alpha,beta,kappa,2);
% Form the sigma points
X = formSigmaPoints(measurement, R, c);
% Perform Unscented Transform to get predicted measurement mean,
% covariance and cross-covariance
[ukf_mean,ukf_covar,Pxy] = unscentedTransform(@(x)obs.finv(x),X,Wmean,Wcov,OOM);

%% PF
% Generate 50 random noise samples from the dynamic model
noise = obs.random(50000000);

% Add noise to the measurements
Yk1 = noise + measurement;
Xk1 = obs.finv(Yk1);

Xk = [Xk1];%,Xk2];
Yk = [Yk1];%,Yk2];

figure;
[bandwidth,density,X,Y]=kde2d(Xk([1,3],:)');
xLim = [1985,2005];
yLim = [-200,200];

% EKF
subplot(1,3,1);
%contour(X,Y,density,1);
h = surf(X,Y,density);
z = get(h,'ZData');
% set(h,'ZData',z-10^6)
hold on;
shading interp
x1 = xLim(1):0.1:xLim(2); x2 = yLim(1):0.1:yLim(2);
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],ekf_mean([1,3])',ekf_covar([1,3],[1,3]));
F = reshape(F,length(x2),length(x1));
contour(x1,x2,F,4,'r','LineWidth',2);
% plot_gaussian_ellipsoid(ekf_mean([1,3]),ekf_covar([1,3],[1,3]),'r',1);
% plot_gaussian_ellipsoid(ekf_mean([1,3]),ekf_covar([1,3],[1,3]),'r',2);
% shading interp
%colormap(jet(3000))
title('Taylor-Series Expansion (EKF)')
xlabel('X (m)');
ylabel('Y (m)');
xlim(xLim);
ylim(yLim);
view(0,-90);

% UKF
subplot(1,3,2);
h = surf(X,Y,density);
z = get(h,'ZData');
%set(h,'ZData',z-10^6)
hold on;
shading interp
x1 = xLim(1):0.1:xLim(2); x2 = yLim(1):0.1:yLim(2);
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],ukf_mean([1,3])',ukf_covar([1,3],[1,3]));
F = reshape(F,length(x2),length(x1));
contour(x1,x2,F,4,'r','LineWidth',2);
%plot_gaussian_ellipsoid(ukf_mean([1,3]),ukf_covar([1,3],[1,3]),'r',1);
%plot_gaussian_ellipsoid(ukf_mean([1,3]),ukf_covar([1,3],[1,3]),'r',2);
%colormap(jet(3000))
title('Unscented Transform (UKF)')
xlabel('X (m)');
ylabel('Y (m)');
xlim(xLim);
ylim(yLim);
view(0,-90);


% PF
subplot(1,3,3);
%contour(X,Y,density,1);
h = surf(X,Y,density);
z = get(h,'ZData');
%set(h,'ZData',z-10^6)
hold on;    
shading interp
%plot(Xk(1,1:1000),Xk(3,1:1000),'r.');
%[bandwidth,density,X,Y]=kde2d(Xk([1,3],1:10000)');
contour(X,Y,density,4,'r','LineWidth',2)
%colormap(jet(3000))
title('Sampling (PF)')
%set(h,'ZData',z-10)
xlabel('X (m)');
ylabel('Y (m)');
%set(gca,'Xdir','reverse')
set(gca,'Ydir','reverse')
xlim(xLim);
ylim(yLim);
view(0,-90);

% model = PositionalObsModelX(config);%PositionalObsModelX(config);%Polar2CartGaussianModelX(config);
% 
% xkm1 = [0.5; 0.3; sqrt(2); sqrt(2)];
% pkm1 = [0 0.5; 0 0.3; sqrt(2) sqrt(2); sqrt(2) sqrt(2)];
% Pkm1 = [0 1 0 1; 0 1 0 1; 0 1 0 1; 1 1 0 1];
% yk = model.obs(1,xkm1);
% Yk = model.sample(1, yk, 500000);
% pk = model.obs_cov();
% Pk = model.eval(1,Yk,xkm1);
% [bandwidth,density,X,Y]=kde2d(Yk');
% figure;
% %contour3(X,Y,density,50);
% surf(X,Y,density);
% hold on;