% Assume a target positioned at x = 1, travelling with speed v = 0.1
state1 = [1000;0.1;1000;0];    % Assume state with four dimensions [x_pos, x_vel, y_pos, y_vel]
state2 = [1500;0.1;1500;0];

% Create an instance of a 1D Constant Velocity model
obs = Polar2CartGaussModelX('NumStateDims',4,'RangeErrVariance',10^2,...
                          'ThetaErrVariance',(pi/45)^2,'Mapping',[1 3]);

% View the transition matrix and process covariance matrices
Q = obs.covariance();

% Predict the target's position and velocity after the interval has passed
measurement1  = obs.heval(state);
measurement2 = obs.heval(state2);

% Generate 50 random noise samples from the dynamic model
noise1 = obs.random(500000);
noise2 = obs.random(500000);

% Add noise to the measurements
Yk1 = noise1 + measurement1;
Xk1 = obs.heval_inv(Yk1);

Yk2 = noise2 + measurement2;
Xk2 = obs.heval_inv(Yk2);

Xk = [Xk1,Xk2];
Yk = [Yk1,Yk2];

figure;

% Left plot
subplot(1,2,1);
[bandwidth,density,X,Y]=kde2d(Yk([1,2],:)');
surf(X,Y,density);
hold on;
shading interp
colormap(jet(3000))
title('Radar measurement noise (Polar)')
xlabel('Bearing (rad)');
ylabel('Range (m)');

% Right plot
subplot(1,2,2);
[bandwidth,density,X,Y]=kde2d(Xk([1,3],:)');
%contour3(X,Y,density,50);
surf(X,Y,density);
hold on;
shading interp
colormap(jet(3000))
title('Radar measurement noise (Cartesian)')
xlabel('X (m)');
ylabel('Y (m)');

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