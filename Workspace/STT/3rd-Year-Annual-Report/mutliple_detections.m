% Assume a target positioned at x = 1, travelling with speed v = 0.1
state1 = [2000;0.1;0;0];    % Assume state with four dimensions [x_pos, x_vel, y_pos, y_vel]
state2 = [1500;0.1;1500;0];

% Create an instance of a 1D Constant Velocity model
obs = Polar2CartGaussModelX('NumStateDims',4,'RangeErrVariance',1^2,...
                          'ThetaErrVariance',(pi/90)^2,'Mapping',[1 3]);

% View the transition matrix and process covariance matrices
Q = obs.covariance();

% Predict the target's position and velocity after the interval has passed
measurement1  = obs.heval(state1);
measurement2 = obs.heval(state2);

% Generate 50 random noise samples from the dynamic model
noise1 = obs.random(5000000);
noise2 = obs.random(5000000);

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
plot3([measurement1(1,:),measurement2(1,:)],[measurement1(2,:),measurement2(2,:)],...
       ones(1,2).*10^5,'rs','MarkerSize',50,'LineWidth',5);
plot3([measurement1(1,:),measurement2(1,:)],[measurement1(2,:),measurement2(2,:)],...
       ones(1,2).*10^5,'r.','MarkerSize',30,'LineWidth',3);
shading interp
title('Radar measurement noise (Polar)')
xlabel('Bearing (rad)');
ylabel('Range (m)');
xlim([-0.2847,1.0669])
ylim([1962.64950677731,2158.10672174215])
view(2)

% Right plot
subplot(1,2,2);
[bandwidth,density,X,Y]=kde2d(Xk([1,3],:)');
%contour3(X,Y,density,50);
surf(X,Y,density);
hold on;
plot3([state1(1,:),state2(1,:)],[state1(3,:),state2(3,:)],...
       ones(1,2).*10^5,'rs','MarkerSize',50,'LineWidth',5);
plot3([state1(1,:),state2(1,:)],[state1(3,:),state2(3,:)],...
       ones(1,2).*10^5,'r.','MarkerSize',30,'LineWidth',3);
shading interp
title('Radar measurement noise (Cartesian)')
xlabel('X (m)');
ylabel('Y (m)');
view(2)
xlim([1173.81974248927,2049.35622317596]);
ylim([-448.684210526315,1817.10526315790]);

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