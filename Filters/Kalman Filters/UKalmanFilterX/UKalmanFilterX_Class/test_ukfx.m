N = 500; % Simulate 500 seconds/iterations

% Constant Velocity Model
config_cv.dim = 2;
config_cv.q = 0.0001;
CVmodel = ConstantVelocityModel(config_cv);

% Constant Heading Model
config.q_vel = 0.01;
config.q_head = 0.16;
CHmodel = ConstantHeadingModel(config);

% Positional Observation Model
config_meas.dim = 2;
config_meas.r = 0.2;
obs_model = PositionalObsModel(config_meas);

% Set initial true target state
s = zeros(config_cv.dim*2,N);
s(:,1) = [0; 0; 0.1; 0.3]; % (Position target at 0,0 with velocity of 1 m/s on each axis

% Initiate Kalman Filter
config_kf.k = 1;  % Use k as 1 sec Dt interval for CV model
config_kf.x = s(:,1); 
config_kf.P = CHmodel.config.Q(1);
ukf = UKalmanFilterX(config_kf, CHmodel, obs_model);

% Containers
filtered_estimates = cell(1,N);


% START OF SIMULATION
% ===================>
for k = 1:N

    % Generate new state and measurement
    if(k~=1)
        s(:,k) = CHmodel.propagate_parts(s(:,k-1));
    end
    y(:,k) = obs_model.sample_from(obs_model.transform_mean(s(:,k)));

    % Iterate Kalman Filter
    ukf.config.y = y(:,k);
    ukf.Predict(); 
    ukf.Update();

    % Store filtered estimates
    filtered_estimates{k} = ukf.config;
end
% END OF SIMULATION
% ===================>

% Compute smoothed estimates
%smoothed_estimates = ukf.Smooth(filtered_estimates);

% Extract estimates
x_filtered = zeros(config_cv.dim*2,N);
x_smoothed = zeros(config_cv.dim*2,N);
for k = 1:N
    x_filtered(:,k) = filtered_estimates{k}.x;
    %x_smoothed(:,k) = smoothed_estimates{k}.x;
end

figure
plot(s(1,:),s(2,:),'.-k', x_filtered(1,:), x_filtered(2,:), 'o-b', x_smoothed(1,:), x_smoothed(2,:), 'x-r');