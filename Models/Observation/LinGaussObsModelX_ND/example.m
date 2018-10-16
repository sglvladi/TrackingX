% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [1;0.1;2;0];    % Assume state with four dimensions [x_pos, x_vel, y_pos, y_vel]
measurement = [1.2;1.8;3];% Measurement state [x_pos,y_pos]

% Create an instance of a 1D Constant Velocity model
obs = LinGaussObsModelX_ND(3,'NumStateDims',4,'ObsErrVariance',50,'Mapping',[1 3 2]);

% View the transition matrix and process covariance matrices
H = obs.heval();
R = obs.covariance();

% Predict the target's position and velocity after the interval has passed
projectedState = obs.heval(state);

% Do the same as above, but this time add process noise to the prediction
projectedState2 = obs.heval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = obs.random(50);

% Check how likely the predictions we made are
lik = obs.pdf(projectedState,state);
lik2 = obs.pdf(projectedState2,state); % HINT: newState2 should be less likely