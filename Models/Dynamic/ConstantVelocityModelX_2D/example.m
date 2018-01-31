% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [1;0.1;2;0];

% Create an instance of a 1D Constant Velocity model with observation noise
% variance equal to 0.1 and specifying a mapping between the 1st index of
% the measurement vector to the 1st index of the state vector and from the
% 2nd index of the measurement vector to the 3rd index of the state vector.
obs = LinGaussObsModel_2D('NumStateDims',4,'ObsErrVariance',0.1,'Mapping',[1 3]);

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