% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [1;2;0.1;pi/3];

% Create an instance of a 1D Constant Velocity model with observation noise
% variance equal to 0.1 and specifying a mapping between the 1st index of
% the measurement vector to the 1st index of the state vector and from the
% 2nd index of the measurement vector to the 3rd index of the state vector.
ch = ConstantHeadingModelX('VelocityErrVariance',0.0001,'HeadingErrVariance',0.09);

% View the transition matrix and process covariance matrices
Q = ch.covariance();

% Predict the target's position and velocity after the interval has passed
newState  = ch.feval(state);

% Do the same as above, but this time add process noise to the prediction
newState2 = ch.feval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = ch.random(50);

% Check how likely the predictions we made are
lik = ch.pdf(newState,state);
lik2 = ch.pdf(newState2,state); % HINT: newState2 should be less likely