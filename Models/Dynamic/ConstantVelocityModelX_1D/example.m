% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [1;0.1];

% Create an instance of a 1D Constant Velocity model
cv = ConstantVelocityModel_1D();

% Set the velocity variance to 0.05 and assume a 1 sec interval
cv.VelocityErrVariance = 0.05;
cv.TimeVariant = 1;

% View the transition matrix and process covariance matrices
F = cv.feval();
Q = cv.covariance();

% Predict the target's position and velocity after the interval has passed
newState  = cv.feval(state);

% Do the same as above, but this time add process noise to the prediction
newState2 = cv.feval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = cv.random(50);

% Check how likely the predictions we made are
lik = cv.pdf(newState,state);
lik2 = cv.pdf(newState2,state); % HINT: newState2 should be less likely