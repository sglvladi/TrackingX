% Assume a target positioned at x = 1, travelling with speed v = 0.1
state = [1;0.1;2;0;3;0.2;4;0.3];

% Create an instance of a 2D Ornstein-Uhlenbeck model
ou = OrnsteinUhlenbeckModelX('NumDims',4,...
                             'VelocityErrVariance',0.1,...
                             'DampingCoefficient',0.1,...
                             'TimestepDuration',duration(1,0,1));

% View the transition matrix and process covariance matrices
F = ou.feval();
Q = ou.covar();

% Predict the target's position and velocity after the interval has passed
newState  = ou.feval(state);

% Do the same as above, but this time add process noise to the prediction
newState2 = ou.feval(state,true);

% Generate 50 random noise samples from the dynamic model
noise = ou.random(50);

% Check how likely the predictions we made are
lik = ou.pdf(newState,state);
lik2 = ou.pdf(newState2,state); % HINT: newState2 should be less likely